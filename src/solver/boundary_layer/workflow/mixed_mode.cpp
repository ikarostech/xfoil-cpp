#include "solver/boundary_layer/workflow/mixed_mode.hpp"

#include <algorithm>
#include <cmath>

#include "model/boundary_layer/physics.hpp"

namespace {
constexpr double kMixedModeConvergenceTolerance = 5.0e-6;
using BoundaryContext = BoundaryLayerMixedModeStationContext;
} // namespace

void BoundaryLayerMixedModeOps::storeStationStateCommon(
    int side, int stationIndex, const BoundaryContext &ctx) {
  if (stationIndex < context_.lattice.get(side).profiles.transitionIndex) {
    context_.lattice.get(side).profiles.skinFrictionCoeff[stationIndex] = ctx.ami;
  } else {
    context_.lattice.get(side).profiles.skinFrictionCoeff[stationIndex] = ctx.cti;
  }
  context_.lattice.get(side).profiles.momentumThickness[stationIndex] = ctx.thi;
  context_.lattice.get(side).profiles.displacementThickness[stationIndex] = ctx.dsi;
  context_.lattice.get(side).profiles.edgeVelocity[stationIndex] = ctx.uei;
  context_.lattice.get(side).profiles.massFlux[stationIndex] = ctx.dsi * ctx.uei;
  context_.lattice.get(side).profiles.skinFrictionCoeffHistory[stationIndex] =
      context_.state.station2.cqz.scalar;

  BoundaryLayerPhysics::refreshCurrentStation(
      context_.state, context_.blCompressibility, context_.blReynolds,
      ctx.xsi, ctx.ami, ctx.cti, ctx.thi, ctx.dsi, ctx.dswaki, ctx.uei);
  context_.state.stepbl();

  if (context_.flowRegime == FlowRegimeEnum::Wake) {
    return;
  }
  if (context_.flowRegime == FlowRegimeEnum::Transition ||
      stationIndex == context_.lattice.get(side).trailingEdgeIndex) {
    context_.flowRegime = FlowRegimeEnum::Turbulent;
  } else {
    context_.flowRegime = FlowRegimeEnum::Laminar;
  }
}

double BoundaryLayerMixedModeOps::fallbackEdgeVelocity(
    int side, int stationIndex,
    BoundaryLayerEdgeVelocityFallbackMode edgeMode) const {
  switch (edgeMode) {
  case BoundaryLayerEdgeVelocityFallbackMode::UsePreviousStation:
    return context_.lattice.get(side).profiles.edgeVelocity[stationIndex - 1];
  case BoundaryLayerEdgeVelocityFallbackMode::AverageNeighbors: {
    double uei = context_.lattice.get(side).profiles.edgeVelocity[stationIndex];
    if (stationIndex < context_.lattice.get(side).stationCount - 1) {
      uei = 0.5 * (context_.lattice.get(side).profiles.edgeVelocity[stationIndex - 1] +
                   context_.lattice.get(side).profiles.edgeVelocity[stationIndex + 1]);
    }
    return uei;
  }
  }
  return context_.lattice.get(side).profiles.edgeVelocity[stationIndex];
}

void BoundaryLayerMixedModeOps::syncStationRegimeStates(
    int side, int stationIndex, FlowRegimeEnum stationRegime) {
  if (stationIndex < context_.lattice.get(side).profiles.transitionIndex) {
    context_.state.station2 = context_.boundaryLayerVariablesSolver.solve(
        context_.state.station2, FlowRegimeEnum::Laminar);
    context_.solverOps.blmid(FlowRegimeEnum::Laminar);
  }
  if (stationIndex >= context_.lattice.get(side).profiles.transitionIndex) {
    context_.state.station2 = context_.boundaryLayerVariablesSolver.solve(
        context_.state.station2, FlowRegimeEnum::Turbulent);
    context_.solverOps.blmid(FlowRegimeEnum::Turbulent);
  }
  if (stationRegime == FlowRegimeEnum::Wake) {
    context_.state.station2 = context_.boundaryLayerVariablesSolver.solve(
        context_.state.station2, FlowRegimeEnum::Wake);
    context_.solverOps.blmid(FlowRegimeEnum::Wake);
  }
  context_.flowRegime = determineRegimeForStation(side, stationIndex);
}

FlowRegimeEnum BoundaryLayerMixedModeOps::determineRegimeForStation(
    int side, int stationIndex) const {
  if (stationIndex > context_.lattice.get(side).trailingEdgeIndex) {
    return FlowRegimeEnum::Wake;
  }
  const int transition_index = context_.lattice.get(side).profiles.transitionIndex;
  if (stationIndex == transition_index) {
    return FlowRegimeEnum::Transition;
  }
  if (stationIndex > transition_index) {
    return FlowRegimeEnum::Turbulent;
  }
  if (stationIndex == 0) {
    return FlowRegimeEnum::Similarity;
  }
  return FlowRegimeEnum::Laminar;
}

void BoundaryLayerMixedModeOps::updateSystemMatricesForStation(
    const Edge &edge, int side, int stationIndex, BoundaryContext &ctx) {
  if (stationIndex == context_.lattice.get(side).trailingEdgeIndex + 1) {
    ctx.tte = context_.lattice.get(1)
                  .profiles.momentumThickness[context_.lattice.top.trailingEdgeIndex] +
              context_.lattice.get(2)
                  .profiles.momentumThickness[context_.lattice.bottom.trailingEdgeIndex];
    ctx.dte = context_.lattice.get(1)
                  .profiles.displacementThickness[context_.lattice.top.trailingEdgeIndex] +
              context_.lattice.get(2)
                  .profiles.displacementThickness[context_.lattice.bottom.trailingEdgeIndex] +
              edge.ante;
    ctx.cte = (context_.lattice.get(1)
                   .profiles.skinFrictionCoeff[context_.lattice.top.trailingEdgeIndex] *
                   context_.lattice.get(1)
                       .profiles.momentumThickness[context_.lattice.top.trailingEdgeIndex] +
               context_.lattice.get(2)
                   .profiles.skinFrictionCoeff[context_.lattice.bottom.trailingEdgeIndex] *
                   context_.lattice.get(2)
                       .profiles.momentumThickness[context_.lattice.bottom.trailingEdgeIndex]) /
              ctx.tte;
    context_.solverOps.tesys(context_.lattice.top.profiles,
                             context_.lattice.bottom.profiles, edge);
  } else {
    context_.solverOps.blsys();
  }
}

void BoundaryLayerMixedModeOps::initializeFirstIterationState(
    int side, int stationIndex, int previousTransition, BoundaryContext &ctx,
    double &ueref, double &hkref) {
  ueref = context_.state.station2.param.uz;
  hkref = context_.state.station2.hkz.scalar;

  if (stationIndex < context_.lattice.get(side).profiles.transitionIndex &&
      stationIndex >= previousTransition) {
    const double uem = context_.lattice.get(side)
                           .profiles.edgeVelocity[std::max(0, stationIndex - 1)];
    const double dsm = context_.lattice.get(side)
                           .profiles.displacementThickness[std::max(0, stationIndex - 1)];
    const double thm = context_.lattice.get(side)
                           .profiles.momentumThickness[std::max(0, stationIndex - 1)];
    const double uem_sq = uem * uem;
    const double msq = uem_sq * context_.blCompressibility.hstinv() /
                       (context_.blCompressibility.gm1bl() *
                        (1.0 - 0.5 * uem_sq * context_.blCompressibility.hstinv()));
    hkref = boundary_layer::hkin(dsm / thm, msq).hk;
  }

  if (stationIndex < previousTransition) {
    if (context_.flowRegime == FlowRegimeEnum::Transition) {
      context_.lattice.get(side).profiles.skinFrictionCoeff[stationIndex] = 0.03;
    } else {
      context_.lattice.get(side).profiles.skinFrictionCoeff[stationIndex] =
          context_.lattice.get(side)
              .profiles.skinFrictionCoeff[std::max(0, stationIndex - 1)];
    }
    ctx.cti = context_.lattice.get(side)
                  .profiles.skinFrictionCoeff[std::max(0, stationIndex - 1)];
    context_.state.station2.param.sz = ctx.cti;
  }
}

void BoundaryLayerMixedModeOps::configureSimilarityRow(double ueref) {
  context_.blc.a2(3, 0) = 0.0;
  context_.blc.a2(3, 1) = 0.0;
  context_.blc.a2(3, 2) = 0.0;
  context_.blc.a2(3, 3) = context_.state.station2.param.uz_uei;
  context_.blc.rhs[3] = ueref - context_.state.station2.param.uz;
}

void BoundaryLayerMixedModeOps::configureViscousRow(
    double hkref, double ueref, double senswt, bool resetSensitivity,
    bool averageSensitivity, double &sens, double &sennew) {
  context_.blc.a2(3, 0) = 0.0;
  context_.blc.a2(3, 1) = context_.state.station2.hkz.t();
  context_.blc.a2(3, 2) = context_.state.station2.hkz.d();
  context_.blc.a2(3, 3) =
      context_.state.station2.hkz.u() * context_.state.station2.param.uz_uei;
  context_.blc.rhs[3] = 1.0;

  const double delta_sen =
      context_.blc.a2.block(0, 0, 4, 4).fullPivLu().solve(context_.blc.rhs)[3];
  sennew = senswt * delta_sen * hkref / ueref;
  if (resetSensitivity) {
    sens = sennew;
  } else if (averageSensitivity) {
    sens = 0.5 * (sens + sennew);
  }

  context_.blc.a2(3, 1) = context_.state.station2.hkz.t() * hkref;
  context_.blc.a2(3, 2) = context_.state.station2.hkz.d() * hkref;
  context_.blc.a2(3, 3) =
      (context_.state.station2.hkz.u() * hkref + sens / ueref) *
      context_.state.station2.param.uz_uei;
  context_.blc.rhs[3] =
      -(hkref * hkref) * (context_.state.station2.hkz.scalar / hkref - 1.0) -
      sens * (context_.state.station2.param.uz / ueref - 1.0);
}

bool BoundaryLayerMixedModeOps::applyMixedModeNewtonStep(
    int side, int stationIndex, double &ami, BoundaryContext &ctx) {
  context_.blc.rhs =
      context_.blc.a2.block(0, 0, 4, 4).fullPivLu().solve(context_.blc.rhs);

  ctx.dmax = std::max(std::fabs(context_.blc.rhs[1] / ctx.thi),
                      std::fabs(context_.blc.rhs[2] / ctx.dsi));
  if (stationIndex >= context_.lattice.get(side).profiles.transitionIndex) {
    ctx.dmax =
        std::max(ctx.dmax, std::fabs(context_.blc.rhs[0] / (10.0 * ctx.cti)));
  }

  double rlx = 1.0;
  if (ctx.dmax > 0.3) {
    rlx = 0.3 / ctx.dmax;
  }

  if (stationIndex < context_.lattice.get(side).profiles.transitionIndex) {
    ami += rlx * context_.blc.rhs[0];
    ctx.ami = ami;
  }
  if (stationIndex >= context_.lattice.get(side).profiles.transitionIndex) {
    ctx.cti += rlx * context_.blc.rhs[0];
  }
  ctx.thi += rlx * context_.blc.rhs[1];
  ctx.dsi += rlx * context_.blc.rhs[2];
  ctx.uei += rlx * context_.blc.rhs[3];

  if (stationIndex >= context_.lattice.get(side).profiles.transitionIndex) {
    ctx.cti = std::clamp(ctx.cti, 0.0000001, 0.30);
  }

  const double hklim =
      (stationIndex <= context_.lattice.get(side).trailingEdgeIndex) ? 1.02
                                                                     : 1.00005;
  const double uei_sq = ctx.uei * ctx.uei;
  const double msq = uei_sq * context_.blCompressibility.hstinv() /
                     (context_.blCompressibility.gm1bl() *
                      (1.0 - 0.5 * uei_sq * context_.blCompressibility.hstinv()));
  double dsw = ctx.dsi - ctx.dswaki;
  dsw = BoundaryLayerPhysics::adjustDisplacementForHkLimit(dsw, ctx.thi, msq,
                                                           hklim);
  ctx.dsi = dsw + ctx.dswaki;

  return ctx.dmax <= kMixedModeConvergenceTolerance;
}

void BoundaryLayerMixedModeOps::checkTransitionIfNeeded(
    int side, int stationIndex, bool skipCheck, int laminarAdvance,
    double &ami) {
  if (skipCheck || context_.flowRegime == FlowRegimeEnum::Turbulent ||
      context_.flowRegime == FlowRegimeEnum::Wake) {
    return;
  }

  context_.transitionSolver.trchek();
  ami = context_.state.station2.param.amplz;
  if (context_.flowRegime == FlowRegimeEnum::Transition) {
    context_.lattice.get(side).profiles.transitionIndex = stationIndex;
  } else {
    context_.lattice.get(side).profiles.transitionIndex = stationIndex + laminarAdvance;
  }
}

void BoundaryLayerMixedModeOps::resetStationKinematicsAfterFailure(
    int side, int stationIndex, BoundaryContext &ctx,
    BoundaryLayerEdgeVelocityFallbackMode edgeMode) {
  if (ctx.dmax <= 0.1 || stationIndex < 2) {
    return;
  }

  if (stationIndex <= context_.lattice.get(side).trailingEdgeIndex) {
    const double ratio = context_.lattice.get(side).arcLengthCoordinates[stationIndex] /
                         context_.lattice.get(side).arcLengthCoordinates[stationIndex - 1];
    const double scale = std::sqrt(ratio);
    ctx.thi = context_.lattice.get(side).profiles.momentumThickness[stationIndex - 1] * scale;
    ctx.dsi = context_.lattice.get(side).profiles.displacementThickness[stationIndex - 1] * scale;
  } else if (stationIndex == context_.lattice.get(side).trailingEdgeIndex + 1) {
    ctx.cti = ctx.cte;
    ctx.thi = ctx.tte;
    ctx.dsi = ctx.dte;
  } else {
    ctx.thi = context_.lattice.get(side).profiles.momentumThickness[stationIndex - 1];
    const double ratlen =
        (context_.lattice.get(side).arcLengthCoordinates[stationIndex] -
         context_.lattice.get(side).arcLengthCoordinates[stationIndex - 1]) /
        (10.0 * context_.lattice.get(side).profiles.displacementThickness[stationIndex - 1]);
    ctx.dsi =
        (context_.lattice.get(side).profiles.displacementThickness[stationIndex - 1] +
         ctx.thi * ratlen) /
        (1.0 + ratlen);
  }

  ctx.uei = fallbackEdgeVelocity(side, stationIndex, edgeMode);
  if (stationIndex == context_.lattice.get(side).profiles.transitionIndex) {
    ctx.cti = 0.05;
  }
  if (stationIndex > context_.lattice.get(side).profiles.transitionIndex) {
    ctx.cti =
        context_.lattice.get(side).profiles.skinFrictionCoeff[stationIndex - 1];
  }
}

void BoundaryLayerMixedModeOps::recoverStationAfterFailure(
    int side, int stationIndex, BoundaryContext &ctx, double &ami,
    BoundaryLayerEdgeVelocityFallbackMode edgeMode,
    int laminarAdvance) {
  context_.flowRegime = ctx.flowRegime;
  resetStationKinematicsAfterFailure(side, stationIndex, ctx, edgeMode);

  BoundaryLayerPhysics::refreshCurrentStation(
      context_.state, context_.blCompressibility, context_.blReynolds,
      ctx.xsi, ctx.ami, ctx.cti, ctx.thi, ctx.dsi, ctx.dswaki, ctx.uei);
  ctx.flowRegime = context_.flowRegime;

  checkTransitionIfNeeded(side, stationIndex, ctx.isSimilarity(), laminarAdvance,
                          ami);
  ctx.flowRegime = context_.flowRegime;

  syncStationRegimeStates(side, stationIndex, ctx.flowRegime);
  ctx.flowRegime = context_.flowRegime;
  ctx.ami = ami;
}
