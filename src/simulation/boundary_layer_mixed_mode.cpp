#include "simulation/boundary_layer_mixed_mode.hpp"

#include <algorithm>
#include <cmath>

#include "domain/boundary_layer.hpp"

namespace {
constexpr double kMixedModeConvergenceTolerance = 5.0e-6;
using BoundaryContext = BoundaryLayerWorkflow::MixedModeStationContext;

double adjustDisplacementForHkLimit(double displacementThickness,
                                    double momentumThickness, double msq,
                                    double hklim) {
  const double h = displacementThickness / momentumThickness;
  const auto hkin_result = boundary_layer::hkin(h, msq);
  const double dh = std::max(0.0, hklim - hkin_result.hk) / hkin_result.hk_h;
  return displacementThickness + dh * momentumThickness;
}
} // namespace

void BoundaryLayerMixedModeOps::storeStationStateCommon(
    int side, int stationIndex, const BoundaryContext &ctx) {
  if (stationIndex < workflow_.lattice.get(side).profiles.transitionIndex) {
    workflow_.lattice.get(side).profiles.skinFrictionCoeff[stationIndex] =
        ctx.ami;
  } else {
    workflow_.lattice.get(side).profiles.skinFrictionCoeff[stationIndex] =
        ctx.cti;
  }
  workflow_.lattice.get(side).profiles.momentumThickness[stationIndex] =
      ctx.thi;
  workflow_.lattice.get(side).profiles.displacementThickness[stationIndex] =
      ctx.dsi;
  workflow_.lattice.get(side).profiles.edgeVelocity[stationIndex] = ctx.uei;
  workflow_.lattice.get(side).profiles.massFlux[stationIndex] = ctx.dsi * ctx.uei;
  workflow_.lattice.get(side).profiles.skinFrictionCoeffHistory[stationIndex] =
      workflow_.state.station2.cqz.scalar;

  workflow_.state.current() = workflow_.blprv(
      workflow_.state.current(), ctx.xsi, ctx.ami, ctx.cti, ctx.thi, ctx.dsi,
      ctx.dswaki, ctx.uei);
  workflow_.blkin(workflow_.state);
  workflow_.state.stepbl();

  if (workflow_.flowRegime == FlowRegimeEnum::Wake) {
    return;
  }
  if (workflow_.flowRegime == FlowRegimeEnum::Transition ||
      stationIndex == workflow_.lattice.get(side).trailingEdgeIndex) {
    workflow_.flowRegime = FlowRegimeEnum::Turbulent;
  } else {
    workflow_.flowRegime = FlowRegimeEnum::Laminar;
  }
}

double BoundaryLayerMixedModeOps::fallbackEdgeVelocity(
    int side, int stationIndex,
    BoundaryLayerWorkflow::EdgeVelocityFallbackMode edgeMode) const {
  switch (edgeMode) {
  case BoundaryLayerWorkflow::EdgeVelocityFallbackMode::UsePreviousStation:
    return workflow_.lattice.get(side).profiles.edgeVelocity[stationIndex - 1];
  case BoundaryLayerWorkflow::EdgeVelocityFallbackMode::AverageNeighbors: {
    double uei = workflow_.lattice.get(side).profiles.edgeVelocity[stationIndex];
    if (stationIndex < workflow_.lattice.get(side).stationCount - 1) {
      uei = 0.5 * (workflow_.lattice.get(side).profiles.edgeVelocity[stationIndex - 1] +
                   workflow_.lattice.get(side).profiles.edgeVelocity[stationIndex + 1]);
    }
    return uei;
  }
  }
  return workflow_.lattice.get(side).profiles.edgeVelocity[stationIndex];
}

void BoundaryLayerMixedModeOps::syncStationRegimeStates(
    int side, int stationIndex, FlowRegimeEnum stationRegime) {
  if (stationIndex < workflow_.lattice.get(side).profiles.transitionIndex) {
    workflow_.state.station2 = workflow_.boundaryLayerVariablesSolver.solve(
        workflow_.state.station2, FlowRegimeEnum::Laminar);
    workflow_.blmid(FlowRegimeEnum::Laminar);
  }
  if (stationIndex >= workflow_.lattice.get(side).profiles.transitionIndex) {
    workflow_.state.station2 = workflow_.boundaryLayerVariablesSolver.solve(
        workflow_.state.station2, FlowRegimeEnum::Turbulent);
    workflow_.blmid(FlowRegimeEnum::Turbulent);
  }
  if (stationRegime == FlowRegimeEnum::Wake) {
    workflow_.state.station2 = workflow_.boundaryLayerVariablesSolver.solve(
        workflow_.state.station2, FlowRegimeEnum::Wake);
    workflow_.blmid(FlowRegimeEnum::Wake);
  }
  workflow_.flowRegime = determineRegimeForStation(side, stationIndex);
}

FlowRegimeEnum BoundaryLayerMixedModeOps::determineRegimeForStation(
    int side, int stationIndex) const {
  if (stationIndex > workflow_.lattice.get(side).trailingEdgeIndex) {
    return FlowRegimeEnum::Wake;
  }
  const int transition_index = workflow_.lattice.get(side).profiles.transitionIndex;
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
  if (workflow_.isStartOfWake(side, stationIndex)) {
    ctx.tte = workflow_.lattice.get(1)
                  .profiles.momentumThickness[workflow_.lattice.top.trailingEdgeIndex] +
              workflow_.lattice.get(2)
                  .profiles.momentumThickness[workflow_.lattice.bottom.trailingEdgeIndex];
    ctx.dte = workflow_.lattice.get(1)
                  .profiles.displacementThickness[workflow_.lattice.top.trailingEdgeIndex] +
              workflow_.lattice.get(2)
                  .profiles.displacementThickness[workflow_.lattice.bottom.trailingEdgeIndex] +
              edge.ante;
    ctx.cte = (workflow_.lattice.get(1)
                   .profiles.skinFrictionCoeff[workflow_.lattice.top.trailingEdgeIndex] *
                   workflow_.lattice.get(1)
                       .profiles.momentumThickness[workflow_.lattice.top.trailingEdgeIndex] +
               workflow_.lattice.get(2)
                   .profiles.skinFrictionCoeff[workflow_.lattice.bottom.trailingEdgeIndex] *
                   workflow_.lattice.get(2)
                       .profiles.momentumThickness[workflow_.lattice.bottom.trailingEdgeIndex]) /
              ctx.tte;
    workflow_.tesys(workflow_.lattice.top.profiles,
                    workflow_.lattice.bottom.profiles, edge);
  } else {
    workflow_.blsys();
  }
}

void BoundaryLayerMixedModeOps::initializeFirstIterationState(
    int side, int stationIndex, int previousTransition, BoundaryContext &ctx,
    double &ueref, double &hkref) {
  ueref = workflow_.state.station2.param.uz;
  hkref = workflow_.state.station2.hkz.scalar;

  if (stationIndex < workflow_.lattice.get(side).profiles.transitionIndex &&
      stationIndex >= previousTransition) {
    const double uem = workflow_.lattice.get(side)
                           .profiles.edgeVelocity[std::max(0, stationIndex - 1)];
    const double dsm = workflow_.lattice.get(side)
                           .profiles.displacementThickness[std::max(0, stationIndex - 1)];
    const double thm = workflow_.lattice.get(side)
                           .profiles.momentumThickness[std::max(0, stationIndex - 1)];
    const double uem_sq = uem * uem;
    const double msq = uem_sq * workflow_.blCompressibility.hstinv /
                       (workflow_.blCompressibility.gm1bl *
                        (1.0 - 0.5 * uem_sq * workflow_.blCompressibility.hstinv));
    hkref = boundary_layer::hkin(dsm / thm, msq).hk;
  }

  if (stationIndex < previousTransition) {
    if (workflow_.flowRegime == FlowRegimeEnum::Transition) {
      workflow_.lattice.get(side).profiles.skinFrictionCoeff[stationIndex] = 0.03;
    } else {
      workflow_.lattice.get(side).profiles.skinFrictionCoeff[stationIndex] =
          workflow_.lattice.get(side)
              .profiles.skinFrictionCoeff[std::max(0, stationIndex - 1)];
    }
    ctx.cti = workflow_.lattice.get(side)
                  .profiles.skinFrictionCoeff[std::max(0, stationIndex - 1)];
    workflow_.state.station2.param.sz = ctx.cti;
  }
}

void BoundaryLayerMixedModeOps::configureSimilarityRow(double ueref) {
  workflow_.blc.a2(3, 0) = 0.0;
  workflow_.blc.a2(3, 1) = 0.0;
  workflow_.blc.a2(3, 2) = 0.0;
  workflow_.blc.a2(3, 3) = workflow_.state.station2.param.uz_uei;
  workflow_.blc.rhs[3] = ueref - workflow_.state.station2.param.uz;
}

void BoundaryLayerMixedModeOps::configureViscousRow(
    double hkref, double ueref, double senswt, bool resetSensitivity,
    bool averageSensitivity, double &sens, double &sennew) {
  workflow_.blc.a2(3, 0) = 0.0;
  workflow_.blc.a2(3, 1) = workflow_.state.station2.hkz.t();
  workflow_.blc.a2(3, 2) = workflow_.state.station2.hkz.d();
  workflow_.blc.a2(3, 3) =
      workflow_.state.station2.hkz.u() * workflow_.state.station2.param.uz_uei;
  workflow_.blc.rhs[3] = 1.0;

  const double delta_sen =
      workflow_.blc.a2.block(0, 0, 4, 4).fullPivLu().solve(workflow_.blc.rhs)[3];
  sennew = senswt * delta_sen * hkref / ueref;
  if (resetSensitivity) {
    sens = sennew;
  } else if (averageSensitivity) {
    sens = 0.5 * (sens + sennew);
  }

  workflow_.blc.a2(3, 1) = workflow_.state.station2.hkz.t() * hkref;
  workflow_.blc.a2(3, 2) = workflow_.state.station2.hkz.d() * hkref;
  workflow_.blc.a2(3, 3) =
      (workflow_.state.station2.hkz.u() * hkref + sens / ueref) *
      workflow_.state.station2.param.uz_uei;
  workflow_.blc.rhs[3] =
      -(hkref * hkref) * (workflow_.state.station2.hkz.scalar / hkref - 1.0) -
      sens * (workflow_.state.station2.param.uz / ueref - 1.0);
}

bool BoundaryLayerMixedModeOps::applyMixedModeNewtonStep(
    int side, int stationIndex, double &ami, BoundaryContext &ctx) {
  workflow_.blc.rhs =
      workflow_.blc.a2.block(0, 0, 4, 4).fullPivLu().solve(workflow_.blc.rhs);

  ctx.dmax = std::max(std::fabs(workflow_.blc.rhs[1] / ctx.thi),
                      std::fabs(workflow_.blc.rhs[2] / ctx.dsi));
  if (stationIndex >= workflow_.lattice.get(side).profiles.transitionIndex) {
    ctx.dmax =
        std::max(ctx.dmax, std::fabs(workflow_.blc.rhs[0] / (10.0 * ctx.cti)));
  }

  double rlx = 1.0;
  if (ctx.dmax > 0.3) {
    rlx = 0.3 / ctx.dmax;
  }

  if (stationIndex < workflow_.lattice.get(side).profiles.transitionIndex) {
    ami += rlx * workflow_.blc.rhs[0];
    ctx.ami = ami;
  }
  if (stationIndex >= workflow_.lattice.get(side).profiles.transitionIndex) {
    ctx.cti += rlx * workflow_.blc.rhs[0];
  }
  ctx.thi += rlx * workflow_.blc.rhs[1];
  ctx.dsi += rlx * workflow_.blc.rhs[2];
  ctx.uei += rlx * workflow_.blc.rhs[3];

  if (stationIndex >= workflow_.lattice.get(side).profiles.transitionIndex) {
    ctx.cti = std::clamp(ctx.cti, 0.0000001, 0.30);
  }

  const double hklim =
      (stationIndex <= workflow_.lattice.get(side).trailingEdgeIndex) ? 1.02 : 1.00005;
  const double uei_sq = ctx.uei * ctx.uei;
  const double msq = uei_sq * workflow_.blCompressibility.hstinv /
                     (workflow_.blCompressibility.gm1bl *
                      (1.0 - 0.5 * uei_sq * workflow_.blCompressibility.hstinv));
  double dsw = ctx.dsi - ctx.dswaki;
  dsw = adjustDisplacementForHkLimit(dsw, ctx.thi, msq, hklim);
  ctx.dsi = dsw + ctx.dswaki;

  return ctx.dmax <= kMixedModeConvergenceTolerance;
}

void BoundaryLayerMixedModeOps::checkTransitionIfNeeded(
    int side, int stationIndex, bool skipCheck, int laminarAdvance,
    double &ami) {
  if (skipCheck || workflow_.flowRegime == FlowRegimeEnum::Turbulent ||
      workflow_.flowRegime == FlowRegimeEnum::Wake) {
    return;
  }

  workflow_.transitionSolver.trchek();
  ami = workflow_.state.station2.param.amplz;
  if (workflow_.flowRegime == FlowRegimeEnum::Transition) {
    workflow_.lattice.get(side).profiles.transitionIndex = stationIndex;
  } else {
    workflow_.lattice.get(side).profiles.transitionIndex =
        stationIndex + laminarAdvance;
  }
}

void BoundaryLayerMixedModeOps::resetStationKinematicsAfterFailure(
    int side, int stationIndex, BoundaryContext &ctx,
    BoundaryLayerWorkflow::EdgeVelocityFallbackMode edgeMode) {
  if (ctx.dmax <= 0.1 || stationIndex < 2) {
    return;
  }

  if (stationIndex <= workflow_.lattice.get(side).trailingEdgeIndex) {
    const double ratio = workflow_.lattice.get(side).arcLengthCoordinates[stationIndex] /
                         workflow_.lattice.get(side).arcLengthCoordinates[stationIndex - 1];
    const double scale = std::sqrt(ratio);
    ctx.thi = workflow_.lattice.get(side).profiles.momentumThickness[stationIndex - 1] * scale;
    ctx.dsi = workflow_.lattice.get(side).profiles.displacementThickness[stationIndex - 1] * scale;
  } else if (stationIndex == workflow_.lattice.get(side).trailingEdgeIndex + 1) {
    ctx.cti = ctx.cte;
    ctx.thi = ctx.tte;
    ctx.dsi = ctx.dte;
  } else {
    ctx.thi = workflow_.lattice.get(side).profiles.momentumThickness[stationIndex - 1];
    const double ratlen =
        (workflow_.lattice.get(side).arcLengthCoordinates[stationIndex] -
         workflow_.lattice.get(side).arcLengthCoordinates[stationIndex - 1]) /
        (10.0 * workflow_.lattice.get(side).profiles.displacementThickness[stationIndex - 1]);
    ctx.dsi =
        (workflow_.lattice.get(side).profiles.displacementThickness[stationIndex - 1] +
         ctx.thi * ratlen) /
        (1.0 + ratlen);
  }

  ctx.uei = fallbackEdgeVelocity(side, stationIndex, edgeMode);
  if (stationIndex == workflow_.lattice.get(side).profiles.transitionIndex) {
    ctx.cti = 0.05;
  }
  if (stationIndex > workflow_.lattice.get(side).profiles.transitionIndex) {
    ctx.cti =
        workflow_.lattice.get(side).profiles.skinFrictionCoeff[stationIndex - 1];
  }
}

void BoundaryLayerMixedModeOps::recoverStationAfterFailure(
    int side, int stationIndex, BoundaryContext &ctx, double &ami,
    BoundaryLayerWorkflow::EdgeVelocityFallbackMode edgeMode,
    int laminarAdvance) {
  ctx.flowRegime = workflow_.applyFlowRegimeCandidate(ctx.flowRegime);
  resetStationKinematicsAfterFailure(side, stationIndex, ctx, edgeMode);

  workflow_.state.current() = workflow_.blprv(
      workflow_.state.current(), ctx.xsi, ctx.ami, ctx.cti, ctx.thi, ctx.dsi,
      ctx.dswaki, ctx.uei);
  workflow_.blkin(workflow_.state);
  ctx.flowRegime = workflow_.currentFlowRegime();

  checkTransitionIfNeeded(side, stationIndex, ctx.isSimilarity(), laminarAdvance,
                          ami);
  ctx.flowRegime = workflow_.currentFlowRegime();

  syncStationRegimeStates(side, stationIndex, ctx.flowRegime);
  ctx.flowRegime = workflow_.currentFlowRegime();
  ctx.ami = ami;
}
