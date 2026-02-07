#include "BoundaryLayer_march.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "XFoil.h"
#include "infrastructure/logger.hpp"

using BoundaryContext = BoundaryLayerWorkflow::MixedModeStationContext;
using EdgeVelocityFallbackMode = BoundaryLayerWorkflow::EdgeVelocityFallbackMode;
using EdgeVelocityDistribution = BoundaryLayerWorkflow::EdgeVelocityDistribution;
using QtanResult = BoundaryLayerWorkflow::QtanResult;
using ClContributions = BoundaryLayerWorkflow::ClContributions;
using BoundaryLayerDelta = BoundaryLayerWorkflow::BoundaryLayerDelta;
using BoundaryLayerMetrics = BoundaryLayerWorkflow::BoundaryLayerMetrics;
using MrchueStationContext = BoundaryLayerMarcher::MrchueStationContext;

int BoundaryLayerMarcher::resetSideState(BoundaryLayerWorkflow& workflow, int side, XFoil& xfoil) {
  const int previousTransition = workflow.lattice.get(side).profiles.transitionIndex;
  workflow.blTransition.xiforc = workflow.xifset(xfoil.foil, xfoil.stagnation, side);
  workflow.flowRegime = FlowRegimeEnum::Laminar;
  workflow.lattice.get(side).profiles.transitionIndex = workflow.lattice.get(side).trailingEdgeIndex;
  return previousTransition;
}

void BoundaryLayerMarcher::storeStationStateCommon(
    BoundaryLayerWorkflow& workflow, int side, int stationIndex, double ami,
    double cti, double thi, double dsi, double uei, double xsi,
    double dswaki) {
  if (stationIndex < workflow.lattice.get(side).profiles.transitionIndex) {
    workflow.lattice.get(side).profiles.skinFrictionCoeff[stationIndex] = ami;
  } else {
    workflow.lattice.get(side).profiles.skinFrictionCoeff[stationIndex] = cti;
  }
  workflow.lattice.get(side).profiles.momentumThickness[stationIndex] = thi;
  workflow.lattice.get(side).profiles.displacementThickness[stationIndex] = dsi;
  workflow.lattice.get(side).profiles.edgeVelocity[stationIndex] = uei;
  workflow.lattice.get(side).profiles.massFlux[stationIndex] = dsi * uei;
  workflow.lattice.get(side).profiles.skinFrictionCoeffHistory[stationIndex] =
      workflow.state.station2.cqz.scalar;

  {
    blData updatedCurrent =
        workflow.blprv(workflow.state.current(), xsi, ami, cti, thi, dsi,
                        dswaki, uei);
    workflow.state.current() = updatedCurrent;
  }
  workflow.blkin(workflow.state);
  workflow.state.stepbl();

  if (workflow.flowRegime == FlowRegimeEnum::Wake) {
    // Keep wake state.
  } else if (workflow.flowRegime == FlowRegimeEnum::Transition ||
             stationIndex == workflow.lattice.get(side).trailingEdgeIndex) {
    workflow.flowRegime = FlowRegimeEnum::Turbulent;
  } else {
    workflow.flowRegime = FlowRegimeEnum::Laminar;
  }
}

BoundaryLayerWorkflow::MixedModeStationContext
BoundaryLayerMarcher::prepareMixedModeStation(BoundaryLayerWorkflow& workflow, int side, int stationIndex,
                                               int previousTransition,
                                               double& ami) {
  BoundaryContext ctx;

  ctx.simi = (stationIndex == 0);
  ctx.wake = stationIndex > workflow.lattice.get(side).trailingEdgeIndex;
  ctx.xsi = workflow.lattice.get(side).arcLengthCoordinates[stationIndex];
  ctx.uei = workflow.lattice.get(side).profiles.edgeVelocity[stationIndex];
  ctx.thi = workflow.lattice.get(side).profiles.momentumThickness[stationIndex];
  ctx.dsi = workflow.lattice.get(side).profiles.displacementThickness[stationIndex];

  if (stationIndex < previousTransition) {
    ami = workflow.lattice.get(side).profiles.skinFrictionCoeff[stationIndex];
    ctx.cti = 0.03;
  } else {
    ctx.cti = workflow.lattice.get(side).profiles.skinFrictionCoeff[stationIndex];
    if (ctx.cti <= 0.0) {
      ctx.cti = 0.03;
    }
  }
  ctx.ami = ami;

  if (ctx.wake) {
    int iw = stationIndex - workflow.lattice.get(side).trailingEdgeIndex;
    ctx.dswaki = workflow.wgap[iw - 1];
  } else {
    ctx.dswaki = 0.0;
  }

  double thickness_limit =
      (stationIndex <= workflow.lattice.get(side).trailingEdgeIndex) ? 1.02 : 1.00005;
  ctx.dsi = std::max(ctx.dsi - ctx.dswaki, thickness_limit * ctx.thi) +
            ctx.dswaki;

  workflow.flowRegime =
      determineRegimeForStation(workflow, side, stationIndex, ctx.simi, ctx.wake);

  return ctx;
}

double BoundaryLayerMarcher::fallbackEdgeVelocity(
    const BoundaryLayerWorkflow& workflow, int side, int stationIndex,
    EdgeVelocityFallbackMode edgeMode) const {
  switch (edgeMode) {
    case EdgeVelocityFallbackMode::UsePreviousStation:
      return workflow.lattice.get(side).profiles.edgeVelocity[stationIndex - 1];
    case EdgeVelocityFallbackMode::AverageNeighbors: {
      double uei = workflow.lattice.get(side).profiles.edgeVelocity[stationIndex];
      if (stationIndex < workflow.lattice.get(side).stationCount - 1) {
        uei = 0.5 * (workflow.lattice.get(side).profiles.edgeVelocity[stationIndex - 1] +
                     workflow.lattice.get(side).profiles.edgeVelocity[stationIndex + 1]);
      }
      return uei;
    }
  }
  return workflow.lattice.get(side).profiles.edgeVelocity[stationIndex];
}

namespace {

void applyRelaxationLimit(const Eigen::VectorXd& dn, double dhi, double dlo,
                          double& relaxation) {
  double max_pos = 0.0;
  double min_neg = 0.0;
  for (const double value : dn) {
    if (value > max_pos) {
      max_pos = value;
    }
    if (value < min_neg) {
      min_neg = value;
    }
  }
  if (max_pos > 0.0) {
    relaxation = std::min(relaxation, dhi / max_pos);
  }
  if (min_neg < 0.0) {
    relaxation = std::min(relaxation, dlo / min_neg);
  }
}

double adjustDisplacementForHkLimit(double displacementThickness,
                                    double momentumThickness, double msq,
                                    double hklim) {
  const double h = displacementThickness / momentumThickness;
  const auto hkin_result = boundary_layer::hkin(h, msq);
  const double dh =
      std::max(0.0, hklim - hkin_result.hk) / hkin_result.hk_h;
  return displacementThickness + dh * momentumThickness;
}

}  // namespace

BoundaryLayerWorkflow::EdgeVelocityDistribution
BoundaryLayerMarcher::computeNewUeDistribution(const BoundaryLayerWorkflow& workflow, const XFoil& xfoil) const {
  EdgeVelocityDistribution distribution;
  distribution.unew.top =
      Eigen::VectorXd::Zero(workflow.lattice.top.stationCount);
  distribution.unew.bottom =
      Eigen::VectorXd::Zero(workflow.lattice.bottom.stationCount);
  distribution.u_ac.top =
      Eigen::VectorXd::Zero(workflow.lattice.top.stationCount);
  distribution.u_ac.bottom =
      Eigen::VectorXd::Zero(workflow.lattice.bottom.stationCount);

  for (int side = 1; side <= 2; ++side) {
    for (int station = 0; station < workflow.lattice.get(side).stationCount - 1; ++station) {
      const int panelIndex = workflow.lattice.get(side).stationToPanel[station];
      double dui = 0.0;
      double dui_ac = 0.0;
      for (int otherSide = 1; otherSide <= 2; ++otherSide) {
        for (int otherStation = 0; otherStation < workflow.lattice.get(otherSide).stationCount - 1; ++otherStation) {
          const int otherPanel = workflow.lattice.get(otherSide).stationToPanel[otherStation];
          const int systemIndex = workflow.lattice.get(otherSide).stationToSystem[otherStation];
          const double influence =
              -workflow.lattice.get(side).panelInfluenceFactor[station] *
              workflow.lattice.get(otherSide).panelInfluenceFactor[otherStation] *
              xfoil.aerodynamicCache.dij(panelIndex, otherPanel);
          dui += influence *
                 (workflow.lattice.get(otherSide).profiles.massFlux[otherStation] +
                  xfoil.vdel[systemIndex](2, 0));
          dui_ac += influence * (-xfoil.vdel[systemIndex](2, 1));
        }
      }

      const double inviscidDerivative =
          xfoil.analysis_state_.controlByAlpha
              ? 0.0
              : workflow.lattice.get(side).inviscidEdgeVelocityMatrix(1, station);
      distribution.unew.get(side)[station] =
          workflow.lattice.get(side).inviscidEdgeVelocityMatrix(0, station) + dui;
      distribution.u_ac.get(side)[station] =
          inviscidDerivative + dui_ac;
    }
  }
  return distribution;
}

BoundaryLayerWorkflow::QtanResult BoundaryLayerMarcher::computeQtan(
    const BoundaryLayerWorkflow& workflow, const EdgeVelocityDistribution& distribution, int point_count) const {
  QtanResult result;
  result.qnew = Eigen::VectorXd::Zero(point_count);
  result.q_ac = Eigen::VectorXd::Zero(point_count);
  for (int side = 1; side <= 2; ++side) {
    const Eigen::VectorXd& unew_vec = distribution.unew.get(side);
    const Eigen::VectorXd& uac_vec = distribution.u_ac.get(side);
    const int limit = workflow.lattice.get(side).trailingEdgeIndex;
    for (int station = 0; station < limit; ++station) {
      const int panelIndex = workflow.lattice.get(side).stationToPanel[station];
      result.qnew[panelIndex] =
          workflow.lattice.get(side).panelInfluenceFactor[station] * unew_vec[station];
      result.q_ac[panelIndex] =
          workflow.lattice.get(side).panelInfluenceFactor[station] * uac_vec[station];
    }
  }
  return result;
}

BoundaryLayerWorkflow::ClContributions
BoundaryLayerMarcher::computeClFromEdgeVelocityDistribution(
    const BoundaryLayerWorkflow& workflow, const XFoil& xfoil,
    const EdgeVelocityDistribution& distribution) const {
  ClContributions contributions;
  const int point_count = xfoil.foil.foil_shape.n;
  const QtanResult qtan = computeQtan(workflow, distribution, point_count);

  const Eigen::VectorXd& qnew = qtan.qnew;
  const Eigen::VectorXd& q_ac = qtan.q_ac;
  const auto compressibility = xfoil.buildCompressibilityParams();
  if (point_count == 0) {
    return contributions;
  }

  const auto cp_first =
      xfoil.computePressureCoefficient(qnew[0], q_ac[0], compressibility);

  double cpg1 = cp_first.cp;
  double cpg1_ms = cp_first.cp_msq;
  double cpg1_ac = cp_first.cp_velocity_derivative;

  for (int i = 0; i < point_count; i++) {
    const int ip = (i + 1) % point_count;
    const auto cp_next =
        xfoil.computePressureCoefficient(qnew[ip], q_ac[ip], compressibility);

    const double cpg2 = cp_next.cp;
    const double cpg2_ms = cp_next.cp_msq;
    const double cpg2_ac = cp_next.cp_velocity_derivative;

    const Eigen::Vector2d dpoint =
        MathUtil::getRotateMatrix(xfoil.analysis_state_.alpha) * (xfoil.foil.foil_shape.points.col(ip) -
                        xfoil.foil.foil_shape.points.col(i));

    const double ag = 0.5 * (cpg2 + cpg1);
    const double ag_ms = 0.5 * (cpg2_ms + cpg1_ms);
    const double ag_ac = 0.5 * (cpg2_ac + cpg1_ac);

    contributions.cl += dpoint.x() * ag;
    contributions.cl_a += dpoint.y() * ag;
    contributions.cl_ms += dpoint.x() * ag_ms;
    contributions.cl_ac += dpoint.x() * ag_ac;

    cpg1 = cpg2;
    cpg1_ms = cpg2_ms;
    cpg1_ac = cpg2_ac;
  }

  return contributions;
}

BoundaryLayerWorkflow::BoundaryLayerDelta
BoundaryLayerMarcher::buildBoundaryLayerDelta(
    const BoundaryLayerWorkflow& workflow, int side, const Eigen::VectorXd& unew_side,
    const Eigen::VectorXd& u_ac_side, double dac,
    const XFoil& xfoil) const {
  BoundaryLayerDelta delta;
  const int len = workflow.lattice.get(side).stationCount - 1;
  if (len <= 0) {
    return delta;
  }

  delta.dskinFrictionCoeff = Eigen::VectorXd(len);
  delta.dmomentumThickness = Eigen::VectorXd(len);
  delta.ddisplacementThickness = Eigen::VectorXd(len);
  delta.dedgeVelocity = Eigen::VectorXd(len);

  const auto iv = workflow.lattice.get(side).stationToSystem.segment(0, len);
  Eigen::VectorXd dmass(len);
  for (int j = 0; j < len; ++j) {
    const int idx = iv[j];
    delta.dskinFrictionCoeff[j] =
        xfoil.vdel[idx](0, 0) - dac * xfoil.vdel[idx](0, 1);
    delta.dmomentumThickness[j] =
        xfoil.vdel[idx](1, 0) - dac * xfoil.vdel[idx](1, 1);
    dmass[j] = xfoil.vdel[idx](2, 0) - dac * xfoil.vdel[idx](2, 1);
  }

  const Eigen::VectorXd edgeVelocity_segment =
      workflow.lattice.get(side).profiles.edgeVelocity.head(len);
  const Eigen::VectorXd displacementThickness_segment =
      workflow.lattice.get(side).profiles.displacementThickness.head(len);
  const Eigen::VectorXd unew_segment = unew_side.head(len);
  const Eigen::VectorXd uac_segment = u_ac_side.head(len);

  delta.dedgeVelocity =
      unew_segment + dac * uac_segment - edgeVelocity_segment;
  delta.ddisplacementThickness =
      (dmass - displacementThickness_segment.cwiseProduct(delta.dedgeVelocity))
          .cwiseQuotient(edgeVelocity_segment);

  return delta;
}

BoundaryLayerWorkflow::BoundaryLayerMetrics
BoundaryLayerMarcher::evaluateSegmentRelaxation(
    const BoundaryLayerWorkflow& workflow, int side, const BoundaryLayerDelta& delta, double dhi, double dlo,
    double& relaxation) const {
  BoundaryLayerMetrics metrics;
  const int len = delta.dskinFrictionCoeff.size();
  if (len <= 0) {
    return metrics;
  }

  const Eigen::VectorXd skinFrictionCoeff_segment =
      workflow.lattice.get(side).profiles.skinFrictionCoeff.head(len);
  const Eigen::VectorXd momentumThickness_segment =
      workflow.lattice.get(side).profiles.momentumThickness.head(len);
  const Eigen::VectorXd displacementThickness_segment =
      workflow.lattice.get(side).profiles.displacementThickness.head(len);

  Eigen::VectorXd dn1(len);
  const int transition_index = workflow.lattice.get(side).profiles.transitionIndex;
  for (int idx = 0; idx < len; ++idx) {
    dn1[idx] = (idx < transition_index)
                   ? delta.dskinFrictionCoeff[idx] / 10.0
                   : delta.dskinFrictionCoeff[idx] /
                         skinFrictionCoeff_segment[idx];
  }
  const Eigen::VectorXd dn2 =
      delta.dmomentumThickness.cwiseQuotient(momentumThickness_segment);
  const Eigen::VectorXd dn3 =
      delta.ddisplacementThickness.cwiseQuotient(displacementThickness_segment);
  const Eigen::VectorXd dn4 = delta.dedgeVelocity.array().abs() / 0.25;

  applyRelaxationLimit(dn1, dhi, dlo, relaxation);
  applyRelaxationLimit(dn2, dhi, dlo, relaxation);
  applyRelaxationLimit(dn3, dhi, dlo, relaxation);
  applyRelaxationLimit(dn4, dhi, dlo, relaxation);

  metrics.rmsContribution =
      (dn1.array().square() + dn2.array().square() + dn3.array().square() +
       dn4.array().square())
          .sum();

  double local_max = dn1.cwiseAbs().maxCoeff();
  local_max = std::max(local_max, dn2.cwiseAbs().maxCoeff());
  local_max = std::max(local_max, dn3.cwiseAbs().maxCoeff());
  local_max = std::max(local_max, dn4.cwiseAbs().maxCoeff());
  metrics.maxChange = local_max;

  return metrics;
}

BoundaryLayerSideProfiles BoundaryLayerMarcher::applyBoundaryLayerDelta(
    const BoundaryLayerWorkflow& workflow, int side, const BoundaryLayerDelta& delta, double relaxation,
    double hstinv, double gamm1) const {
  BoundaryLayerSideProfiles state;
  state.skinFrictionCoeff =
      workflow.lattice.get(side).profiles.skinFrictionCoeff;
  state.momentumThickness =
      workflow.lattice.get(side).profiles.momentumThickness;
  state.displacementThickness =
      workflow.lattice.get(side).profiles.displacementThickness;
  state.edgeVelocity = workflow.lattice.get(side).profiles.edgeVelocity;
  state.massFlux = workflow.lattice.get(side).profiles.massFlux;
  state.skinFrictionCoeffHistory =
      workflow.lattice.get(side).profiles.skinFrictionCoeffHistory;
  state.transitionIndex = workflow.lattice.get(side).profiles.transitionIndex;

  const int len = delta.dskinFrictionCoeff.size();
  if (len <= 0) {
    return state;
  }

  state.skinFrictionCoeff.head(len) += relaxation * delta.dskinFrictionCoeff;
  state.momentumThickness.head(len) += relaxation * delta.dmomentumThickness;
  state.displacementThickness.head(len) +=
      relaxation * delta.ddisplacementThickness;
  state.edgeVelocity.head(len) += relaxation * delta.dedgeVelocity;

  const int transition_index = std::max(0, workflow.lattice.get(side).profiles.transitionIndex);
  for (int idx = transition_index; idx < len; ++idx) {
    state.skinFrictionCoeff[idx] =
        std::min(state.skinFrictionCoeff[idx], 0.25);
  }

  for (int ibl = 0; ibl < len; ++ibl) {
    double dswaki = 0.0;
    if (ibl > workflow.lattice.get(side).trailingEdgeIndex) {
      const int wake_index =
          ibl - (workflow.lattice.get(side).trailingEdgeIndex + 1);
      dswaki = workflow.wgap[wake_index];
    }

    const double hklim =
        (ibl <= workflow.lattice.get(side).trailingEdgeIndex) ? 1.02 : 1.00005;
    const double edgeVelocity_val = state.edgeVelocity[ibl];
    const double edgeVelocity_sq = edgeVelocity_val * edgeVelocity_val;
    const double denom = 1.0 - 0.5 * edgeVelocity_sq * hstinv;
    const double msq = edgeVelocity_sq * hstinv / (gamm1 * denom);
    double dsw = state.displacementThickness[ibl] - dswaki;
    dsw = adjustDisplacementForHkLimit(
        dsw, state.momentumThickness[ibl], msq, hklim);
    state.displacementThickness[ibl] = dsw + dswaki;
    state.massFlux[ibl] =
        state.displacementThickness[ibl] * state.edgeVelocity[ibl];
  }

  return state;
}

void BoundaryLayerMarcher::syncStationRegimeStates(BoundaryLayerWorkflow& workflow, int side,
                                                    int stationIndex,
                                                    bool wake) {
  if (stationIndex < workflow.lattice.get(side).profiles.transitionIndex) {
    workflow.state.station2 =
        workflow.boundaryLayerVariablesSolver.solve(workflow.state.station2, FlowRegimeEnum::Laminar);
    workflow.blmid(FlowRegimeEnum::Laminar);
  }
  if (stationIndex >= workflow.lattice.get(side).profiles.transitionIndex) {
    workflow.state.station2 = workflow.boundaryLayerVariablesSolver.solve(workflow.state.station2,
                                                        FlowRegimeEnum::Turbulent);
    workflow.blmid(FlowRegimeEnum::Turbulent);
  }
  if (wake) {
    workflow.state.station2 =
        workflow.boundaryLayerVariablesSolver.solve(workflow.state.station2, FlowRegimeEnum::Wake);
    workflow.blmid(FlowRegimeEnum::Wake);
  }
  const bool similarity = (stationIndex == 0);
  workflow.flowRegime = determineRegimeForStation(workflow, side, stationIndex,
                                                  similarity, wake);
}

FlowRegimeEnum BoundaryLayerMarcher::determineRegimeForStation(
    const BoundaryLayerWorkflow& workflow, int side, int stationIndex, bool similarity, bool wake) const {
  if (wake) {
    return FlowRegimeEnum::Wake;
  }
  const int transitionIndex = workflow.lattice.get(side).profiles.transitionIndex;
  if (stationIndex == transitionIndex) {
    return FlowRegimeEnum::Transition;
  }
  if (stationIndex > transitionIndex) {
    return FlowRegimeEnum::Turbulent;
  }
  if (similarity) {
    return FlowRegimeEnum::Similarity;
  }
  return FlowRegimeEnum::Laminar;
}

bool BoundaryLayerMarcher::mrchdu(BoundaryLayerWorkflow& workflow, XFoil& xfoil) {
  return mrchdu(workflow, workflow.state, xfoil);
}

bool BoundaryLayerMarcher::mrchdu(BoundaryLayerWorkflow& workflow, BoundaryLayerState& state, XFoil& xfoil) {
  const double deps = 0.000005;
  const double senswt = 1000.0;

  double sens = 0.0;
  double sennew = 0.0;
  double ami = 0.0;

  for (int side = 1; side <= 2; ++side) {
    if (!marchBoundaryLayerSide(workflow, state, side, deps, senswt, sens,
                                sennew, ami, xfoil)) {
      return false;
    }
  }
  return true;
}

bool BoundaryLayerMarcher::marchBoundaryLayerSide(
    BoundaryLayerWorkflow& workflow, BoundaryLayerState& state, int side, double deps, double senswt,
    double& sens, double& sennew, double& ami, XFoil& xfoil) {
  const int previousTransition = resetSideState(workflow, side, xfoil);

  for (int stationIndex = 0;
       stationIndex < workflow.lattice.get(side).stationCount - 1; ++stationIndex) {
    if (!processBoundaryLayerStation(workflow, state, side, stationIndex,
                                     previousTransition, deps, senswt, sens,
                                     sennew, ami, xfoil)) {
      return false;
    }
  }

  return true;
}

bool BoundaryLayerMarcher::processBoundaryLayerStation(
    BoundaryLayerWorkflow& workflow, BoundaryLayerState& state, int side, int stationIndex,
    int previousTransition, double deps, double senswt, double& sens,
    double& sennew, double& ami, XFoil& xfoil) {
  BoundaryContext ctx =
      prepareMixedModeStation(workflow, side, stationIndex, previousTransition, ami);

  bool converged =
      xfoil.performMixedModeNewtonIteration(side, stationIndex,
                                            previousTransition, ctx, deps,
                                            senswt, sens, sennew, ami);
  if (!converged) {
    xfoil.handleMixedModeNonConvergence(side, stationIndex, ctx, ami);
  }

  sens = sennew;
  storeStationStateCommon(workflow, side, stationIndex, ctx.ami, ctx.cti, ctx.thi,
                          ctx.dsi, ctx.uei, ctx.xsi, ctx.dswaki);
  return true;
}

bool BoundaryLayerMarcher::mrchue(BoundaryLayerWorkflow& workflow, XFoil& xfoil) {
  return mrchue(workflow, workflow.state, xfoil);
}

bool BoundaryLayerMarcher::mrchue(BoundaryLayerWorkflow& workflow, BoundaryLayerState& state, XFoil& xfoil) {
  std::stringstream ss;
  for (int side = 1; side <= 2; ++side) {
    if (!marchMrchueSide(workflow, state, side, xfoil, ss)) {
      return false;
    }
  }
  return true;
}

bool BoundaryLayerMarcher::marchMrchueSide(BoundaryLayerWorkflow& workflow, BoundaryLayerState& state,
                                            int side, XFoil& xfoil,
                                            std::stringstream& ss) {
  ss << "    Side " << side << " ...\n";
  Logger::instance().write(ss.str());
  ss.str("");

  resetSideState(workflow, side, xfoil);

  double thi = 0.0;
  double dsi = 0.0;
  double ami = 0.0;
  double cti = 0.0;
  initializeMrchueSide(workflow, side, thi, dsi, ami, cti);

  for (int stationIndex = 0;
       stationIndex < workflow.lattice.get(side).stationCount - 1; ++stationIndex) {
    MrchueStationContext ctx;
    prepareMrchueStationContext(workflow, side, stationIndex, ctx, thi, dsi, ami, cti);
    bool converged =
        performMrchueNewtonLoop(workflow, side, stationIndex, ctx, xfoil, ss);
    if (!converged) {
      handleMrchueStationFailure(workflow, side, stationIndex, ctx, xfoil, ss);
    }
    storeMrchueStationState(workflow, side, stationIndex, ctx);

    ami = ctx.ami;
    cti = ctx.cti;
    thi = ctx.thi;
    dsi = ctx.dsi;

    if (stationIndex == workflow.lattice.get(side).trailingEdgeIndex) {
      thi = workflow.lattice.get(1).profiles.momentumThickness[workflow.lattice.top.trailingEdgeIndex] +
            workflow.lattice.get(2).profiles.momentumThickness[workflow.lattice.bottom.trailingEdgeIndex];
      dsi = workflow.lattice.get(1).profiles.displacementThickness[workflow.lattice.top.trailingEdgeIndex] +
            workflow.lattice.get(2).profiles.displacementThickness[workflow.lattice.bottom.trailingEdgeIndex] +
            xfoil.foil.edge.ante;
    }
  }

  return true;
}

void BoundaryLayerMarcher::initializeMrchueSide(BoundaryLayerWorkflow& workflow, int side, double& thi,
                                                 double& dsi, double& ami,
                                                 double& cti) {
  const double xsi = workflow.lattice.get(side).arcLengthCoordinates[0];
  const double uei = workflow.lattice.get(side).profiles.edgeVelocity[0];
  const double ucon = uei / xsi;
  const double tsq = 0.45 / (ucon * 6.0 * workflow.blReynolds.reybl);
  thi = std::sqrt(tsq);
  dsi = 2.2 * thi;
  ami = 0.0;
  cti = 0.03;
}

void BoundaryLayerMarcher::prepareMrchueStationContext(
    BoundaryLayerWorkflow& workflow, int side, int stationIndex, MrchueStationContext& ctx, double thi,
    double dsi, double ami, double cti) {
  ctx.simi = (stationIndex == 0);
  ctx.wake = stationIndex > workflow.lattice.get(side).trailingEdgeIndex;
  ctx.xsi = workflow.lattice.get(side).arcLengthCoordinates[stationIndex];
  ctx.uei = workflow.lattice.get(side).profiles.edgeVelocity[stationIndex];
  ctx.thi = thi;
  ctx.dsi = dsi;
  ctx.ami = ami;
  ctx.cti = cti;
  ctx.direct = true;
  ctx.dmax = 0.0;
  ctx.hmax = 0.0;
  ctx.htarg = 0.0;
  if (ctx.wake) {
    const int iw = stationIndex - workflow.lattice.get(side).trailingEdgeIndex;
    ctx.dswaki = workflow.wgap[iw - 1];
  } else {
    ctx.dswaki = 0.0;
  }
  workflow.flowRegime = determineRegimeForStation(workflow, side, stationIndex, ctx.simi,
                                                ctx.wake);
}

bool BoundaryLayerMarcher::performMrchueNewtonLoop(
    BoundaryLayerWorkflow& workflow, int side, int stationIndex, MrchueStationContext& ctx, XFoil& xfoil,
    std::stringstream& ss) {
  constexpr double kHlmax = 3.8;
  constexpr double kHtmax = 2.5;

  bool converged = false;
  bool direct = true;
  double htarg = 0.0;
  double dmax_local = 0.0;
  workflow.flowRegime = determineRegimeForStation(workflow, side, stationIndex, ctx.simi,
                                                ctx.wake);

  for (int itbl = 1; itbl <= 25; ++itbl) { 
    workflow.state.current() = workflow.blprv(workflow.state.current(), ctx.xsi, ctx.ami, ctx.cti,
                        ctx.thi, ctx.dsi, ctx.dswaki, ctx.uei);
    
    workflow.blkin(workflow.state);

    if ((!ctx.simi) && (!(workflow.flowRegime == FlowRegimeEnum::Turbulent || workflow.flowRegime == FlowRegimeEnum::Wake))) {
      workflow.transitionSolver.trchek(xfoil);
      ctx.ami = workflow.state.station2.param.amplz;

      if (workflow.flowRegime == FlowRegimeEnum::Transition) {
        workflow.lattice.get(side).profiles.transitionIndex = stationIndex;
        if (ctx.cti <= 0.0) {
          ctx.cti = 0.03;
          workflow.state.station2.param.sz = ctx.cti;
        }
      } else {
        workflow.lattice.get(side).profiles.transitionIndex = stationIndex + 2;
      }
    }

    if (stationIndex ==
        workflow.lattice.get(side).trailingEdgeIndex + 1) {
      ctx.tte = workflow.lattice.get(1).profiles.momentumThickness[workflow.lattice.top.trailingEdgeIndex] +
                workflow.lattice.get(2).profiles.momentumThickness[workflow.lattice.bottom.trailingEdgeIndex];
      ctx.dte = workflow.lattice.get(1).profiles.displacementThickness[workflow.lattice.top.trailingEdgeIndex] +
                workflow.lattice.get(2).profiles.displacementThickness[workflow.lattice.bottom.trailingEdgeIndex] +
                xfoil.foil.edge.ante;
      ctx.cte =
          (workflow.lattice.get(1).profiles.skinFrictionCoeff[workflow.lattice.top.trailingEdgeIndex] *
               workflow.lattice.get(1).profiles.momentumThickness[workflow.lattice.top.trailingEdgeIndex] +
           workflow.lattice.get(2).profiles.skinFrictionCoeff[workflow.lattice.bottom.trailingEdgeIndex] *
               workflow.lattice.get(2).profiles.momentumThickness[workflow.lattice.bottom.trailingEdgeIndex]) /
          ctx.tte;
      workflow.tesys(workflow.lattice.top.profiles, workflow.lattice.bottom.profiles, xfoil.foil.edge);
    } else {
      workflow.blsys();
    }

    if (direct) {
      workflow.blc.a2(3, 0) = 0.0;
      workflow.blc.a2(3, 1) = 0.0;
      workflow.blc.a2(3, 2) = 0.0;
      workflow.blc.a2(3, 3) = 1.0;
      workflow.blc.rhs[3] = 0.0;
      workflow.blc.rhs =
          workflow.blc.a2.block(0, 0, 4, 4).fullPivLu().solve(workflow.blc.rhs);

      dmax_local =
          std::max(std::fabs(workflow.blc.rhs[1] / ctx.thi),
                   std::fabs(workflow.blc.rhs[2] / ctx.dsi));
      if (stationIndex < workflow.lattice.get(side).profiles.transitionIndex) {
        dmax_local =
            std::max(dmax_local, std::fabs(workflow.blc.rhs[0] / 10.0));
      }
      if (stationIndex >= workflow.lattice.get(side).profiles.transitionIndex) {
        dmax_local =
            std::max(dmax_local, std::fabs(workflow.blc.rhs[0] / ctx.cti));
      }

      double rlx = 1.0;
      if (dmax_local > 0.3) {
        rlx = 0.3 / dmax_local;
      }

      if (stationIndex != workflow.lattice.get(side).trailingEdgeIndex + 1) {
        const double msq =
            ctx.uei * ctx.uei * workflow.blCompressibility.hstinv /
            (workflow.blCompressibility.gm1bl *
             (1.0 - 0.5 * ctx.uei * ctx.uei * workflow.blCompressibility.hstinv));
        const double htest =
            (ctx.dsi + rlx * workflow.blc.rhs[2]) /
            (ctx.thi + rlx * workflow.blc.rhs[1]);
        const auto hkin_result = boundary_layer::hkin(htest, msq);
        const double hktest = hkin_result.hk;

        const double hmax =
            (stationIndex < workflow.lattice.get(side).profiles.transitionIndex)
                ? kHlmax
                : kHtmax;
        direct = (hktest < hmax);
        if (!direct) {
          htarg = workflow.calcHtarg(stationIndex, side, ctx.wake);
          if (ctx.wake) {
            htarg = std::max(htarg, 1.01);
          } else {
            htarg = std::max(htarg, hmax);
          }

          ss << "     mrchue: inverse mode at " << stationIndex
             << "    hk=" << std::fixed << std::setprecision(3) << htarg
             << "\n";
          Logger::instance().write(ss.str());
          ss.str("");
          continue;
        }
      }

      if (stationIndex >= workflow.lattice.get(side).profiles.transitionIndex) {
        ctx.cti += rlx * workflow.blc.rhs[0];
      }
      ctx.thi += rlx * workflow.blc.rhs[1];
      ctx.dsi += rlx * workflow.blc.rhs[2];
      ctx.uei += rlx * workflow.blc.rhs[3];
    } else {
      workflow.blc.a2(3, 0) = 0.0;
      workflow.blc.a2(3, 1) = workflow.state.station2.hkz.t();
      workflow.blc.a2(3, 2) = workflow.state.station2.hkz.d();
      workflow.blc.a2(3, 3) = workflow.state.station2.hkz.u();
      workflow.blc.rhs[3] = ctx.htarg - workflow.state.station2.hkz.scalar;
      workflow.blc.rhs =
          workflow.blc.a2.block(0, 0, 4, 4).fullPivLu().solve(workflow.blc.rhs);

      dmax_local =
          std::max(std::fabs(workflow.blc.rhs[1] / ctx.thi),
                   std::fabs(workflow.blc.rhs[2] / ctx.dsi));
      if (stationIndex >= workflow.lattice.get(side).profiles.transitionIndex) {
        dmax_local =
            std::max(dmax_local, std::fabs(workflow.blc.rhs[0] / ctx.cti));
      }

      double rlx = 1.0;
      if (dmax_local > 0.3) {
        rlx = 0.3 / dmax_local;
      }

      if (stationIndex >= workflow.lattice.get(side).profiles.transitionIndex) {
        ctx.cti += rlx * workflow.blc.rhs[0];
      }
      ctx.thi += rlx * workflow.blc.rhs[1];
      ctx.dsi += rlx * workflow.blc.rhs[2];
      ctx.uei += rlx * workflow.blc.rhs[3];
    }

    if (stationIndex >= workflow.lattice.get(side).profiles.transitionIndex) {
      ctx.cti = std::clamp(ctx.cti, 0.0000001, 0.30);
    }

    const double hklim =
        (stationIndex <= workflow.lattice.get(side).trailingEdgeIndex) ? 1.02
                                                              : 1.00005;
    const double msq =
        ctx.uei * ctx.uei * workflow.blCompressibility.hstinv /
        (workflow.blCompressibility.gm1bl *
         (1.0 - 0.5 * ctx.uei * ctx.uei * workflow.blCompressibility.hstinv));
    double dsw = ctx.dsi - ctx.dswaki;
    dsw = adjustDisplacementForHkLimit(
        dsw, ctx.thi, msq, hklim);
    ctx.dsi = dsw + ctx.dswaki;

    if (dmax_local <= 0.00001) {
      converged = true;
      break;
    }
  }

  ctx.dmax = dmax_local;
  ctx.htarg = htarg;
  ctx.direct = direct;
  return converged;
}

void BoundaryLayerMarcher::handleMrchueStationFailure(
    BoundaryLayerWorkflow& workflow, int side, int stationIndex, MrchueStationContext& ctx, XFoil& xfoil,
    std::stringstream& ss) {
  workflow.flowRegime = determineRegimeForStation(workflow, side, stationIndex, ctx.simi,
                                                ctx.wake);

  ss << "     mrchue: convergence failed at " << stationIndex << ",  side "
     << side << ", res =" << std::fixed << std::setprecision(3) << ctx.dmax
     << "\n";
  Logger::instance().write(ss.str());
  ss.str("");

  resetStationKinematicsAfterFailure(
      workflow, side, stationIndex, ctx,
      EdgeVelocityFallbackMode::AverageNeighbors);

  {
    blData updatedCurrent =
        workflow.blprv(workflow.state.current(), ctx.xsi, ctx.ami, ctx.cti,
                        ctx.thi, ctx.dsi, ctx.dswaki, ctx.uei);
    workflow.state.current() = updatedCurrent;
  }
  workflow.blkin(workflow.state);

  workflow.checkTransitionIfNeeded(
      xfoil, side, stationIndex, ctx.simi, 2, ctx.ami);
  syncStationRegimeStates(workflow, side, stationIndex, ctx.wake);
}

void BoundaryLayerMarcher::storeMrchueStationState(
    BoundaryLayerWorkflow& workflow, int side, int stationIndex, const MrchueStationContext& ctx) {
  storeStationStateCommon(workflow, side, stationIndex, ctx.ami, ctx.cti, ctx.thi,
                          ctx.dsi, ctx.uei, ctx.xsi, ctx.dswaki);
}
