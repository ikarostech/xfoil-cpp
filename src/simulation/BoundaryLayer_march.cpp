#include "BoundaryLayer.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "XFoil.h"

using BoundaryContext = BoundaryLayerWorkflow::MixedModeStationContext;

int BoundaryLayerWorkflow::resetSideState(int side, XFoil& xfoil) {
  const int previousTransition = lattice.get(side).profiles.transitionIndex;
  xfoil.blTransition.xiforc = xifset(xfoil.foil, xfoil.stagnation, side);
  flowRegime = FlowRegimeEnum::Laminar;
  lattice.get(side).profiles.transitionIndex = lattice.get(side).trailingEdgeIndex;
  return previousTransition;
}

void BoundaryLayerWorkflow::storeStationStateCommon(
    int side, int stationIndex, double ami, double cti, double thi,
    double dsi, double uei, double xsi, double dswaki, XFoil& xfoil) {
  if (stationIndex < lattice.get(side).profiles.transitionIndex) {
    lattice.get(side).profiles.skinFrictionCoeff[stationIndex] = ami;
  } else {
    lattice.get(side).profiles.skinFrictionCoeff[stationIndex] = cti;
  }
  lattice.get(side).profiles.momentumThickness[stationIndex] = thi;
  lattice.get(side).profiles.displacementThickness[stationIndex] = dsi;
  lattice.get(side).profiles.edgeVelocity[stationIndex] = uei;
  lattice.get(side).profiles.massFlux[stationIndex] = dsi * uei;
  lattice.get(side).profiles.skinFrictionCoeffHistory[stationIndex] = state.station2.cqz.scalar;

  {
    blData updatedCurrent =
        blprv(xfoil, state.current(), xsi, ami, cti, thi, dsi, dswaki, uei);
    state.current() = updatedCurrent;
  }
  blkin(xfoil, state);
  state.stepbl();

  if (flowRegime == FlowRegimeEnum::Wake) {
    // Keep wake state.
  } else if (flowRegime == FlowRegimeEnum::Transition ||
             stationIndex == lattice.get(side).trailingEdgeIndex) {
    flowRegime = FlowRegimeEnum::Turbulent;
  } else {
    flowRegime = FlowRegimeEnum::Laminar;
  }
}

BoundaryLayerWorkflow::MixedModeStationContext
BoundaryLayerWorkflow::prepareMixedModeStation(int side, int stationIndex,
                                               int previousTransition,
                                               double& ami) {
  MixedModeStationContext ctx;

  ctx.simi = (stationIndex == 0);
  ctx.wake = stationIndex > lattice.get(side).trailingEdgeIndex;
  ctx.xsi = lattice.get(side).arcLengthCoordinates[stationIndex];
  ctx.uei = lattice.get(side).profiles.edgeVelocity[stationIndex];
  ctx.thi = lattice.get(side).profiles.momentumThickness[stationIndex];
  ctx.dsi = lattice.get(side).profiles.displacementThickness[stationIndex];

  if (stationIndex < previousTransition) {
    ami = lattice.get(side).profiles.skinFrictionCoeff[stationIndex];
    ctx.cti = 0.03;
  } else {
    ctx.cti = lattice.get(side).profiles.skinFrictionCoeff[stationIndex];
    if (ctx.cti <= 0.0) {
      ctx.cti = 0.03;
    }
  }
  ctx.ami = ami;

  if (ctx.wake) {
    int iw = stationIndex - lattice.get(side).trailingEdgeIndex;
    ctx.dswaki = wgap[iw - 1];
  } else {
    ctx.dswaki = 0.0;
  }

  double thickness_limit =
      (stationIndex <= lattice.get(side).trailingEdgeIndex) ? 1.02 : 1.00005;
  ctx.dsi = std::max(ctx.dsi - ctx.dswaki, thickness_limit * ctx.thi) +
            ctx.dswaki;

  flowRegime =
      determineRegimeForStation(side, stationIndex, ctx.simi, ctx.wake);

  return ctx;
}

double BoundaryLayerWorkflow::fallbackEdgeVelocity(
    int side, int stationIndex,
    EdgeVelocityFallbackMode edgeMode) const {
  switch (edgeMode) {
    case EdgeVelocityFallbackMode::UsePreviousStation:
      return lattice.get(side).profiles.edgeVelocity[stationIndex - 1];
    case EdgeVelocityFallbackMode::AverageNeighbors: {
      double uei = lattice.get(side).profiles.edgeVelocity[stationIndex];
      if (stationIndex < lattice.get(side).stationCount - 1) {
        uei = 0.5 * (lattice.get(side).profiles.edgeVelocity[stationIndex - 1] +
                     lattice.get(side).profiles.edgeVelocity[stationIndex + 1]);
      }
      return uei;
    }
  }
  return lattice.get(side).profiles.edgeVelocity[stationIndex];
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

}  // namespace

BoundaryLayerWorkflow::EdgeVelocityDistribution
BoundaryLayerWorkflow::computeNewUeDistribution(const XFoil& xfoil) const {
  EdgeVelocityDistribution distribution;
  distribution.unew.top =
      Eigen::VectorXd::Zero(lattice.top.stationCount);
  distribution.unew.bottom =
      Eigen::VectorXd::Zero(lattice.bottom.stationCount);
  distribution.u_ac.top =
      Eigen::VectorXd::Zero(lattice.top.stationCount);
  distribution.u_ac.bottom =
      Eigen::VectorXd::Zero(lattice.bottom.stationCount);

  for (int side = 1; side <= 2; ++side) {
    for (int station = 0; station < lattice.get(side).stationCount - 1; ++station) {
      const int panelIndex = lattice.get(side).stationToPanel[station];
      double dui = 0.0;
      double dui_ac = 0.0;
      for (int otherSide = 1; otherSide <= 2; ++otherSide) {
        for (int otherStation = 0; otherStation < lattice.get(otherSide).stationCount - 1; ++otherStation) {
          const int otherPanel = lattice.get(otherSide).stationToPanel[otherStation];
          const int systemIndex = lattice.get(otherSide).stationToSystem[otherStation];
          const double influence =
              -lattice.get(side).panelInfluenceFactor[station] *
              lattice.get(otherSide).panelInfluenceFactor[otherStation] *
              xfoil.aerodynamicCache.dij(panelIndex, otherPanel);
          dui += influence *
                 (lattice.get(otherSide).profiles.massFlux[otherStation] +
                  xfoil.vdel[systemIndex](2, 0));
          dui_ac += influence * (-xfoil.vdel[systemIndex](2, 1));
        }
      }

      const double inviscidDerivative =
          xfoil.analysis_state_.controlByAlpha
              ? 0.0
              : lattice.get(side).inviscidEdgeVelocityMatrix(1, station);
      distribution.unew.get(side)[station] =
          lattice.get(side).inviscidEdgeVelocityMatrix(0, station) + dui;
      distribution.u_ac.get(side)[station] =
          inviscidDerivative + dui_ac;
    }
  }
  return distribution;
}

BoundaryLayerWorkflow::QtanResult BoundaryLayerWorkflow::computeQtan(
    const EdgeVelocityDistribution& distribution, int point_count) const {
  QtanResult result;
  result.qnew = Eigen::VectorXd::Zero(point_count);
  result.q_ac = Eigen::VectorXd::Zero(point_count);
  for (int side = 1; side <= 2; ++side) {
    const Eigen::VectorXd& unew_vec = distribution.unew.get(side);
    const Eigen::VectorXd& uac_vec = distribution.u_ac.get(side);
    const int limit = lattice.get(side).trailingEdgeIndex;
    for (int station = 0; station < limit; ++station) {
      const int panelIndex = lattice.get(side).stationToPanel[station];
      result.qnew[panelIndex] =
          lattice.get(side).panelInfluenceFactor[station] * unew_vec[station];
      result.q_ac[panelIndex] =
          lattice.get(side).panelInfluenceFactor[station] * uac_vec[station];
    }
  }
  return result;
}

BoundaryLayerWorkflow::ClContributions
BoundaryLayerWorkflow::computeClFromEdgeVelocityDistribution(
    const XFoil& xfoil, const EdgeVelocityDistribution& distribution) const {
  ClContributions contributions;
  const int point_count = xfoil.foil.foil_shape.n;
  const QtanResult qtan = computeQtan(distribution, point_count);

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
BoundaryLayerWorkflow::buildBoundaryLayerDelta(
    int side, const Eigen::VectorXd& unew_side,
    const Eigen::VectorXd& u_ac_side, double dac,
    const XFoil& xfoil) const {
  BoundaryLayerDelta delta;
  const int len = lattice.get(side).stationCount - 1;
  if (len <= 0) {
    return delta;
  }

  delta.dskinFrictionCoeff = Eigen::VectorXd(len);
  delta.dmomentumThickness = Eigen::VectorXd(len);
  delta.ddisplacementThickness = Eigen::VectorXd(len);
  delta.dedgeVelocity = Eigen::VectorXd(len);

  const auto iv = lattice.get(side).stationToSystem.segment(0, len);
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
      lattice.get(side).profiles.edgeVelocity.head(len);
  const Eigen::VectorXd displacementThickness_segment =
      lattice.get(side).profiles.displacementThickness.head(len);
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
BoundaryLayerWorkflow::evaluateSegmentRelaxation(
    int side, const BoundaryLayerDelta& delta, double dhi, double dlo,
    double& relaxation) const {
  BoundaryLayerMetrics metrics;
  const int len = delta.dskinFrictionCoeff.size();
  if (len <= 0) {
    return metrics;
  }

  const Eigen::VectorXd skinFrictionCoeff_segment =
      lattice.get(side).profiles.skinFrictionCoeff.head(len);
  const Eigen::VectorXd momentumThickness_segment =
      lattice.get(side).profiles.momentumThickness.head(len);
  const Eigen::VectorXd displacementThickness_segment =
      lattice.get(side).profiles.displacementThickness.head(len);

  Eigen::VectorXd dn1(len);
  const int transition_index = lattice.get(side).profiles.transitionIndex;
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

BoundaryLayerSideProfiles BoundaryLayerWorkflow::applyBoundaryLayerDelta(
    int side, const BoundaryLayerDelta& delta, double relaxation,
    double hstinv, double gamm1) const {
  BoundaryLayerSideProfiles state;
  state.skinFrictionCoeff =
      lattice.get(side).profiles.skinFrictionCoeff;
  state.momentumThickness =
      lattice.get(side).profiles.momentumThickness;
  state.displacementThickness =
      lattice.get(side).profiles.displacementThickness;
  state.edgeVelocity = lattice.get(side).profiles.edgeVelocity;
  state.massFlux = lattice.get(side).profiles.massFlux;
  state.skinFrictionCoeffHistory =
      lattice.get(side).profiles.skinFrictionCoeffHistory;
  state.transitionIndex = lattice.get(side).profiles.transitionIndex;

  const int len = delta.dskinFrictionCoeff.size();
  if (len <= 0) {
    return state;
  }

  state.skinFrictionCoeff.head(len) += relaxation * delta.dskinFrictionCoeff;
  state.momentumThickness.head(len) += relaxation * delta.dmomentumThickness;
  state.displacementThickness.head(len) +=
      relaxation * delta.ddisplacementThickness;
  state.edgeVelocity.head(len) += relaxation * delta.dedgeVelocity;

  const int transition_index = std::max(0, lattice.get(side).profiles.transitionIndex);
  for (int idx = transition_index; idx < len; ++idx) {
    state.skinFrictionCoeff[idx] =
        std::min(state.skinFrictionCoeff[idx], 0.25);
  }

  for (int ibl = 0; ibl < len; ++ibl) {
    double dswaki = 0.0;
    if (ibl > lattice.get(side).trailingEdgeIndex) {
      const int wake_index =
          ibl - (lattice.get(side).trailingEdgeIndex + 1);
      dswaki = wgap[wake_index];
    }

    const double hklim =
        (ibl <= lattice.get(side).trailingEdgeIndex) ? 1.02 : 1.00005;
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

void BoundaryLayerWorkflow::syncStationRegimeStates(int side,
                                                    int stationIndex,
                                                    bool wake,
                                                    XFoil& xfoil) {
  if (stationIndex < lattice.get(side).profiles.transitionIndex) {
    state.station2 = blvar(state.station2, FlowRegimeEnum::Laminar);
    blmid(FlowRegimeEnum::Laminar);
  }
  if (stationIndex >= lattice.get(side).profiles.transitionIndex) {
    state.station2 = blvar(state.station2, FlowRegimeEnum::Turbulent);
    blmid(FlowRegimeEnum::Turbulent);
  }
  if (wake) {
    state.station2 = blvar(state.station2, FlowRegimeEnum::Wake);
    blmid(FlowRegimeEnum::Wake);
  }
  const bool similarity = (stationIndex == 0);
  flowRegime = determineRegimeForStation(side, stationIndex,
                                                similarity, wake);
}

FlowRegimeEnum BoundaryLayerWorkflow::determineRegimeForStation(
    int side, int stationIndex, bool similarity, bool wake) const {
  if (wake) {
    return FlowRegimeEnum::Wake;
  }
  const int transitionIndex = lattice.get(side).profiles.transitionIndex;
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

bool BoundaryLayerWorkflow::mrchdu(XFoil& xfoil) {
  return mrchdu(state, xfoil);
}

bool BoundaryLayerWorkflow::mrchdu(BoundaryLayerState& state, XFoil& xfoil) {
  const double deps = 0.000005;
  const double senswt = 1000.0;

  double sens = 0.0;
  double sennew = 0.0;
  double ami = 0.0;

  for (int side = 1; side <= 2; ++side) {
    if (!marchBoundaryLayerSide(state, side, deps, senswt, sens, sennew, ami,
                                xfoil)) {
      return false;
    }
  }
  return true;
}

bool BoundaryLayerWorkflow::marchBoundaryLayerSide(
    BoundaryLayerState& state, int side, double deps, double senswt,
    double& sens, double& sennew, double& ami, XFoil& xfoil) {
  const int previousTransition = resetSideState(side, xfoil);

  for (int stationIndex = 0;
       stationIndex < lattice.get(side).stationCount - 1; ++stationIndex) {
    if (!processBoundaryLayerStation(state, side, stationIndex,
                                     previousTransition, deps, senswt, sens,
                                     sennew, ami, xfoil)) {
      return false;
    }
  }

  return true;
}

bool BoundaryLayerWorkflow::processBoundaryLayerStation(
    BoundaryLayerState& state, int side, int stationIndex,
    int previousTransition, double deps, double senswt, double& sens,
    double& sennew, double& ami, XFoil& xfoil) {
  BoundaryContext ctx =
      prepareMixedModeStation(side, stationIndex, previousTransition, ami);

  bool converged =
      xfoil.performMixedModeNewtonIteration(side, stationIndex,
                                            previousTransition, ctx, deps,
                                            senswt, sens, sennew, ami);
  if (!converged) {
    xfoil.handleMixedModeNonConvergence(side, stationIndex, ctx, ami);
  }

  sens = sennew;
  storeStationStateCommon(side, stationIndex, ctx.ami, ctx.cti, ctx.thi,
                          ctx.dsi, ctx.uei, ctx.xsi, ctx.dswaki, xfoil);

  if (XFoil::isCancelled()) {
    return false;
  }

  return true;
}

bool BoundaryLayerWorkflow::mrchue(XFoil& xfoil) {
  return mrchue(state, xfoil);
}

bool BoundaryLayerWorkflow::mrchue(BoundaryLayerState& state, XFoil& xfoil) {
  std::stringstream ss;
  for (int side = 1; side <= 2; ++side) {
    if (!marchMrchueSide(state, side, xfoil, ss)) {
      return false;
    }
  }
  return true;
}

bool BoundaryLayerWorkflow::marchMrchueSide(BoundaryLayerState& state,
                                            int side, XFoil& xfoil,
                                            std::stringstream& ss) {
  ss << "    Side " << side << " ...\n";
  xfoil.writeString(ss.str());
  ss.str("");

  resetSideState(side, xfoil);

  double thi = 0.0;
  double dsi = 0.0;
  double ami = 0.0;
  double cti = 0.0;
  initializeMrchueSide(side, thi, dsi, ami, cti, xfoil);

  for (int stationIndex = 0;
       stationIndex < lattice.get(side).stationCount - 1; ++stationIndex) {
    MrchueStationContext ctx;
    prepareMrchueStationContext(side, stationIndex, ctx, thi, dsi, ami, cti,
                                xfoil);
    bool converged =
        performMrchueNewtonLoop(side, stationIndex, ctx, xfoil, ss);
    if (!converged) {
      handleMrchueStationFailure(side, stationIndex, ctx, xfoil, ss);
    }
    storeMrchueStationState(side, stationIndex, ctx, xfoil);

    ami = ctx.ami;
    cti = ctx.cti;
    thi = ctx.thi;
    dsi = ctx.dsi;

    if (stationIndex == lattice.get(side).trailingEdgeIndex) {
      thi = lattice.get(1).profiles.momentumThickness[lattice.top.trailingEdgeIndex] +
            lattice.get(2).profiles.momentumThickness[lattice.bottom.trailingEdgeIndex];
      dsi = lattice.get(1).profiles.displacementThickness[lattice.top.trailingEdgeIndex] +
            lattice.get(2).profiles.displacementThickness[lattice.bottom.trailingEdgeIndex] +
            xfoil.foil.edge.ante;
    }

    if (XFoil::isCancelled()) {
      return false;
    }
  }

  return true;
}

void BoundaryLayerWorkflow::initializeMrchueSide(int side, double& thi,
                                                 double& dsi, double& ami,
                                                 double& cti, XFoil& xfoil) {
  const double xsi = lattice.get(side).arcLengthCoordinates[0];
  const double uei = lattice.get(side).profiles.edgeVelocity[0];
  const double ucon = uei / xsi;
  const double tsq = 0.45 / (ucon * 6.0 * xfoil.blReynolds.reybl);
  thi = std::sqrt(tsq);
  dsi = 2.2 * thi;
  ami = 0.0;
  cti = 0.03;
}

void BoundaryLayerWorkflow::prepareMrchueStationContext(
    int side, int stationIndex, MrchueStationContext& ctx, double thi,
    double dsi, double ami, double cti, XFoil& xfoil) {
  ctx.simi = (stationIndex == 0);
  ctx.wake = stationIndex > lattice.get(side).trailingEdgeIndex;
  ctx.xsi = lattice.get(side).arcLengthCoordinates[stationIndex];
  ctx.uei = lattice.get(side).profiles.edgeVelocity[stationIndex];
  ctx.thi = thi;
  ctx.dsi = dsi;
  ctx.ami = ami;
  ctx.cti = cti;
  ctx.direct = true;
  ctx.dmax = 0.0;
  ctx.hmax = 0.0;
  ctx.htarg = 0.0;
  if (ctx.wake) {
    const int iw = stationIndex - lattice.get(side).trailingEdgeIndex;
    ctx.dswaki = wgap[iw - 1];
  } else {
    ctx.dswaki = 0.0;
  }
  flowRegime = determineRegimeForStation(side, stationIndex, ctx.simi,
                                                ctx.wake);
}

bool BoundaryLayerWorkflow::performMrchueNewtonLoop(
    int side, int stationIndex, MrchueStationContext& ctx, XFoil& xfoil,
    std::stringstream& ss) {
  constexpr double kHlmax = 3.8;
  constexpr double kHtmax = 2.5;

  bool converged = false;
  bool direct = true;
  double htarg = 0.0;
  double dmax_local = 0.0;
  flowRegime = determineRegimeForStation(side, stationIndex, ctx.simi,
                                                ctx.wake);

  for (int itbl = 1; itbl <= 25; ++itbl) {
    {
      blData updatedCurrent =
          blprv(xfoil, state.current(), ctx.xsi, ctx.ami, ctx.cti, ctx.thi,
                ctx.dsi, ctx.dswaki, ctx.uei);
      state.current() = updatedCurrent;
    }
    blkin(xfoil, state);

    if ((!ctx.simi) && (!(flowRegime == FlowRegimeEnum::Turbulent || flowRegime == FlowRegimeEnum::Wake))) {
      trchek(xfoil);
      ctx.ami = state.station2.param.amplz;

      if (flowRegime == FlowRegimeEnum::Transition) {
        lattice.get(side).profiles.transitionIndex = stationIndex;
        if (ctx.cti <= 0.0) {
          ctx.cti = 0.03;
          state.station2.param.sz = ctx.cti;
        }
      } else {
        lattice.get(side).profiles.transitionIndex = stationIndex + 2;
      }
    }

    if (stationIndex ==
        lattice.get(side).trailingEdgeIndex + 1) {
      ctx.tte = lattice.get(1).profiles.momentumThickness[lattice.top.trailingEdgeIndex] +
                lattice.get(2).profiles.momentumThickness[lattice.bottom.trailingEdgeIndex];
      ctx.dte = lattice.get(1).profiles.displacementThickness[lattice.top.trailingEdgeIndex] +
                lattice.get(2).profiles.displacementThickness[lattice.bottom.trailingEdgeIndex] +
                xfoil.foil.edge.ante;
      ctx.cte =
          (lattice.get(1).profiles.skinFrictionCoeff[lattice.top.trailingEdgeIndex] *
               lattice.get(1).profiles.momentumThickness[lattice.top.trailingEdgeIndex] +
           lattice.get(2).profiles.skinFrictionCoeff[lattice.bottom.trailingEdgeIndex] *
               lattice.get(2).profiles.momentumThickness[lattice.bottom.trailingEdgeIndex]) /
          ctx.tte;
      tesys(lattice.top.profiles, lattice.bottom.profiles, xfoil.foil.edge);
    } else {
      blsys(xfoil);
    }

    if (direct) {
      blc.a2(3, 0) = 0.0;
      blc.a2(3, 1) = 0.0;
      blc.a2(3, 2) = 0.0;
      blc.a2(3, 3) = 1.0;
      blc.rhs[3] = 0.0;
      blc.rhs =
          blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);

      dmax_local =
          std::max(std::fabs(blc.rhs[1] / ctx.thi),
                   std::fabs(blc.rhs[2] / ctx.dsi));
      if (stationIndex < lattice.get(side).profiles.transitionIndex) {
        dmax_local =
            std::max(dmax_local, std::fabs(blc.rhs[0] / 10.0));
      }
      if (stationIndex >= lattice.get(side).profiles.transitionIndex) {
        dmax_local =
            std::max(dmax_local, std::fabs(blc.rhs[0] / ctx.cti));
      }

      double rlx = 1.0;
      if (dmax_local > 0.3) {
        rlx = 0.3 / dmax_local;
      }

      if (stationIndex != lattice.get(side).trailingEdgeIndex + 1) {
        const double msq =
            ctx.uei * ctx.uei * xfoil.blCompressibility.hstinv /
            (xfoil.blCompressibility.gm1bl *
             (1.0 - 0.5 * ctx.uei * ctx.uei * xfoil.blCompressibility.hstinv));
        const double htest =
            (ctx.dsi + rlx * blc.rhs[2]) /
            (ctx.thi + rlx * blc.rhs[1]);
        const auto hkin_result = boundary_layer::hkin(htest, msq);
        const double hktest = hkin_result.hk;

        const double hmax =
            (stationIndex < lattice.get(side).profiles.transitionIndex)
                ? kHlmax
                : kHtmax;
        direct = (hktest < hmax);
        if (!direct) {
          htarg = calcHtarg(stationIndex, side, ctx.wake);
          if (ctx.wake) {
            htarg = std::max(htarg, 1.01);
          } else {
            htarg = std::max(htarg, hmax);
          }

          ss << "     mrchue: inverse mode at " << stationIndex
             << "    hk=" << std::fixed << std::setprecision(3) << htarg
             << "\n";
          xfoil.writeString(ss.str());
          ss.str("");
          continue;
        }
      }

      if (stationIndex >= lattice.get(side).profiles.transitionIndex) {
        ctx.cti += rlx * blc.rhs[0];
      }
      ctx.thi += rlx * blc.rhs[1];
      ctx.dsi += rlx * blc.rhs[2];
      ctx.uei += rlx * blc.rhs[3];
    } else {
      blc.a2(3, 0) = 0.0;
      blc.a2(3, 1) = state.station2.hkz.t();
      blc.a2(3, 2) = state.station2.hkz.d();
      blc.a2(3, 3) = state.station2.hkz.u();
      blc.rhs[3] = ctx.htarg - state.station2.hkz.scalar;
      blc.rhs =
          blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);

      dmax_local =
          std::max(std::fabs(blc.rhs[1] / ctx.thi),
                   std::fabs(blc.rhs[2] / ctx.dsi));
      if (stationIndex >= lattice.get(side).profiles.transitionIndex) {
        dmax_local =
            std::max(dmax_local, std::fabs(blc.rhs[0] / ctx.cti));
      }

      double rlx = 1.0;
      if (dmax_local > 0.3) {
        rlx = 0.3 / dmax_local;
      }

      if (stationIndex >= lattice.get(side).profiles.transitionIndex) {
        ctx.cti += rlx * blc.rhs[0];
      }
      ctx.thi += rlx * blc.rhs[1];
      ctx.dsi += rlx * blc.rhs[2];
      ctx.uei += rlx * blc.rhs[3];
    }

    if (stationIndex >= lattice.get(side).profiles.transitionIndex) {
      ctx.cti = std::clamp(ctx.cti, 0.0000001, 0.30);
    }

    const double hklim =
        (stationIndex <= lattice.get(side).trailingEdgeIndex) ? 1.02
                                                              : 1.00005;
    const double msq =
        ctx.uei * ctx.uei * xfoil.blCompressibility.hstinv /
        (xfoil.blCompressibility.gm1bl *
         (1.0 - 0.5 * ctx.uei * ctx.uei * xfoil.blCompressibility.hstinv));
    double dsw = ctx.dsi - ctx.dswaki;
    dsw = adjustDisplacementForHkLimit(dsw, ctx.thi, msq, hklim);
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

void BoundaryLayerWorkflow::handleMrchueStationFailure(
    int side, int stationIndex, MrchueStationContext& ctx, XFoil& xfoil,
    std::stringstream& ss) {
  flowRegime = determineRegimeForStation(side, stationIndex, ctx.simi,
                                                ctx.wake);

  ss << "     mrchue: convergence failed at " << stationIndex << ",  side "
     << side << ", res =" << std::fixed << std::setprecision(3) << ctx.dmax
     << "\n";
  xfoil.writeString(ss.str());
  ss.str("");

  resetStationKinematicsAfterFailure(
      side, stationIndex, ctx,
      EdgeVelocityFallbackMode::AverageNeighbors);

  {
    blData updatedCurrent =
        blprv(xfoil, state.current(), ctx.xsi, ctx.ami, ctx.cti, ctx.thi,
              ctx.dsi, ctx.dswaki, ctx.uei);
    state.current() = updatedCurrent;
  }
  blkin(xfoil, state);

  xfoil.checkTransitionIfNeeded(side, stationIndex, ctx.simi, 2, ctx.ami);
  syncStationRegimeStates(side, stationIndex, ctx.wake, xfoil);
}

void BoundaryLayerWorkflow::storeMrchueStationState(
    int side, int stationIndex, const MrchueStationContext& ctx, XFoil& xfoil) {
  storeStationStateCommon(side, stationIndex, ctx.ami, ctx.cti, ctx.thi,
                          ctx.dsi, ctx.uei, ctx.xsi, ctx.dswaki, xfoil);
}
