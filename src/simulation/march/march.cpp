#include "simulation/march/base.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "XFoil.h"
#include "core/math_util.hpp"
#include "infrastructure/logger.hpp"

using BoundaryContext = BoundaryLayerWorkflow::MixedModeStationContext;
using EdgeVelocityFallbackMode = BoundaryLayerWorkflow::EdgeVelocityFallbackMode;
using EdgeVelocityDistribution = BoundaryLayerWorkflow::EdgeVelocityDistribution;
using QtanResult = BoundaryLayerWorkflow::QtanResult;
using ClContributions = BoundaryLayerWorkflow::ClContributions;
using BoundaryLayerDelta = BoundaryLayerWorkflow::BoundaryLayerDelta;
using BoundaryLayerMetrics = BoundaryLayerWorkflow::BoundaryLayerMetrics;

int Marcher::resetSideState(BoundaryLayerWorkflow& workflow, int side,
                                         const Foil& foil,
                                         const StagnationResult& stagnation) {
  const int previousTransition = workflow.lattice.get(side).profiles.transitionIndex;
  workflow.blTransition.xiforc = workflow.xifset(foil, stagnation, side);
  workflow.flowRegime = FlowRegimeEnum::Laminar;
  workflow.lattice.get(side).profiles.transitionIndex = workflow.lattice.get(side).trailingEdgeIndex;
  return previousTransition;
}

void Marcher::storeStationStateCommon(
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

double Marcher::fallbackEdgeVelocity(
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
Marcher::computeNewUeDistribution(const BoundaryLayerWorkflow& workflow,
                                               const XFoil& xfoil,
                                               const Matrix3x2dVector& vdel) const {
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
                  vdel[systemIndex](2, 0));
          dui_ac += influence * (-vdel[systemIndex](2, 1));
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

BoundaryLayerWorkflow::QtanResult Marcher::computeQtan(
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
Marcher::computeClFromEdgeVelocityDistribution(
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
Marcher::buildBoundaryLayerDelta(
    const BoundaryLayerWorkflow& workflow, int side, const Eigen::VectorXd& unew_side,
    const Eigen::VectorXd& u_ac_side, double dac,
    const XFoil& xfoil, const Matrix3x2dVector& vdel) const {
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
        vdel[idx](0, 0) - dac * vdel[idx](0, 1);
    delta.dmomentumThickness[j] =
        vdel[idx](1, 0) - dac * vdel[idx](1, 1);
    dmass[j] = vdel[idx](2, 0) - dac * vdel[idx](2, 1);
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
Marcher::evaluateSegmentRelaxation(
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

BoundaryLayerSideProfiles Marcher::applyBoundaryLayerDelta(
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

void Marcher::syncStationRegimeStates(BoundaryLayerWorkflow& workflow, int side,
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

FlowRegimeEnum Marcher::determineRegimeForStation(
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
