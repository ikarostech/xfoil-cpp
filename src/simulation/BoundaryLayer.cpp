#include "BoundaryLayer.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

#include "XFoil.h"
#include "core/boundary_layer_util.hpp"
#include "core/math_util.hpp"
#include "domain/boundary_layer/boundary_layer_builder.hpp"
#include "domain/coefficient/aero_coefficients.hpp"
#include "domain/coefficient/skin_friction.hpp"
#include "domain/flow_state.hpp"
#include "infrastructure/logger.hpp"
#include "simulation/boundary_layer_setbl.hpp"

using BoundaryContext = BoundaryLayerWorkflow::MixedModeStationContext;
using Eigen::Matrix;
using Eigen::Vector;
using Eigen::Vector2d;
using Eigen::VectorXd;

namespace {
constexpr double kMixedModeConvergenceTolerance = 5.0e-6;

void applyRelaxationLimit(const Eigen::VectorXd &dn, double dhi, double dlo,
                          double &relaxation) {
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
}

double BoundaryLayerWorkflow::adjustDisplacementForHkLimit(
    double displacementThickness, double momentumThickness, double msq,
    double hklim) {
  const double h = displacementThickness / momentumThickness;

  boundary_layer::KineticShapeParameterResult hkin_result =
      boundary_layer::hkin(h, msq);

  const double dh = std::max(0.0, hklim - hkin_result.hk) / hkin_result.hk_h;
  return displacementThickness + dh * momentumThickness;
}

BoundaryLayerWorkflow::EdgeVelocityDistribution
BoundaryLayerWorkflow::computeNewUeDistribution(
    const XFoil &xfoil, const Matrix3x2dVector &vdel) const {
  EdgeVelocityDistribution distribution;
  distribution.unew.top = Eigen::VectorXd::Zero(lattice.top.stationCount);
  distribution.unew.bottom = Eigen::VectorXd::Zero(lattice.bottom.stationCount);
  distribution.u_ac.top = Eigen::VectorXd::Zero(lattice.top.stationCount);
  distribution.u_ac.bottom = Eigen::VectorXd::Zero(lattice.bottom.stationCount);

  for (int side = 1; side <= 2; ++side) {
    for (int station = 0; station < lattice.get(side).stationCount - 1;
         ++station) {
      const int panelIndex = lattice.get(side).stationToPanel[station];
      double dui = 0.0;
      double dui_ac = 0.0;
      for (int otherSide = 1; otherSide <= 2; ++otherSide) {
        for (int otherStation = 0;
             otherStation < lattice.get(otherSide).stationCount - 1;
             ++otherStation) {
          const int otherPanel =
              lattice.get(otherSide).stationToPanel[otherStation];
          const int systemIndex =
              lattice.get(otherSide).stationToSystem[otherStation];
          const double influence = -lattice.get(side).panelInfluenceFactor[station] *
                                   lattice.get(otherSide)
                                       .panelInfluenceFactor[otherStation] *
                                   xfoil.aerodynamicCache.dij(panelIndex, otherPanel);
          dui += influence * (lattice.get(otherSide).profiles.massFlux[otherStation] +
                              vdel[systemIndex](2, 0));
          dui_ac += influence * (-vdel[systemIndex](2, 1));
        }
      }

      const double inviscidDerivative =
          xfoil.analysis_state_.controlByAlpha
              ? 0.0
              : lattice.get(side).inviscidEdgeVelocityMatrix(1, station);
      distribution.unew.get(side)[station] =
          lattice.get(side).inviscidEdgeVelocityMatrix(0, station) + dui;
      distribution.u_ac.get(side)[station] = inviscidDerivative + dui_ac;
    }
  }
  return distribution;
}

BoundaryLayerWorkflow::ClContributions
BoundaryLayerWorkflow::computeClFromEdgeVelocityDistribution(
    const XFoil &xfoil, const EdgeVelocityDistribution &distribution) const {
  ClContributions contributions;
  const int point_count = xfoil.foil.foil_shape.n;
  if (point_count == 0) {
    return contributions;
  }

  Eigen::VectorXd qnew = Eigen::VectorXd::Zero(point_count);
  Eigen::VectorXd q_ac = Eigen::VectorXd::Zero(point_count);
  for (int side = 1; side <= 2; ++side) {
    const Eigen::VectorXd &unew_vec = distribution.unew.get(side);
    const Eigen::VectorXd &uac_vec = distribution.u_ac.get(side);
    const int limit = lattice.get(side).trailingEdgeIndex;
    for (int station = 0; station < limit; ++station) {
      const int panelIndex = lattice.get(side).stationToPanel[station];
      qnew[panelIndex] = lattice.get(side).panelInfluenceFactor[station] *
                         unew_vec[station];
      q_ac[panelIndex] = lattice.get(side).panelInfluenceFactor[station] *
                         uac_vec[station];
    }
  }

  const auto compressibility = xfoil.buildCompressibilityParams();
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
        MathUtil::getRotateMatrix(xfoil.analysis_state_.alpha) *
        (xfoil.foil.foil_shape.points.col(ip) -
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

int BoundaryLayerWorkflow::resetSideState(int side, const Foil &foil,
                                          const StagnationResult &stagnation) {
  const int previousTransition = lattice.get(side).profiles.transitionIndex;
  blTransition.xiforc = xifset(foil, stagnation, side);
  flowRegime = FlowRegimeEnum::Laminar;
  lattice.get(side).profiles.transitionIndex = lattice.get(side).trailingEdgeIndex;
  return previousTransition;
}

BoundaryLayerWorkflow::StationReadModel
BoundaryLayerWorkflow::readStationModel(int side, int stationIndex) const {
  const auto &sideLattice = lattice.get(side);

  StationReadModel model;
  model.stationCount = sideLattice.stationCount;
  model.trailingEdgeIndex = sideLattice.trailingEdgeIndex;
  model.transitionIndex = sideLattice.profiles.transitionIndex;
  model.arcLength = sideLattice.arcLengthCoordinates[stationIndex];
  model.edgeVelocity = sideLattice.profiles.edgeVelocity[stationIndex];
  model.momentumThickness = sideLattice.profiles.momentumThickness[stationIndex];
  model.displacementThickness =
      sideLattice.profiles.displacementThickness[stationIndex];
  model.skinFrictionCoeff = sideLattice.profiles.skinFrictionCoeff[stationIndex];

  if (stationIndex > sideLattice.trailingEdgeIndex) {
    const int wakeIndex = stationIndex - sideLattice.trailingEdgeIndex;
    model.wakeGap = wgap[wakeIndex - 1];
  }

  return model;
}

int BoundaryLayerWorkflow::readSideStationCount(int side) const {
  return lattice.get(side).stationCount;
}

BoundaryLayerWorkflow::TrailingEdgeReadModel
BoundaryLayerWorkflow::readTrailingEdgeModel() const {
  TrailingEdgeReadModel model;
  model.topTrailingEdgeIndex = lattice.top.trailingEdgeIndex;
  model.bottomTrailingEdgeIndex = lattice.bottom.trailingEdgeIndex;

  model.topMomentumThickness =
      lattice.top.profiles.momentumThickness[model.topTrailingEdgeIndex];
  model.bottomMomentumThickness =
      lattice.bottom.profiles.momentumThickness[model.bottomTrailingEdgeIndex];
  model.topDisplacementThickness =
      lattice.top.profiles.displacementThickness[model.topTrailingEdgeIndex];
  model.bottomDisplacementThickness =
      lattice.bottom.profiles.displacementThickness[model.bottomTrailingEdgeIndex];
  model.topSkinFrictionCoeff =
      lattice.top.profiles.skinFrictionCoeff[model.topTrailingEdgeIndex];
  model.bottomSkinFrictionCoeff =
      lattice.bottom.profiles.skinFrictionCoeff[model.bottomTrailingEdgeIndex];
  return model;
}

FlowRegimeEnum
BoundaryLayerWorkflow::applyFlowRegimeCandidate(FlowRegimeEnum candidate) {
  flowRegime = candidate;
  return flowRegime;
}

FlowRegimeEnum BoundaryLayerWorkflow::currentFlowRegime() const {
  return flowRegime;
}

void BoundaryLayerWorkflow::emitMarchInfoLog(std::string_view message) const {
  Logger::instance().write(std::string(message));
}

void BoundaryLayerWorkflow::emitMarchFailureLog(std::string_view phase, int side,
                                                int stationIndex,
                                                double residual) const {
  std::stringstream ss;
  ss << "     " << phase << ": convergence failed at " << stationIndex
     << ",  side " << side << ", res =" << std::fixed << std::setprecision(3)
     << residual << "\n";
  Logger::instance().write(ss.str());
}

double BoundaryLayerWorkflow::readNewtonRhs(int row) const {
  return blc.rhs[row];
}

void BoundaryLayerWorkflow::solveMrchueDirectNewtonSystem() {
  blc.a2(3, 0) = 0.0;
  blc.a2(3, 1) = 0.0;
  blc.a2(3, 2) = 0.0;
  blc.a2(3, 3) = 1.0;
  blc.rhs[3] = 0.0;
  blc.rhs = blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);
}

void BoundaryLayerWorkflow::solveMrchueInverseNewtonSystem(double htarg) {
  blc.a2(3, 0) = 0.0;
  blc.a2(3, 1) = state.station2.hkz.t();
  blc.a2(3, 2) = state.station2.hkz.d();
  blc.a2(3, 3) = state.station2.hkz.u();
  blc.rhs[3] = htarg - state.station2.hkz.scalar;
  blc.rhs = blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);
}

void BoundaryLayerWorkflow::storeStationStateCommon(
    int side, int stationIndex, const MixedModeStationContext &ctx) {
  if (stationIndex < lattice.get(side).profiles.transitionIndex) {
    lattice.get(side).profiles.skinFrictionCoeff[stationIndex] = ctx.ami;
  } else {
    lattice.get(side).profiles.skinFrictionCoeff[stationIndex] = ctx.cti;
  }
  lattice.get(side).profiles.momentumThickness[stationIndex] = ctx.thi;
  lattice.get(side).profiles.displacementThickness[stationIndex] = ctx.dsi;
  lattice.get(side).profiles.edgeVelocity[stationIndex] = ctx.uei;
  lattice.get(side).profiles.massFlux[stationIndex] = ctx.dsi * ctx.uei;
  lattice.get(side).profiles.skinFrictionCoeffHistory[stationIndex] =
      state.station2.cqz.scalar;

  state.current() =
      blprv(state.current(), ctx.xsi, ctx.ami, ctx.cti, ctx.thi, ctx.dsi,
            ctx.dswaki, ctx.uei);
  blkin(state);
  state.stepbl();

  if (flowRegime == FlowRegimeEnum::Wake) {
    return;
  }
  if (flowRegime == FlowRegimeEnum::Transition ||
      stationIndex == lattice.get(side).trailingEdgeIndex) {
    flowRegime = FlowRegimeEnum::Turbulent;
  } else {
    flowRegime = FlowRegimeEnum::Laminar;
  }
}

double BoundaryLayerWorkflow::fallbackEdgeVelocity(
    int side, int stationIndex, EdgeVelocityFallbackMode edgeMode) const {
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

BoundaryLayerWorkflow::BoundaryLayerDelta
BoundaryLayerWorkflow::buildBoundaryLayerDelta(
    int side, const Eigen::VectorXd &unew_side, const Eigen::VectorXd &u_ac_side,
    double dac, const Matrix3x2dVector &vdel) const {
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
    delta.dskinFrictionCoeff[j] = vdel[idx](0, 0) - dac * vdel[idx](0, 1);
    delta.dmomentumThickness[j] = vdel[idx](1, 0) - dac * vdel[idx](1, 1);
    dmass[j] = vdel[idx](2, 0) - dac * vdel[idx](2, 1);
  }

  const Eigen::VectorXd edgeVelocity_segment =
      lattice.get(side).profiles.edgeVelocity.head(len);
  const Eigen::VectorXd displacementThickness_segment =
      lattice.get(side).profiles.displacementThickness.head(len);
  const Eigen::VectorXd unew_segment = unew_side.head(len);
  const Eigen::VectorXd uac_segment = u_ac_side.head(len);

  delta.dedgeVelocity = unew_segment + dac * uac_segment - edgeVelocity_segment;
  delta.ddisplacementThickness =
      (dmass - displacementThickness_segment.cwiseProduct(delta.dedgeVelocity))
          .cwiseQuotient(edgeVelocity_segment);

  return delta;
}

BoundaryLayerWorkflow::BoundaryLayerMetrics
BoundaryLayerWorkflow::evaluateSegmentRelaxation(
    int side, const BoundaryLayerDelta &delta, double dhi, double dlo,
    double &relaxation) const {
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

  metrics.rmsContribution = (dn1.array().square() + dn2.array().square() +
                             dn3.array().square() + dn4.array().square())
                                .sum();
  double local_max = dn1.cwiseAbs().maxCoeff();
  local_max = std::max(local_max, dn2.cwiseAbs().maxCoeff());
  local_max = std::max(local_max, dn3.cwiseAbs().maxCoeff());
  local_max = std::max(local_max, dn4.cwiseAbs().maxCoeff());
  metrics.maxChange = local_max;
  return metrics;
}

BoundaryLayerSideProfiles BoundaryLayerWorkflow::applyBoundaryLayerDelta(
    int side, const BoundaryLayerDelta &delta, double relaxation, double hstinv,
    double gamm1) const {
  BoundaryLayerSideProfiles updated = lattice.get(side).profiles;
  const int len = delta.dskinFrictionCoeff.size();
  if (len <= 0) {
    return updated;
  }

  updated.skinFrictionCoeff.head(len) += relaxation * delta.dskinFrictionCoeff;
  updated.momentumThickness.head(len) += relaxation * delta.dmomentumThickness;
  updated.displacementThickness.head(len) +=
      relaxation * delta.ddisplacementThickness;
  updated.edgeVelocity.head(len) += relaxation * delta.dedgeVelocity;

  const int transition_index = std::max(0, lattice.get(side).profiles.transitionIndex);
  for (int idx = transition_index; idx < len; ++idx) {
    updated.skinFrictionCoeff[idx] = std::min(updated.skinFrictionCoeff[idx], 0.25);
  }

  for (int ibl = 0; ibl < len; ++ibl) {
    double dswaki = 0.0;
    if (ibl > lattice.get(side).trailingEdgeIndex) {
      const int wake_index = ibl - (lattice.get(side).trailingEdgeIndex + 1);
      dswaki = wgap[wake_index];
    }

    const double hklim = (ibl <= lattice.get(side).trailingEdgeIndex) ? 1.02 : 1.00005;
    const double edgeVelocity_val = updated.edgeVelocity[ibl];
    const double edgeVelocity_sq = edgeVelocity_val * edgeVelocity_val;
    const double denom = 1.0 - 0.5 * edgeVelocity_sq * hstinv;
    const double msq = edgeVelocity_sq * hstinv / (gamm1 * denom);
    double dsw = updated.displacementThickness[ibl] - dswaki;
    dsw = adjustDisplacementForHkLimit(dsw, updated.momentumThickness[ibl], msq,
                                       hklim);
    updated.displacementThickness[ibl] = dsw + dswaki;
    updated.massFlux[ibl] =
        updated.displacementThickness[ibl] * updated.edgeVelocity[ibl];
  }

  return updated;
}

void BoundaryLayerWorkflow::syncStationRegimeStates(int side, int stationIndex,
                                                    FlowRegimeEnum stationRegime) {
  if (stationIndex < lattice.get(side).profiles.transitionIndex) {
    state.station2 =
        boundaryLayerVariablesSolver.solve(state.station2, FlowRegimeEnum::Laminar);
    blmid(FlowRegimeEnum::Laminar);
  }
  if (stationIndex >= lattice.get(side).profiles.transitionIndex) {
    state.station2 = boundaryLayerVariablesSolver.solve(
        state.station2, FlowRegimeEnum::Turbulent);
    blmid(FlowRegimeEnum::Turbulent);
  }
  if (stationRegime == FlowRegimeEnum::Wake) {
    state.station2 =
        boundaryLayerVariablesSolver.solve(state.station2, FlowRegimeEnum::Wake);
    blmid(FlowRegimeEnum::Wake);
  }
  flowRegime = determineRegimeForStation(side, stationIndex);
}

FlowRegimeEnum BoundaryLayerWorkflow::determineRegimeForStation(
    int side, int stationIndex) const {
  if (stationIndex > lattice.get(side).trailingEdgeIndex) {
    return FlowRegimeEnum::Wake;
  }
  const int transitionIndex = lattice.get(side).profiles.transitionIndex;
  if (stationIndex == transitionIndex) {
    return FlowRegimeEnum::Transition;
  }
  if (stationIndex > transitionIndex) {
    return FlowRegimeEnum::Turbulent;
  }
  if (stationIndex == 0) {
    return FlowRegimeEnum::Similarity;
  }
  return FlowRegimeEnum::Laminar;
}

bool BoundaryLayerWorkflow::blkin(BoundaryLayerState &state) {
  //----------------------------------------------------------
  //     calculates turbulence-independent secondary "2"
  //     variables from the primary "2" variables.
  //----------------------------------------------------------
  // BlCompressibilityParams& blCompressibility = this->blCompressibility;
  // BlReynoldsParams& blReynolds = this->blReynolds;
  blData &current = state.current();
  //---- set edge mach number ** 2
  current.param.mz = current.param.uz * current.param.uz *
                     blCompressibility.hstinv /
                     (blCompressibility.gm1bl *
                      (1.0 - 0.5 * current.param.uz * current.param.uz *
                                 blCompressibility.hstinv));
  double tr2 = 1.0 + 0.5 * blCompressibility.gm1bl * current.param.mz;
  current.param.mz_uz = 2.0 * current.param.mz * tr2 / current.param.uz;
  current.param.mz_ms = current.param.uz * current.param.uz * tr2 /
                        (blCompressibility.gm1bl *
                         (1.0 - 0.5 * current.param.uz * current.param.uz *
                                    blCompressibility.hstinv)) *
                        blCompressibility.hstinv_ms;

  //---- set edge density (isentropic relation)
  current.param.rz =
      blCompressibility.rstbl * pow(tr2, (-1.0 / blCompressibility.gm1bl));
  current.param.rz_uz = -current.param.rz / tr2 * 0.5 * current.param.mz_uz;
  current.param.rz_ms =
      -current.param.rz / tr2 * 0.5 * current.param.mz_ms +
      blCompressibility.rstbl_ms * pow(tr2, (-1.0 / blCompressibility.gm1bl));

  //---- set shape parameter
  current.param.hz = current.param.dz / current.param.tz;
  current.param.hz_dz = 1.0 / current.param.tz;
  current.param.hz_tz = -current.param.hz / current.param.tz;

  //---- set edge static/stagnation enthalpy
  double herat = 1.0 - 0.5 * current.param.uz * current.param.uz *
                           blCompressibility.hstinv;
  double he_u2 = -current.param.uz * blCompressibility.hstinv;
  double he_ms =
      -0.5 * current.param.uz * current.param.uz * blCompressibility.hstinv_ms;
  //---- set molecular viscosity
  double v2_he = (1.5 / herat - 1.0 / (herat + kHvrat));

  //---- set kinematic shape parameter
  boundary_layer::KineticShapeParameterResult hkin_result =
      boundary_layer::hkin(current.param.hz, current.param.mz);
  current.hkz.scalar = hkin_result.hk;

  current.hkz.u() = hkin_result.hk_msq * current.param.mz_uz;
  current.hkz.t() = hkin_result.hk_h * current.param.hz_tz;
  current.hkz.d() = hkin_result.hk_h * current.param.hz_dz;
  current.hkz.ms() = hkin_result.hk_msq * current.param.mz_ms;

  //---- set momentum thickness reynolds number
  current.rtz.scalar = current.param.rz * current.param.uz * current.param.tz /
                       (sqrt(herat * herat * herat) * (1.0 + kHvrat) /
                        (herat + kHvrat) / blReynolds.reybl);
  current.rtz.u() = current.rtz.scalar *
                    (1.0 / current.param.uz +
                     current.param.rz_uz / current.param.rz - v2_he * he_u2);
  current.rtz.t() = current.rtz.scalar / current.param.tz;
  current.rtz.ms() =
      current.rtz.scalar *
      (current.param.rz_ms / current.param.rz +
       (1 / blReynolds.reybl * blReynolds.reybl_ms - v2_he * he_ms));
  current.rtz.re() =
      current.rtz.scalar * (blReynolds.reybl_re / blReynolds.reybl);

  return true;
}

bool BoundaryLayerWorkflow::isStartOfWake(int side, int stationIndex) {
  return stationIndex == lattice.get(side).trailingEdgeIndex + 1;
}

void BoundaryLayerWorkflow::updateSystemMatricesForStation(
    const Edge &edge, int side, int stationIndex, BoundaryContext &ctx) {
  if (isStartOfWake(side, stationIndex)) {
    ctx.tte = lattice.get(1)
                  .profiles.momentumThickness[lattice.top.trailingEdgeIndex] +
              lattice.get(2)
                  .profiles.momentumThickness[lattice.bottom.trailingEdgeIndex];
    ctx.dte =
        lattice.get(1)
            .profiles.displacementThickness[lattice.top.trailingEdgeIndex] +
        lattice.get(2)
            .profiles.displacementThickness[lattice.bottom.trailingEdgeIndex] +
        edge.ante;
    ctx.cte =
        (lattice.get(1)
                 .profiles.skinFrictionCoeff[lattice.top.trailingEdgeIndex] *
             lattice.get(1)
                 .profiles.momentumThickness[lattice.top.trailingEdgeIndex] +
         lattice.get(2)
                 .profiles.skinFrictionCoeff[lattice.bottom.trailingEdgeIndex] *
             lattice.get(2)
                 .profiles
                 .momentumThickness[lattice.bottom.trailingEdgeIndex]) /
        ctx.tte;
    tesys(lattice.top.profiles, lattice.bottom.profiles, edge);
  } else {
    blsys();
  }
}

void BoundaryLayerWorkflow::initializeFirstIterationState(
    int side, int stationIndex, int previousTransition, BoundaryContext &ctx,
    double &ueref, double &hkref) {
  ueref = state.station2.param.uz;
  hkref = state.station2.hkz.scalar;

  if (stationIndex < lattice.get(side).profiles.transitionIndex &&
      stationIndex >= previousTransition) {
    const double uem =
        lattice.get(side).profiles.edgeVelocity[std::max(0, stationIndex - 1)];
    const double dsm =
        lattice.get(side)
            .profiles.displacementThickness[std::max(0, stationIndex - 1)];
    const double thm =
        lattice.get(side)
            .profiles.momentumThickness[std::max(0, stationIndex - 1)];

    const double uem_sq = uem * uem;
    const double msq = uem_sq * blCompressibility.hstinv /
                       (blCompressibility.gm1bl *
                        (1.0 - 0.5 * uem_sq * blCompressibility.hstinv));
    const auto hkin_result = boundary_layer::hkin(dsm / thm, msq);
    hkref = hkin_result.hk;
  }

  if (stationIndex < previousTransition) {
    if (flowRegime == FlowRegimeEnum::Transition) {
      lattice.get(side).profiles.skinFrictionCoeff[stationIndex] = 0.03;
    } else {
      lattice.get(side).profiles.skinFrictionCoeff[stationIndex] =
          lattice.get(side)
              .profiles.skinFrictionCoeff[std::max(0, stationIndex - 1)];
    }
    ctx.cti = lattice.get(side)
                  .profiles.skinFrictionCoeff[std::max(0, stationIndex - 1)];
    state.station2.param.sz = ctx.cti;
  }
}

void BoundaryLayerWorkflow::configureSimilarityRow(double ueref) {
  blc.a2(3, 0) = 0.0;
  blc.a2(3, 1) = 0.0;
  blc.a2(3, 2) = 0.0;
  blc.a2(3, 3) = state.station2.param.uz_uei;
  blc.rhs[3] = ueref - state.station2.param.uz;
}

void BoundaryLayerWorkflow::configureViscousRow(double hkref, double ueref,
                                                double senswt,
                                                bool resetSensitivity,
                                                bool averageSensitivity,
                                                double &sens, double &sennew) {
  blc.a2(3, 0) = 0.0;
  blc.a2(3, 1) = state.station2.hkz.t();
  blc.a2(3, 2) = state.station2.hkz.d();
  blc.a2(3, 3) = state.station2.hkz.u() * state.station2.param.uz_uei;
  blc.rhs[3] = 1.0;

  const double delta_sen =
      blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs)[3];

  sennew = senswt * delta_sen * hkref / ueref;
  if (resetSensitivity) {
    sens = sennew;
  } else if (averageSensitivity) {
    sens = 0.5 * (sens + sennew);
  }

  blc.a2(3, 1) = state.station2.hkz.t() * hkref;
  blc.a2(3, 2) = state.station2.hkz.d() * hkref;
  blc.a2(3, 3) = (state.station2.hkz.u() * hkref + sens / ueref) *
                 state.station2.param.uz_uei;
  blc.rhs[3] = -(hkref * hkref) * (state.station2.hkz.scalar / hkref - 1.0) -
               sens * (state.station2.param.uz / ueref - 1.0);
}

bool BoundaryLayerWorkflow::applyMixedModeNewtonStep(int side, int stationIndex,
                                                     double &ami,
                                                     BoundaryContext &ctx) {
  blc.rhs = blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);

  ctx.dmax = std::max(std::fabs(blc.rhs[1] / ctx.thi),
                      std::fabs(blc.rhs[2] / ctx.dsi));
  if (stationIndex >= lattice.get(side).profiles.transitionIndex) {
    ctx.dmax = std::max(ctx.dmax, std::fabs(blc.rhs[0] / (10.0 * ctx.cti)));
  }

  double rlx = 1.0;
  if (ctx.dmax > 0.3) {
    rlx = 0.3 / ctx.dmax;
  }

  if (stationIndex < lattice.get(side).profiles.transitionIndex) {
    ami += rlx * blc.rhs[0];
    ctx.ami = ami;
  }
  if (stationIndex >= lattice.get(side).profiles.transitionIndex) {
    ctx.cti += rlx * blc.rhs[0];
  }
  ctx.thi += rlx * blc.rhs[1];
  ctx.dsi += rlx * blc.rhs[2];
  ctx.uei += rlx * blc.rhs[3];

  if (stationIndex >= lattice.get(side).profiles.transitionIndex) {
    ctx.cti = std::clamp(ctx.cti, 0.0000001, 0.30);
  }

  const double hklim =
      (stationIndex <= lattice.get(side).trailingEdgeIndex) ? 1.02 : 1.00005;
  const double uei_sq = ctx.uei * ctx.uei;
  const double msq = uei_sq * blCompressibility.hstinv /
                     (blCompressibility.gm1bl *
                      (1.0 - 0.5 * uei_sq * blCompressibility.hstinv));
  double dsw = ctx.dsi - ctx.dswaki;
  dsw = adjustDisplacementForHkLimit(dsw, ctx.thi, msq, hklim);
  ctx.dsi = dsw + ctx.dswaki;

  return ctx.dmax <= kMixedModeConvergenceTolerance;
}

SkinFrictionCoefficients
BoundaryLayerWorkflow::blmid(FlowRegimeEnum flowRegimeType) {
  blData &previous = state.previous();
  blData &current = state.current();

  if (flowRegimeType == FlowRegimeEnum::Similarity) {
    previous.hkz = current.hkz;
    previous.rtz = current.rtz;
    previous.param.mz = current.param.mz;
    previous.param.mz_uz = current.param.mz_uz;
    previous.param.mz_ms = current.param.mz_ms;
  }

  const double hka = 0.5 * (previous.hkz.scalar + current.hkz.scalar);
  const double rta = 0.5 * (previous.rtz.scalar + current.rtz.scalar);
  const double ma = 0.5 * (previous.param.mz + current.param.mz);

  skin_friction::C_f cf_res =
      skin_friction::getSkinFriction(hka, rta, ma, flowRegimeType);

  SkinFrictionCoefficients coeffs;
  coeffs.cfm = cf_res.cf;
  const double cfm_hka = cf_res.hk;
  const double cfm_rta = cf_res.rt;
  const double cfm_ma = cf_res.msq;

  coeffs.cfm_u1 =
      0.5 * (cfm_hka * previous.hkz.u() + cfm_ma * previous.param.mz_uz +
             cfm_rta * previous.rtz.u());
  coeffs.cfm_t1 =
      0.5 * (cfm_hka * previous.hkz.t() + cfm_rta * previous.rtz.t());
  coeffs.cfm_d1 = 0.5 * (cfm_hka * previous.hkz.d());

  coeffs.cfm_u2 =
      0.5 * (cfm_hka * current.hkz.u() + cfm_ma * current.param.mz_uz +
             cfm_rta * current.rtz.u());
  coeffs.cfm_t2 = 0.5 * (cfm_hka * current.hkz.t() + cfm_rta * current.rtz.t());
  coeffs.cfm_d2 = 0.5 * (cfm_hka * current.hkz.d());

  coeffs.cfm_ms =
      0.5 * (cfm_hka * previous.hkz.ms() + cfm_ma * previous.param.mz_ms +
             cfm_rta * previous.rtz.ms() + cfm_hka * current.hkz.ms() +
             cfm_ma * current.param.mz_ms + cfm_rta * current.rtz.ms());
  coeffs.cfm_re =
      0.5 * (cfm_rta * previous.rtz.re() + cfm_rta * current.rtz.re());

  return coeffs;
}

blData BoundaryLayerWorkflow::blprv(blData data, double xsi, double ami,
                                    double cti, double thi, double dsi,
                                    double dswaki, double uei) const {
  data.param.xz = xsi;
  data.param.amplz = ami;
  data.param.sz = cti;
  data.param.tz = thi;
  data.param.dz = dsi - dswaki;
  data.param.dwz = dswaki;

  data.param.uz =
      uei * (1.0 - blCompressibility.tkbl) /
      (1.0 - blCompressibility.tkbl * (uei / blCompressibility.qinfbl) *
                 (uei / blCompressibility.qinfbl));
  data.param.uz_uei =
      (1.0 + blCompressibility.tkbl *
                 (2.0 * data.param.uz * uei / blCompressibility.qinfbl /
                      blCompressibility.qinfbl -
                  1.0)) /
      (1.0 - blCompressibility.tkbl * (uei / blCompressibility.qinfbl) *
                 (uei / blCompressibility.qinfbl));
  data.param.uz_ms =
      (data.param.uz * (uei / blCompressibility.qinfbl) *
           (uei / blCompressibility.qinfbl) -
       uei) *
      blCompressibility.tkbl_ms /
      (1.0 - blCompressibility.tkbl * (uei / blCompressibility.qinfbl) *
                 (uei / blCompressibility.qinfbl));
  return data;
}

bool BoundaryLayerWorkflow::blsys() {
  blData &previous = state.previous();
  blData &current = state.current();

  SkinFrictionCoefficients skinFriction = blmid(flowRegime);
  current = boundaryLayerVariablesSolver.solve(current, flowRegime);

  if (flowRegime == FlowRegimeEnum::Similarity) {
    state.stepbl();
  }

  if (flowRegime == FlowRegimeEnum::Transition) {
    transitionSolver.trdif();
  } else {
    blc = blDiffSolver.solve(flowRegime, state, skinFriction,
                             blTransition.amcrit);
  }

  if (flowRegime == FlowRegimeEnum::Similarity) {
    blc.a2 += blc.a1;
    blc.a1.setZero();
  }

  for (int k = 0; k < 4; ++k) {
    double res_u1 = blc.a1(k, 3);
    double res_u2 = blc.a2(k, 3);
    double res_ms = blc.d_msq[k];

    blc.a1(k, 3) *= previous.param.uz_uei;
    blc.a2(k, 3) *= current.param.uz_uei;
    blc.d_msq[k] =
        res_u1 * previous.param.uz_ms + res_u2 * current.param.uz_ms + res_ms;
  }

  return true;
}

bool BoundaryLayerWorkflow::tesys(
    const BoundaryLayerSideProfiles &top_profiles,
    const BoundaryLayerSideProfiles &bottom_profiles, const Edge &edge) {
  blc.clear();

  state.station2 =
      boundaryLayerVariablesSolver.solve(state.station2, FlowRegimeEnum::Wake);

  const int top_te = lattice.top.trailingEdgeIndex;
  const int bottom_te = lattice.bottom.trailingEdgeIndex;
  const double tte = top_profiles.momentumThickness[top_te] +
                     bottom_profiles.momentumThickness[bottom_te];
  const double dte = top_profiles.displacementThickness[top_te] +
                     bottom_profiles.displacementThickness[bottom_te] +
                     edge.ante;
  const double cte = (top_profiles.skinFrictionCoeff[top_te] *
                          top_profiles.momentumThickness[top_te] +
                      bottom_profiles.skinFrictionCoeff[bottom_te] *
                          bottom_profiles.momentumThickness[bottom_te]) /
                     tte;

  blc.a1(0, 0) = -1.0;
  blc.a2(0, 0) = 1.0;
  blc.rhs[0] = cte - state.station2.param.sz;

  blc.a1(1, 1) = -1.0;
  blc.a2(1, 1) = 1.0;
  blc.rhs[1] = tte - state.station2.param.tz;

  blc.a1(2, 2) = -1.0;
  blc.a2(2, 2) = 1.0;
  blc.rhs[2] = dte - state.station2.param.dz - state.station2.param.dwz;

  return true;
}

void BoundaryLayerWorkflow::applySetblOutput(SetblOutputView &output) {
  blCompressibility = output.blCompressibility;
  blReynolds = output.blReynolds;
  lattice.top.profiles = std::move(output.profiles.top);
  lattice.bottom.profiles = std::move(output.profiles.bottom);
  flowRegime = output.flowRegime;
  blTransition = output.blTransition;
}

void BoundaryLayerWorkflow::checkTransitionIfNeeded(int side, int stationIndex,
                                                    bool skipCheck,
                                                    int laminarAdvance,
                                                    double &ami) {
  if (skipCheck || flowRegime == FlowRegimeEnum::Turbulent ||
      flowRegime == FlowRegimeEnum::Wake) {
    return;
  }

  transitionSolver.trchek();
  ami = state.station2.param.amplz;
  if (flowRegime == FlowRegimeEnum::Transition) {
    lattice.get(side).profiles.transitionIndex = stationIndex;
  } else {
    lattice.get(side).profiles.transitionIndex = stationIndex + laminarAdvance;
  }
}

SetblOutputView BoundaryLayerWorkflow::setbl(
    SidePairRef<const BoundaryLayerSideProfiles> profiles,
    const FlowState &analysis_state, const AeroCoefficients &aero_coeffs,
    double acrit, const Foil &foil, const StagnationResult &stagnation,
    const Eigen::MatrixXd &dij, bool bl_initialized) {
  return runBoundaryLayerSetbl(*this, profiles, analysis_state, aero_coeffs,
                               acrit, foil, stagnation, dij, bl_initialized);
}

SidePair<Eigen::VectorXd>
BoundaryLayerWorkflow::ueset(const Eigen::MatrixXd &dij) const {
  //---------------------------------------------------------
  //     sets ue from inviscid ue plus all source influence
  //---------------------------------------------------------
  SidePair<Eigen::VectorXd> edge_velocity;
  edge_velocity.top = lattice.top.profiles.edgeVelocity;
  edge_velocity.bottom = lattice.bottom.profiles.edgeVelocity;
  auto computeInducedVelocity = [&](int side, int station) {
    double dui = 0.0;
    const double side_panel_factor =
        lattice.get(side).panelInfluenceFactor[station];
    const int side_panel_index = lattice.get(side).stationToPanel[station];

    for (int js = 1; js <= 2; ++js) {
      for (int jbl = 0; jbl < lattice.get(js).stationCount - 1; ++jbl) {
        const double ue_m =
            -side_panel_factor * lattice.get(js).panelInfluenceFactor[jbl] *
            dij(side_panel_index, lattice.get(js).stationToPanel[jbl]);
        dui += ue_m * lattice.get(js).profiles.massFlux[jbl];
      }
    }
    return dui;
  };

  for (int is = 1; is <= 2; is++) {
    for (int ibl = 0; ibl < lattice.get(is).stationCount - 1; ++ibl) {
      const double dui = computeInducedVelocity(is, ibl);
      edge_velocity.get(is)[ibl] =
          lattice.get(is).inviscidEdgeVelocityMatrix(0, ibl) + dui;
    }
  }
  return edge_velocity;
}

/** -----------------------------------------------------
 * 	   sets forced-transition bl coordinate locations.
 * ----------------------------------------------------- */
double BoundaryLayerWorkflow::xifset(const Foil &foil,
                                     const StagnationResult &stagnation,
                                     int is) const {
  std::stringstream ss;
  VectorXd w1 = VectorXd::Zero(foil.foil_shape.n);
  double str;

  if (lattice.get(is).transitionLocation >= 1.0) {
    return lattice.get(is)
        .arcLengthCoordinates[lattice.get(is).trailingEdgeIndex];
  }

  Vector2d point_chord = foil.edge.point_te - foil.edge.point_le;

  //---- calculate chord-based x/c, y/c
  for (int i = 0; i < foil.foil_shape.n; i++) {
    w1[i] = (foil.foil_shape.points.col(i) - foil.edge.point_le)
                .dot(point_chord.normalized());
  }

  VectorXd w3 = spline::splind(w1, foil.foil_shape.spline_length);
  if (is == 1) {
    str = foil.edge.sle + (foil.foil_shape.spline_length[0] - foil.edge.sle) *
                              lattice.top.transitionLocation;
  } else {
    str =
        foil.edge.sle +
        (foil.foil_shape.spline_length[foil.foil_shape.n - 1] - foil.edge.sle) *
            lattice.bottom.transitionLocation;
  }
  str = spline::sinvrt(str, lattice.get(is).transitionLocation, w1, w3,
                       foil.foil_shape.spline_length, foil.foil_shape.n);
  double xiforc = std::min(
      (str - stagnation.sst),
      lattice.get(is).arcLengthCoordinates[lattice.get(is).trailingEdgeIndex]);
  if (xiforc < 0.0) {
    std::stringstream ss;
    ss << " ***  stagnation point is past trip on side " << is << "\n";
    Logger::instance().write(ss.str());
    return lattice.get(is)
        .arcLengthCoordinates[lattice.get(is).trailingEdgeIndex];
  }

  return xiforc;
}
