#include "BoundaryLayer.hpp"
#include "simulation/march/march.hpp"

#include <algorithm>
#include <array>
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

void BoundaryLayerWorkflow::runTransitionCheckForMrchue(int side, int stationIndex,
                                                        double &ami,
                                                        double &cti) {
  transitionSolver.trchek();
  ami = state.station2.param.amplz;
  if (flowRegime == FlowRegimeEnum::Transition) {
    lattice.get(side).profiles.transitionIndex = stationIndex;
    if (cti <= 0.0) {
      cti = 0.03;
      state.station2.param.sz = cti;
    }
  } else {
    lattice.get(side).profiles.transitionIndex = stationIndex + 2;
  }
}

bool BoundaryLayerWorkflow::solveTeSystemForCurrentProfiles(const Edge &edge) {
  return tesys(lattice.top.profiles, lattice.bottom.profiles, edge);
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

double BoundaryLayerWorkflow::calcHtarg(int ibl, int is, bool wake) {
  if (ibl < lattice.get(is).profiles.transitionIndex) {
    return state.station1.hkz.scalar +
           0.03 * (state.station2.param.xz - state.station1.param.xz) /
               state.station1.param.tz;
  }

  if (ibl == lattice.get(is).profiles.transitionIndex) {
    return state.station1.hkz.scalar +
           (0.03 * (xt.scalar - state.station1.param.xz) -
            0.15 * (state.station2.param.xz - xt.scalar)) /
               state.station1.param.tz;
  }

  if (wake) {
    const double cst = 0.03 *
                       (state.station2.param.xz - state.station1.param.xz) /
                       state.station1.param.tz;
    auto euler = [](double hk2, double hk1, double cst_local) {
      return hk2 - (hk2 + cst_local * pow(hk2 - 1, 3) - hk1) /
                       (1 + 3 * cst_local * pow(hk2 - 1, 2));
    };
    state.station2.hkz.scalar = state.station1.hkz.scalar;
    for (int i = 0; i < 3; i++) {
      state.station2.hkz.scalar =
          euler(state.station2.hkz.scalar, state.station1.hkz.scalar, cst);
    }
    return state.station2.hkz.scalar;
  }

  return state.station1.hkz.scalar -
         0.15 * (state.station2.param.xz - state.station1.param.xz) /
             state.station1.param.tz;
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

BoundaryLayerWorkflow::BlReferenceParams
BoundaryLayerWorkflow::computeBlReferenceParams(
    const FlowState &analysis_state, const AeroCoefficients &aero_coeffs,
    double acrit) const {
  BlReferenceParams params;
  double clmr = 0.0;
  double ma_clmr = 0.0;
  double herat = 0.0;
  double herat_ms = 0.0;
  //---- set the cl used to define mach, reynolds numbers
  if (analysis_state.controlByAlpha)
    clmr = aero_coeffs.cl;
  else
    clmr = analysis_state.clspec;

  //---- set current minf(cl)
  const double cla = std::max(clmr, 0.000001);
  switch (analysis_state.machType) {
  case FlowState::MachType::CONSTANT:
    params.currentMach = analysis_state.referenceMach;
    ma_clmr = 0.0;
    break;
  case FlowState::MachType::FIXED_LIFT:
    params.currentMach = analysis_state.referenceMach / std::sqrt(cla);
    ma_clmr = -0.5 * params.currentMach / cla;
    break;
  case FlowState::MachType::FIXED_LIFT_AND_DYNAMIC_PRESSURE:
    params.currentMach = analysis_state.referenceMach;
    ma_clmr = 0.0;
    break;
  default:
    params.currentMach = analysis_state.currentMach;
    ma_clmr = 0.0;
    break;
  }

  switch (analysis_state.reynoldsType) {
  case FlowState::ReynoldsType::CONSTANT:
    params.currentRe = analysis_state.referenceRe;
    params.re_clmr = 0.0;
    break;
  case FlowState::ReynoldsType::FIXED_LIFT:
    params.currentRe = analysis_state.referenceRe / std::sqrt(cla);
    params.re_clmr = -0.5 * params.currentRe / cla;
    break;
  case FlowState::ReynoldsType::FIXED_LIFT_AND_DYNAMIC_PRESSURE:
    params.currentRe = analysis_state.referenceRe / cla;
    params.re_clmr = -params.currentRe / cla;
    break;
  default:
    params.currentRe = analysis_state.currentRe;
    params.re_clmr = 0.0;
    break;
  }

  params.msq_clmr = 2.0 * params.currentMach * ma_clmr;

  //---- set compressibility parameter tklam and derivative tk_msq
  const double beta = std::sqrt(1.0 - params.currentMach * params.currentMach);
  const double beta_msq = -0.5 / beta;
  const double karmanTsienFactor =
      (params.currentMach / (1.0 + beta)) * (params.currentMach / (1.0 + beta));
  const double karmanTsienFactor_msq =
      1.0 / ((1.0 + beta) * (1.0 + beta)) -
      2.0 * karmanTsienFactor / (1.0 + beta) * beta_msq;
  //---- set gas constant (= cp/cv)
  params.blCompressibility.gm1bl = 1.4 - 1;

  //---- set parameters for compressibility correction
  params.blCompressibility.qinfbl = analysis_state.qinf;
  params.blCompressibility.tkbl = karmanTsienFactor;
  params.blCompressibility.tkbl_ms = karmanTsienFactor_msq;

  //---- stagnation density and 1/enthalpy
  params.blCompressibility.rstbl =
      pow((1.0 + 0.5 * params.blCompressibility.gm1bl * params.currentMach *
                     params.currentMach),
          (1.0 / params.blCompressibility.gm1bl));
  params.blCompressibility.rstbl_ms =
      0.5 * params.blCompressibility.rstbl /
      (1.0 + 0.5 * params.blCompressibility.gm1bl * params.currentMach *
                 params.currentMach);
  params.blCompressibility.hstinv =
      params.blCompressibility.gm1bl *
      MathUtil::pow(params.currentMach / params.blCompressibility.qinfbl, 2) /
      (1.0 + 0.5 * params.blCompressibility.gm1bl * params.currentMach *
                 params.currentMach);
  params.blCompressibility.hstinv_ms =
      params.blCompressibility.gm1bl *
          MathUtil::pow(1.0 / params.blCompressibility.qinfbl, 2) /
          (1.0 + 0.5 * params.blCompressibility.gm1bl * params.currentMach *
                     params.currentMach) -
      0.5 * params.blCompressibility.gm1bl * params.blCompressibility.hstinv /
          (1.0 + 0.5 * params.blCompressibility.gm1bl * params.currentMach *
                     params.currentMach);

  //---- set reynolds number based on freestream density, velocity, viscosity
  herat = 1.0 - 0.5 * params.blCompressibility.qinfbl *
                    params.blCompressibility.qinfbl *
                    params.blCompressibility.hstinv;
  herat_ms = -0.5 * params.blCompressibility.qinfbl *
             params.blCompressibility.qinfbl *
             params.blCompressibility.hstinv_ms;

  params.blReynolds.reybl = params.currentRe * sqrt(herat * herat * herat) *
                            (1.0 + BoundaryLayerWorkflow::kHvrat) /
                            (herat + BoundaryLayerWorkflow::kHvrat);
  params.blReynolds.reybl_re = sqrt(herat * herat * herat) *
                               (1.0 + BoundaryLayerWorkflow::kHvrat) /
                               (herat + BoundaryLayerWorkflow::kHvrat);
  params.blReynolds.reybl_ms =
      params.blReynolds.reybl *
      (1.5 / herat - 1.0 / (herat + BoundaryLayerWorkflow::kHvrat)) * herat_ms;

  params.amcrit = acrit;
  return params;
}

BoundaryLayerWorkflow::EdgeVelocitySensitivityResult
BoundaryLayerWorkflow::prepareEdgeVelocityAndSensitivities(
    SidePairRef<const BoundaryLayerSideProfiles> profiles,
    const Eigen::MatrixXd &dij, int nsys) const {
  EdgeVelocitySensitivityResult result;

  result.edgeVelocity = ueset(dij);
  result.outputEdgeVelocity.top = profiles.top.edgeVelocity;
  result.outputEdgeVelocity.bottom = profiles.bottom.edgeVelocity;

  result.jvte.top = lattice.top.stationToSystem[lattice.top.trailingEdgeIndex];
  result.jvte.bottom =
      lattice.bottom.stationToSystem[lattice.bottom.trailingEdgeIndex];

  result.dule.top =
      result.outputEdgeVelocity.top[0] - result.edgeVelocity.top[0];
  result.dule.bottom =
      result.outputEdgeVelocity.bottom[0] - result.edgeVelocity.bottom[0];

  //---- set le and te ue sensitivities wrt all m values
  const auto le_te_sensitivities = computeLeTeSensitivities(
      lattice.get(1).stationToPanel[0], lattice.get(2).stationToPanel[0],
      lattice.get(1).stationToPanel[lattice.top.trailingEdgeIndex],
      lattice.get(2).stationToPanel[lattice.bottom.trailingEdgeIndex], nsys,
      dij);
  result.ule_m = le_te_sensitivities.ule_m;
  result.ute_m = le_te_sensitivities.ute_m;

  result.ule_a.top = lattice.get(1).inviscidEdgeVelocityMatrix(1, 0);
  result.ule_a.bottom = lattice.get(2).inviscidEdgeVelocityMatrix(1, 0);
  return result;
}

void BoundaryLayerWorkflow::assembleBlJacobianForStation(
    int is, int iv, int nsys, const std::array<SetblStation, 2> &setblStations,
    const SetblSideData &setblSides, bool controlByAlpha, double re_clmr,
    double msq_clmr, SetblOutputView &output) {

  output.bl_newton_system.vb[iv] = blc.a1.block(0, 0, 3, 2);
  output.bl_newton_system.va[iv] = blc.a2.block(0, 0, 3, 2);

  Eigen::Matrix<double, 3, 4> A;
  A.col(0) = blc.a1.col(3).head<3>();
  A.col(1) = blc.a1.col(2).head<3>();
  A.col(2) = blc.a2.col(3).head<3>();
  A.col(3) = blc.a2.col(2).head<3>();

  Eigen::Matrix<double, 4, 2> B;
  B << setblStations[0].due, setblStations[0].u_a, setblStations[0].dds,
      setblStations[0].d_a, setblStations[1].due, setblStations[1].u_a,
      setblStations[1].dds, setblStations[1].d_a;
  const Eigen::Vector3d ax =
      (blc.a1.col(4) + blc.a2.col(4) + blc.d_xi).head<3>();
  const Eigen::RowVector2d xi =
      setblStations[0].xi_ule *
          Eigen::RowVector2d(setblSides.dule.get(1), setblSides.ule_a.get(1)) +
      setblStations[1].xi_ule *
          Eigen::RowVector2d(setblSides.dule.get(2), setblSides.ule_a.get(2));
  output.bl_newton_system.vdel[iv] = A * B + ax * xi;

  for (int jv = 1; jv < nsys; jv++) {
    const Eigen::Vector4d m(setblStations[0].u_m(jv), setblStations[0].d_m(jv),
                            setblStations[1].u_m(jv), setblStations[1].d_m(jv));
    const double xi_m = setblStations[0].xi_ule * setblSides.ule_m.get(1)(jv) +
                        setblStations[1].xi_ule * setblSides.ule_m.get(2)(jv);
    const Eigen::Vector3d vm = A * m + ax * xi_m;
    output.bl_newton_system.vm.at(0, jv, iv) = vm[0];
    output.bl_newton_system.vm.at(1, jv, iv) = vm[1];
    output.bl_newton_system.vm.at(2, jv, iv) = vm[2];
  }

  if (controlByAlpha) {
    output.bl_newton_system.vdel[iv].col(1).head<3>() =
        (blc.d_re.head(3)) * re_clmr + (blc.d_msq.head(3)) * msq_clmr;
  }
}

BoundaryLayerWorkflow::SimilarityStationCoefficients
BoundaryLayerWorkflow::resetSimilarityStationCoefficients(
    const VectorXd &u_m1, const VectorXd &d_m1) const {
  SimilarityStationCoefficients result;
  result.u_m1 = u_m1;
  result.d_m1 = d_m1;
  for (int js = 1; js <= 2; ++js) {
    for (int jbl = 0; jbl < lattice.get(js).stationCount - 1; ++jbl) {
      const int jv = lattice.get(js).stationToSystem[jbl];
      result.u_m1[jv] = 0.0;
      result.d_m1[jv] = 0.0;
    }
  }
  return result;
}

BoundaryLayerWorkflow::SideSweepInitResult
BoundaryLayerWorkflow::initializeSideSweepState(
    const Foil &foil, const StagnationResult &stagnation, int is) const {
  SideSweepInitResult result;
  result.u_a1 = 0.0;
  result.d_a1 = 0.0;
  result.due1 = 0.0;
  result.dds1 = 0.0;
  result.xiforc = xifset(foil, stagnation, is);
  return result;
}

BoundaryLayerWorkflow::StationPrimaryVars
BoundaryLayerWorkflow::loadStationPrimaryVars(int is, int ibl,
                                              bool stationIsWake,
                                              const SetblOutputView &output,
                                              double ami, double cti) const {
  StationPrimaryVars vars;
  vars.xsi = lattice.get(is).arcLengthCoordinates[ibl];

  if (ibl < output.profiles.get(is).transitionIndex)
    vars.ami = output.profiles.get(is).skinFrictionCoeff[ibl];
  else
    vars.cti = output.profiles.get(is).skinFrictionCoeff[ibl];
  if (ibl < output.profiles.get(is).transitionIndex) {
    vars.cti = cti;
  } else {
    vars.ami = ami;
  }

  vars.uei = output.profiles.get(is).edgeVelocity[ibl];
  vars.thi = output.profiles.get(is).momentumThickness[ibl];
  vars.mdi = output.profiles.get(is).massFlux[ibl];
  vars.dsi = vars.mdi / vars.uei;

  if (stationIsWake) {
    int iw = ibl - lattice.get(is).trailingEdgeIndex;
    vars.dswaki = wgap[iw - 1];
  } else {
    vars.dswaki = 0.0;
  }

  return vars;
}

BoundaryLayerWorkflow::StationUpdateResult
BoundaryLayerWorkflow::updateStationMatricesAndState(
    int is, int ibl, int iv, const StationPrimaryVars &vars,
    const SidePair<VectorXd> &usav, const SetblOutputView &output,
    const BoundaryLayerState &base_state, int system_size,
    const Eigen::MatrixXd &dij) {
  StationUpdateResult result;
  const double d2_m2 = 1.0 / vars.uei;
  const double d2_u2 = -vars.dsi / vars.uei;

  result.u_m2 = VectorXd::Zero(system_size);
  result.d_m2 = VectorXd::Zero(system_size);
  for (int js = 1; js <= 2; js++) {
    for (int jbl = 0; jbl < lattice.get(js).stationCount - 1; ++jbl) {
      int jv = lattice.get(js).stationToSystem[jbl];
      result.u_m2[jv] = -lattice.get(is).panelInfluenceFactor[ibl] *
                        lattice.get(js).panelInfluenceFactor[jbl] *
                        dij(lattice.get(is).stationToPanel[ibl],
                            lattice.get(js).stationToPanel[jbl]);
      result.d_m2[jv] = d2_u2 * result.u_m2[jv];
    }
  }
  result.d_m2[iv] = result.d_m2[iv] + d2_m2;

  result.u_a2 = lattice.get(is).inviscidEdgeVelocityMatrix(1, ibl);
  result.d_a2 = d2_u2 * result.u_a2;

  // "forced" changes from mismatch between edge velocities and usav
  result.due2 = output.profiles.get(is).edgeVelocity[ibl] - usav.get(is)[ibl];
  result.dds2 = d2_u2 * result.due2;

  result.state = base_state;
  result.state.current() =
      blprv(result.state.current(), vars.xsi, vars.ami, vars.cti, vars.thi,
            vars.dsi, vars.dswaki, vars.uei);
  blkin(result.state);
  return result;
}

void BoundaryLayerWorkflow::buildTransitionLog(
    bool stationIsTransitionCandidate, FlowRegimeEnum flowRegime) const {
  if (stationIsTransitionCandidate &&
      flowRegime != FlowRegimeEnum::Transition) {
    std::stringstream ss;
    ss << "setbl: xtr???  n1=" << state.station1.param.amplz
       << " n2=" << state.station2.param.amplz << ":\n";
    Logger::instance().write(ss.str());
  }
}

BoundaryLayerWorkflow::TeWakeUpdateResult
BoundaryLayerWorkflow::computeTeWakeCoefficients(
    int is, int ibl, const SidePair<VectorXd> &usav,
    const SidePair<VectorXd> &ute_m, const SidePair<int> &jvte,
    const VectorXd &d_m1_template, const SetblOutputView &output,
    const Edge &edge) const {
  TeWakeUpdateResult result;
  if (ibl != lattice.get(is).trailingEdgeIndex + 1) {
    return result;
  }
  result.isStartOfWake = true;

  result.coeffs.tte =
      output.profiles.get(1).momentumThickness[lattice.top.trailingEdgeIndex] +
      output.profiles.get(2)
          .momentumThickness[lattice.bottom.trailingEdgeIndex];
  result.coeffs.dte =
      output.profiles.get(1)
          .displacementThickness[lattice.top.trailingEdgeIndex] +
      output.profiles.get(2)
          .displacementThickness[lattice.bottom.trailingEdgeIndex] +
      edge.ante;
  result.coeffs.cte =
      (output.profiles.get(1).skinFrictionCoeff[lattice.top.trailingEdgeIndex] *
           output.profiles.get(1)
               .momentumThickness[lattice.top.trailingEdgeIndex] +
       output.profiles.get(2)
               .skinFrictionCoeff[lattice.bottom.trailingEdgeIndex] *
           output.profiles.get(2)
               .momentumThickness[lattice.bottom.trailingEdgeIndex]) /
      result.coeffs.tte;

  result.coeffs.tte_tte1 = 1.0;
  result.coeffs.tte_tte2 = 1.0;
  result.coeffs.dte_mte1 =
      1.0 / output.profiles.top.edgeVelocity[lattice.top.trailingEdgeIndex];
  result.coeffs.dte_ute1 =
      -output.profiles.get(1)
           .displacementThickness[lattice.top.trailingEdgeIndex] /
      output.profiles.top.edgeVelocity[lattice.top.trailingEdgeIndex];
  result.coeffs.dte_mte2 =
      1.0 /
      output.profiles.bottom.edgeVelocity[lattice.bottom.trailingEdgeIndex];
  result.coeffs.dte_ute2 =
      -output.profiles.get(2)
           .displacementThickness[lattice.bottom.trailingEdgeIndex] /
      output.profiles.bottom.edgeVelocity[lattice.bottom.trailingEdgeIndex];
  result.coeffs.cte_cte1 =
      output.profiles.get(1).momentumThickness[lattice.top.trailingEdgeIndex] /
      result.coeffs.tte;
  result.coeffs.cte_cte2 =
      output.profiles.get(2)
          .momentumThickness[lattice.bottom.trailingEdgeIndex] /
      result.coeffs.tte;
  result.coeffs.cte_tte1 =
      (output.profiles.get(1).skinFrictionCoeff[lattice.top.trailingEdgeIndex] -
       result.coeffs.cte) /
      result.coeffs.tte;
  result.coeffs.cte_tte2 =
      (output.profiles.get(2)
           .skinFrictionCoeff[lattice.bottom.trailingEdgeIndex] -
       result.coeffs.cte) /
      result.coeffs.tte;

  // Re-define d1 sensitivities wrt m since d1 depends on both te ds values.
  result.d_m1 = d_m1_template;
  for (int js = 1; js <= 2; js++) {
    for (int jbl = 0; jbl < lattice.get(js).stationCount - 1; ++jbl) {
      int jv = lattice.get(js).stationToSystem[jbl];
      result.d_m1[jv] = result.coeffs.dte_ute1 * ute_m.get(1)[jv] +
                        result.coeffs.dte_ute2 * ute_m.get(2)[jv];
    }
  }
  result.d_m1[jvte.get(1)] = result.d_m1[jvte.get(1)] + result.coeffs.dte_mte1;
  result.d_m1[jvte.get(2)] = result.d_m1[jvte.get(2)] + result.coeffs.dte_mte2;

  // "forced" changes from edge velocity mismatch
  result.due1 = 0.0;
  result.dds1 =
      result.coeffs.dte_ute1 *
          (output.profiles.top.edgeVelocity[lattice.top.trailingEdgeIndex] -
           usav.top[lattice.top.trailingEdgeIndex]) +
      result.coeffs.dte_ute2 *
          (output.profiles.bottom
               .edgeVelocity[lattice.bottom.trailingEdgeIndex] -
           usav.bottom[lattice.bottom.trailingEdgeIndex]);

  return result;
}

BoundaryLayerWorkflow::TeWakeJacobianAdjustments
BoundaryLayerWorkflow::computeTeWakeJacobianAdjustments(
    const TeWakeCoefficients &coeffs) const {
  TeWakeJacobianAdjustments result;
  result.vz[0][0] = blc.a1(0, 0) * coeffs.cte_cte1;
  result.vz[0][1] =
      blc.a1(0, 0) * coeffs.cte_tte1 + blc.a1(0, 1) * coeffs.tte_tte1;
  result.vb(0, 0) = blc.a1(0, 0) * coeffs.cte_cte2;
  result.vb(0, 1) =
      blc.a1(0, 0) * coeffs.cte_tte2 + blc.a1(0, 1) * coeffs.tte_tte2;

  result.vz[1][0] = blc.a1(1, 0) * coeffs.cte_cte1;
  result.vz[1][1] =
      blc.a1(1, 0) * coeffs.cte_tte1 + blc.a1(1, 1) * coeffs.tte_tte1;
  result.vb(1, 0) = blc.a1(1, 0) * coeffs.cte_cte2;
  result.vb(1, 1) =
      blc.a1(1, 0) * coeffs.cte_tte2 + blc.a1(1, 1) * coeffs.tte_tte2;

  result.vz[2][0] = blc.a1(2, 0) * coeffs.cte_cte1;
  result.vz[2][1] =
      blc.a1(2, 0) * coeffs.cte_tte1 + blc.a1(2, 1) * coeffs.tte_tte1;
  result.vb(2, 0) = blc.a1(2, 0) * coeffs.cte_cte2;
  result.vb(2, 1) =
      blc.a1(2, 0) * coeffs.cte_tte2 + blc.a1(2, 1) * coeffs.tte_tte2;
  return result;
}

BoundaryLayerWorkflow::StationArraysAdvanceResult
BoundaryLayerWorkflow::advanceStationArrays(const VectorXd &u_m2,
                                            const VectorXd &d_m2, double u_a2,
                                            double d_a2, double due2,
                                            double dds2) const {
  StationArraysAdvanceResult result;
  result.u_m1 = u_m2;
  result.d_m1 = d_m2;
  result.u_a1 = u_a2;
  result.d_a1 = d_a2;
  result.due1 = due2;
  result.dds1 = dds2;
  return result;
}

void BoundaryLayerWorkflow::initializeSetblReferenceParams(
    const FlowState &analysis_state, const AeroCoefficients &aero_coeffs,
    double acrit, SetblOutputView &output, double &re_clmr, double &msq_clmr,
    double &currentMach, double &currentRe) {
  auto reference_params =
      computeBlReferenceParams(analysis_state, aero_coeffs, acrit);
  currentMach = reference_params.currentMach;
  currentRe = reference_params.currentRe;
  re_clmr = reference_params.re_clmr;
  msq_clmr = reference_params.msq_clmr;
  output.blCompressibility = reference_params.blCompressibility;
  output.blReynolds = reference_params.blReynolds;
  output.blTransition.amcrit = reference_params.amcrit;
  blCompressibility = output.blCompressibility;
  blReynolds = output.blReynolds;
  blTransition.amcrit = output.blTransition.amcrit;
}

void BoundaryLayerWorkflow::initializeSetblSystemStorage(
    SetblOutputView &output, std::array<SetblStation, 2> &stations,
    SetblSideData &sideData) const {
  const auto zero = Eigen::Matrix<double, 3, 2>::Zero();
  output.bl_newton_system.vm.resize(nsys);
  output.bl_newton_system.va.resize(nsys, zero);
  output.bl_newton_system.vb.resize(nsys, zero);
  output.bl_newton_system.vdel.resize(nsys, zero);
  stations[0].resizeSystem(nsys);
  stations[1].resizeSystem(nsys);
  sideData.resizeSystem(nsys);
}

void BoundaryLayerWorkflow::initializeSetblProfiles(
    SetblOutputView &output) const {
  output.profiles.top = lattice.top.profiles;
  output.profiles.bottom = lattice.bottom.profiles;
}

void BoundaryLayerWorkflow::initializeSetblEdgeVelocityState(
    SidePairRef<const BoundaryLayerSideProfiles> profiles,
    const Eigen::MatrixXd &dij, SetblOutputView &output,
    SetblSideData &sideData) const {
  auto edge_result = prepareEdgeVelocityAndSensitivities(profiles, dij, nsys);
  sideData.usav = edge_result.edgeVelocity;
  sideData.jvte = edge_result.jvte;
  sideData.dule = edge_result.dule;
  sideData.ule_m = edge_result.ule_m;
  sideData.ute_m = edge_result.ute_m;
  sideData.ule_a = edge_result.ule_a;
  output.profiles.top.edgeVelocity = edge_result.outputEdgeVelocity.top;
  output.profiles.bottom.edgeVelocity = edge_result.outputEdgeVelocity.bottom;
}

void BoundaryLayerWorkflow::processSetblSide(
    int side, const Foil &foil, const StagnationResult &stagnation,
    bool controlByAlpha, const Eigen::MatrixXd &dij, SetblOutputView &output,
    std::array<SetblStation, 2> &stations, SetblSideData &sideData, double &cti,
    double &ami, double re_clmr, double msq_clmr) {
  auto similarity_coeffs =
      resetSimilarityStationCoefficients(stations[0].u_m, stations[0].d_m);
  stations[0].u_m = similarity_coeffs.u_m1;
  stations[0].d_m = similarity_coeffs.d_m1;

  auto sweep_init = initializeSideSweepState(foil, stagnation, side);
  stations[0].u_a = sweep_init.u_a1;
  stations[0].d_a = sweep_init.d_a1;
  stations[0].due = sweep_init.due1;
  stations[0].dds = sweep_init.dds1;
  output.blTransition.xiforc = sweep_init.xiforc;
  blTransition.xiforc = output.blTransition.xiforc;

  for (int station = 0; station < lattice.get(side).stationCount - 1;
       ++station) {
    const int iv = lattice.get(side).stationToSystem[station];
    const bool station_is_wake =
        (station > lattice.get(side).trailingEdgeIndex);
    const bool station_is_transition_candidate =
        (station == output.profiles.get(side).transitionIndex);

    output.flowRegime = determineRegimeForStation(side, station);
    flowRegime = output.flowRegime;

    auto vars = loadStationPrimaryVars(side, station, station_is_wake, output,
                                       ami, cti);
    ami = vars.ami;
    cti = vars.cti;

    auto station_update = updateStationMatricesAndState(
        side, station, iv, vars, sideData.usav, output, state,
        stations[1].u_m.size(), dij);
    stations[1].u_m = station_update.u_m2;
    stations[1].d_m = station_update.d_m2;
    stations[1].u_a = station_update.u_a2;
    stations[1].d_a = station_update.d_a2;
    stations[1].due = station_update.due2;
    stations[1].dds = station_update.dds2;
    state = station_update.state;

    if (station_is_transition_candidate) {
      transitionSolver.trchek();
      ami = state.station2.param.amplz;
    }
    buildTransitionLog(station_is_transition_candidate, output.flowRegime);

    auto te_update = computeTeWakeCoefficients(
        side, station, sideData.usav, sideData.ute_m, sideData.jvte,
        stations[0].d_m, output, foil.edge);
    if (te_update.isStartOfWake) {
      tesys(lattice.top.profiles, lattice.bottom.profiles, foil.edge);
      stations[0].d_m = te_update.d_m1;
      stations[0].due = te_update.due1;
      stations[0].dds = te_update.dds1;
    } else {
      blsys();
    }

    output.profiles.get(side).skinFrictionCoeffHistory[station] =
        state.station2.cqz.scalar;

    if (side == 1) {
      stations[0].xi_ule = stagnation.sst_go;
      stations[1].xi_ule = -stagnation.sst_gp;
    } else {
      stations[0].xi_ule = -stagnation.sst_go;
      stations[1].xi_ule = stagnation.sst_gp;
    }

    assembleBlJacobianForStation(side, iv, nsys, stations, sideData,
                                 controlByAlpha, re_clmr, msq_clmr, output);

    if (te_update.isStartOfWake) {
      auto te_jacobian = computeTeWakeJacobianAdjustments(te_update.coeffs);
      for (int row = 0; row < 3; ++row) {
        output.bl_newton_system.vz[row][0] = te_jacobian.vz[row][0];
        output.bl_newton_system.vz[row][1] = te_jacobian.vz[row][1];
      }
      output.bl_newton_system.vb[iv] = te_jacobian.vb;
    }

    if (output.flowRegime == FlowRegimeEnum::Transition) {
      output.profiles.get(side).transitionIndex = station;
      lattice.get(side).profiles.transitionIndex = station;
      output.flowRegime = FlowRegimeEnum::Turbulent;
      flowRegime = output.flowRegime;
    }

    if (station == lattice.get(side).trailingEdgeIndex) {
      output.flowRegime = FlowRegimeEnum::Wake;
      flowRegime = output.flowRegime;
      state.station2 = boundaryLayerVariablesSolver.solve(state.station2,
                                                          FlowRegimeEnum::Wake);
      blmid(FlowRegimeEnum::Wake);
    }

    auto advance =
        advanceStationArrays(stations[1].u_m, stations[1].d_m, stations[1].u_a,
                             stations[1].d_a, stations[1].due, stations[1].dds);
    stations[0].u_m = advance.u_m1;
    stations[0].d_m = advance.d_m1;
    stations[0].u_a = advance.u_a1;
    stations[0].d_a = advance.d_a1;
    stations[0].due = advance.due1;
    stations[0].dds = advance.dds1;
    state.stepbl();
  }
}

namespace {

struct SetblWorkingState {
  std::array<BoundaryLayerWorkflow::SetblStation, 2> stations{};
  BoundaryLayerWorkflow::SetblSideData sides;
  double cti = 0.0;
  double ami = 0.0;
  double re_clmr = 0.0;
  double msq_clmr = 0.0;
};

} // namespace

SetblOutputView BoundaryLayerWorkflow::setbl(
    SidePairRef<const BoundaryLayerSideProfiles> profiles,
    const FlowState &analysis_state, const AeroCoefficients &aero_coeffs,
    double acrit, const Foil &foil, const StagnationResult &stagnation,
    const Eigen::MatrixXd &dij, bool bl_initialized) {
  //-------------------------------------------------
  //	   sets up the bl newton system coefficients for the current bl
  // variables
  //     and the edge velocities received from setup. the local bl system
  //     coefficients are then incorporated into the global newton system.
  //-------------------------------------------------
  SetblOutputView output{};
  MarcherUe marcher_ue;
  MarcherDu marcher_du;
  SetblWorkingState state;
  initializeSetblSystemStorage(output, state.stations, state.sides);

  double current_mach = analysis_state.currentMach;
  double current_re = analysis_state.currentRe;
  initializeSetblReferenceParams(analysis_state, aero_coeffs, acrit, output,
                                 state.re_clmr, state.msq_clmr, current_mach,
                                 current_re);

  if (!bl_initialized) {
    Logger::instance().write("   Initializing bl ...\n");
    marcher_ue.mrchue(*this, foil, stagnation);
  }
  marcher_du.mrchdu(*this, foil, stagnation);
  initializeSetblProfiles(output);

  initializeSetblEdgeVelocityState(profiles, dij, output, state.sides);

  for (int side = 1; side <= 2; ++side) {
    processSetblSide(side, foil, stagnation, analysis_state.controlByAlpha, dij,
                     output, state.stations, state.sides, state.cti, state.ami,
                     state.re_clmr, state.msq_clmr);
  }

  return output;
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

BoundaryLayerWorkflow::LeTeSensitivities
BoundaryLayerWorkflow::computeLeTeSensitivities(
    int ile1, int ile2, int ite1, int ite2, int nsys,
    const Eigen::MatrixXd &dij) const {
  LeTeSensitivities sensitivities;
  sensitivities.ule_m.top = VectorXd::Zero(nsys);
  sensitivities.ule_m.bottom = VectorXd::Zero(nsys);
  sensitivities.ute_m.top = VectorXd::Zero(nsys);
  sensitivities.ute_m.bottom = VectorXd::Zero(nsys);
  for (int js = 1; js <= 2; ++js) {
    for (int jbl = 0; jbl < lattice.get(js).stationCount - 1; ++jbl) {
      const int j = lattice.get(js).stationToPanel[jbl];
      const int jv = lattice.get(js).stationToSystem[jbl];
      const double panelInfluenceFactor_js =
          lattice.get(js).panelInfluenceFactor[jbl];
      sensitivities.ule_m.top[jv] = -lattice.top.panelInfluenceFactor[0] *
                                    panelInfluenceFactor_js * dij(ile1, j);
      sensitivities.ule_m.bottom[jv] = -lattice.bottom.panelInfluenceFactor[0] *
                                       panelInfluenceFactor_js * dij(ile2, j);
      sensitivities.ute_m.top[jv] =
          -lattice.top.panelInfluenceFactor[lattice.top.trailingEdgeIndex] *
          panelInfluenceFactor_js * dij(ite1, j);
      sensitivities.ute_m.bottom[jv] =
          -lattice.bottom
               .panelInfluenceFactor[lattice.bottom.trailingEdgeIndex] *
          panelInfluenceFactor_js * dij(ite2, j);
    }
  }
  return sensitivities;
}
