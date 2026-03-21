#include "application/xfoil/XFoilAnalysis.hpp"

#include "Eigen/Core"
#include <cmath>
#include <limits>
#include <numbers>

#include "model/coefficient/xfoil_postprocess.hpp"
#include "solver/inviscid/InviscidSolver.hpp"
#include "solver/xfoil/xfoil_flowfield.hpp"

using Eigen::Matrix2Xd;

XFoilAnalysis::XFoilAnalysis(FlowState &analysis_state, XFoilResult &result,
                             double &acrit, XFoilWorkspace &workspace)
    : analysis_state_(analysis_state), result_(result),
      aero_coeffs_(result.aeroCoefficients), acrit_(acrit),
      workspace_(workspace), state_(workspace.state), foil(workspace.foil),
      boundaryLayer(workspace.boundaryLayer), cpi_(result.inviscidCp),
      cpv_(result.viscousCp) {
  state_.viscous.convergedMach = 0.0;
  acrit_ = 9.0;
  boundaryLayer.setTransitionLocations(1.0, 1.0);
  analysis_state_.machType = MachType::CONSTANT;
  analysis_state_.referenceMach = 0.0;
  setVAccel(0.01);
  analysis_state_.reynoldsType = ReynoldsType::CONSTANT;
  analysis_state_.referenceRe = 0.0;
  publishTransitionLocations();
}

XFoilAnalysis::CompressibilityParams
XFoilAnalysis::buildCompressibilityParams() const {
  return xfoil_postprocess::buildCompressibilityParams(analysis_state_);
}

XFoilAnalysis::PressureCoefficientResult
XFoilAnalysis::computePressureCoefficient(
    double tangential_velocity, double velocity_derivative,
    const CompressibilityParams &params) const {
  return xfoil_postprocess::computePressureCoefficient(
      tangential_velocity, velocity_derivative, analysis_state_.qinf, params);
}

bool XFoilAnalysis::isBLInitialized() const {
  if (!hasPanelMap()) {
    return false;
  }
  const auto top = boundaryLayer.readSideModel(1);
  const auto bottom = boundaryLayer.readSideModel(2);
  if (!top.hasStations || !bottom.hasStations) {
    return false;
  }
  return top.hasFiniteThickness && bottom.hasFiniteThickness;
}

void XFoilAnalysis::setBLInitialized(bool bInitialized) {
  if (bInitialized) {
    return;
  }
  boundaryLayer.zeroProfiles();
  invalidateConvergedSolution();
}

void XFoilAnalysis::invalidateConvergedSolution() {
  state_.viscous.convergedAlpha = std::numeric_limits<double>::quiet_NaN();
  state_.viscous.convergedMach = std::numeric_limits<double>::quiet_NaN();
  result_.converged = false;
}

void XFoilAnalysis::invalidateWakeGeometry() {
  const int point_count = foil.foil_shape.n;
  const int wake_point_count = foil.wake_shape.n;
  if (const int total_nodes = point_count + wake_point_count;
      wake_point_count > 0 && state_.inviscid.cache.dij.rows() >= total_nodes &&
      state_.inviscid.cache.dij.cols() >= total_nodes) {
    state_.inviscid.cache.dij
        .block(point_count, point_count, wake_point_count, wake_point_count)
        .setConstant(std::numeric_limits<double>::quiet_NaN());
  }
  foil.wake_shape.points.resize(0, 0);
  foil.wake_shape.normal_vector.resize(0, 0);
  foil.wake_shape.spline_length.resize(0);
  foil.wake_shape.angle_panel.resize(0);
}

void XFoilAnalysis::invalidatePanelMap() { boundaryLayer.clearPanelMap(); }

bool XFoilAnalysis::hasPanelMap() const {
  const int point_count = foil.foil_shape.n;
  const int total_nodes = point_count + foil.wake_shape.n;
  return boundaryLayer.hasValidPanelMap(total_nodes);
}

bool XFoilAnalysis::hasAirfoilInfluenceMatrix() const {
  const int point_count = foil.foil_shape.n;
  if (const int total_nodes = point_count + foil.wake_shape.n;
      state_.inviscid.cache.dij.rows() < total_nodes ||
      state_.inviscid.cache.dij.cols() < total_nodes) {
    return false;
  }
  const auto block = state_.inviscid.cache.dij.block(0, 0, point_count, point_count);
  return block.allFinite() && block.cwiseAbs().maxCoeff() > 0.0;
}

bool XFoilAnalysis::hasWakeInfluenceMatrix() const {
  const int point_count = foil.foil_shape.n;
  const int wake_point_count = foil.wake_shape.n;
  if (const int total_nodes = point_count + wake_point_count;
      state_.inviscid.cache.dij.rows() < total_nodes ||
      state_.inviscid.cache.dij.cols() < total_nodes) {
    return false;
  }
  const auto block = state_.inviscid.cache.dij.block(
      point_count, point_count, wake_point_count, wake_point_count);
  return block.allFinite() && block.cwiseAbs().maxCoeff() > 0.0;
}

bool XFoilAnalysis::hasConvergedSolution() const {
  if (!std::isfinite(state_.viscous.convergedAlpha) ||
      !std::isfinite(state_.viscous.convergedMach)) {
    return false;
  }
  const double alpha_tol = 1.0e-12;
  if (const double mach_tol = 1.0e-12;
      std::fabs(analysis_state_.alpha - state_.viscous.convergedAlpha) >
              alpha_tol ||
      std::fabs(analysis_state_.currentMach - state_.viscous.convergedMach) >
              mach_tol) {
    return false;
  }

  const int total_nodes_with_wake = foil.foil_shape.n + foil.wake_shape.n;
  return state_.viscous.qvis.size() >= total_nodes_with_wake &&
         state_.viscous.qvis.head(total_nodes_with_wake).allFinite();
}

void XFoilAnalysis::publishPressureCoefficients(const ViscalEndResult &result) {
  cpi_ = result.inviscidCp;
  cpv_ = result.viscousCp;
}

void XFoilAnalysis::publishTransitionLocations() {
  result_.transitionLocations.top = boundaryLayer.transitionLocation(1);
  result_.transitionLocations.bottom = boundaryLayer.transitionLocation(2);
}

double XFoilAnalysis::getXcp() const { return aero_coeffs_.xcp; }

double XFoilAnalysis::cdcalc() const {
  return xfoil_postprocess::computeCd(analysis_state_, boundaryLayer,
                                      isBLInitialized());
}

XFoilAnalysis::ClComputation XFoilAnalysis::clcalc(Vector2d ref) const {
  return xfoil_postprocess::computeCl(foil, analysis_state_, state_.inviscid,
                                      ref);
}

void XFoilAnalysis::applyClComputation(const ClComputation &result) {
  xfoil_postprocess::applyClComputation(aero_coeffs_, result);
}

XFoilAnalysis::Matrix2Xd XFoilAnalysis::gamqv() const {
  return xfoil_flowfield::buildSurfaceVortex(foil, state_);
}

FoilAerodynamicCache XFoilAnalysis::ggcalc() {
  return xfoil_flowfield::buildUnitVorticityDistributions(foil, analysis_state_,
                                                          state_.inviscid);
}

bool XFoilAnalysis::specal() { return InviscidSolver::specal(*this); }

bool XFoilAnalysis::speccl() { return InviscidSolver::speccl(*this); }

bool XFoilAnalysis::abcopy(Matrix2Xd copyFrom) {
  constexpr double kPointMergeTolerance = 1.0e-14;
  int point_count = static_cast<int>(copyFrom.cols());

  int r = 1;
  while (r < point_count) {
    const double delta_norm = (copyFrom.col(r - 1) - copyFrom.col(r)).norm();
    if (delta_norm <= kPointMergeTolerance) {
      for (int j = r; j < point_count - 1; j++) {
        copyFrom.col(j) = copyFrom.col(j + 1);
      }
      point_count -= 1;
    } else {
      r++;
    }
  }

  const int wake_point_count = point_count / 8 + 2;
  foil.wake_shape.n = wake_point_count;
  Matrix2Xd foil_points = Matrix2Xd::Zero(2, point_count + wake_point_count);
  foil_points.leftCols(point_count) = copyFrom.leftCols(point_count);

  foil.foil_shape.n = point_count;
  initialize();

  foil = Foil(foil_points, point_count);
  foil.wake_shape.n = wake_point_count;
  updateTrailingEdgeState();

  invalidateWakeGeometry();
  invalidatePanelMap();
  setBLInitialized(false);
  invalidateConvergedSolution();

  return true;
}

void XFoilAnalysis::updateTrailingEdgeState() {
  xfoil_flowfield::updateTrailingEdgeState(foil, state_.inviscid);
}

bool XFoilAnalysis::initialize() {
  initializeDataStructures();
  resetVariables();
  setVAccel(0.01);
  state_.operatingPointCoupling.machPerLift =
      getActualMach(1.0, analysis_state_.machType);
  state_.operatingPointCoupling.reynoldsPerLift =
      getActualReynolds(1.0, analysis_state_.reynoldsType);
  return true;
}

void XFoilAnalysis::initializeDataStructures() {
  const int point_count = foil.foil_shape.n;
  const int wake_nodes = foil.wake_shape.n;
  const int total_nodes_with_wake = point_count + wake_nodes;
  const int surface_buffer_nodes = point_count + 6;
  const int bl_node_count = point_count + wake_nodes;

  boundaryLayer.initializeLattices(bl_node_count);

  state_.inviscid.cache.bij =
      MatrixXd::Zero(point_count + 1, total_nodes_with_wake);
  state_.inviscid.cache.dij =
      MatrixXd::Zero(total_nodes_with_wake, total_nodes_with_wake);
  state_.inviscid.cache.gamu = Matrix2Xd::Zero(2, point_count + 1);
  state_.inviscid.surfaceVortex = Matrix2Xd::Zero(2, point_count);
  state_.init.qf0.assign(surface_buffer_nodes + 1, 0.0);
  state_.init.qf1.assign(surface_buffer_nodes + 1, 0.0);
  state_.init.qf2.assign(surface_buffer_nodes + 1, 0.0);
  state_.init.qf3.assign(surface_buffer_nodes + 1, 0.0);
  state_.inviscid.qinvu = Matrix2Xd::Zero(2, total_nodes_with_wake);
  state_.inviscid.qinvMatrix = Matrix2Xd::Zero(2, total_nodes_with_wake);
  state_.viscous.qvis = VectorXd::Zero(total_nodes_with_wake);

  boundaryLayer.initializeWakeGap(wake_nodes);
  boundaryLayer.systemCoefficients().clear();

  state_.inviscid.qgamm = VectorXd::Zero(point_count);
}

void XFoilAnalysis::resetFlags() {
  invalidateWakeGeometry();
  invalidatePanelMap();
  invalidateConvergedSolution();
  analysis_state_.viscous = false;
  analysis_state_.controlByAlpha = false;
  foil.edge.sharp = false;
  boundaryLayer.setFlowRegime(FlowRegimeEnum::Laminar);
}

void XFoilAnalysis::resetVariables() {
  boundaryLayer.resetForReinitialization();

  analysis_state_.qinf = 1.0;
  aero_coeffs_.cl = 0.0;
  aero_coeffs_.cm = 0.0;
  aero_coeffs_.cd = 0.0;
  aero_coeffs_.xcp = 0.0;
  state_.inviscid.sigmaTe = state_.inviscid.gammaTe = 0.0;
  state_.viscous.convergedAlpha = std::numeric_limits<double>::quiet_NaN();
  resetFlags();
  state_.init.amax = 0.0;
  analysis_state_.alpha = 0.0;
  analysis_state_.clspec = 0.0;
  foil.edge.ante = 0.0;
  foil.edge.aste = 0.0;
  foil.edge.dste = 0.0;
  analysis_state_.currentMach = 0.0;
  analysis_state_.currentRe = 0.0;
  state_.operatingPointCoupling.machPerLift = 0.0;
  state_.operatingPointCoupling.reynoldsPerLift = 0.0;
  aero_coeffs_.cl_alf = 0.0;
  aero_coeffs_.cl_msq = 0.0;
  state_.viscous.stagnation = StagnationResult{};
}

double XFoilAnalysis::getActualMach(double cls, MachType mach_control) {
  FlowState &state = analysis_state_;
  const double cla = std::max(cls, 0.000001);
  switch (mach_control) {
  case MachType::CONSTANT:
    state.currentMach = state.referenceMach;
    return 0.0;
  case MachType::FIXED_LIFT:
    state.currentMach = state.referenceMach / std::sqrt(cla);
    return -0.5 * state.currentMach / cla;
  case MachType::FIXED_LIFT_AND_DYNAMIC_PRESSURE:
    state.currentMach = state.referenceMach;
    return 0.0;
  default:
    return 0;
  }
}

double XFoilAnalysis::getActualReynolds(double cls,
                                        ReynoldsType reynolds_control) {
  FlowState &state = analysis_state_;
  const double cla = std::max(cls, 0.000001);
  switch (reynolds_control) {
  case ReynoldsType::CONSTANT:
    state.currentRe = state.referenceRe;
    return 0.0;
  case ReynoldsType::FIXED_LIFT:
    state.currentRe = state.referenceRe / std::sqrt(cla);
    return -0.5 * state.currentRe / cla;
  case ReynoldsType::FIXED_LIFT_AND_DYNAMIC_PRESSURE:
    state.currentRe = state.referenceRe / cla;
    return -state.currentRe / cla;
  default:
    return 0;
  }
}

bool XFoilAnalysis::setMach() {
  state_.operatingPointCoupling.machPerLift =
      getActualMach(1.0, analysis_state_.machType);
  state_.operatingPointCoupling.reynoldsPerLift =
      getActualReynolds(1.0, analysis_state_.reynoldsType);
  const int point_count = foil.foil_shape.n;
  cpi_ = InviscidSolver::cpcalc(point_count,
                                state_.inviscid.qinvMatrix.row(0).transpose(),
                                analysis_state_.qinf,
                                analysis_state_.currentMach);
  if (analysis_state_.viscous) {
    cpv_ = InviscidSolver::cpcalc(point_count + foil.wake_shape.n,
                                  state_.viscous.qvis, analysis_state_.qinf,
                                  analysis_state_.currentMach);
  }
  const auto cl_result = clcalc(cmref);
  applyClComputation(cl_result);
  aero_coeffs_.cd = cdcalc();
  invalidateConvergedSolution();
  return true;
}

double XFoilAnalysis::VAccel() { return vaccel_; }

void XFoilAnalysis::setVAccel(double accel) { vaccel_ = accel; }
