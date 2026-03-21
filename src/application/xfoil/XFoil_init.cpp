#include "application/xfoil/XFoil.h"
#include "solver/inviscid/InviscidSolver.hpp"
#include <algorithm>
#include <limits>
#include <numbers>

// Initialization and global state related member functions split from XFoil.cpp

bool XFoil::initialize() {

  // allocate arrays and clear containers
  initializeDataStructures();

  // reset numerical and physical variables
  resetVariables();

  //---- drop tolerance for bl system solver
  setVAccel(0.01);

  //---- set minf, reinf, based on current cl-dependence
  state_.operatingPointCoupling.machPerLift =
      getActualMach(1.0, analysis_state_.machType);
  state_.operatingPointCoupling.reynoldsPerLift =
      getActualReynolds(1.0, analysis_state_.reynoldsType);

  return true;
}

void XFoil::initializeDataStructures() {
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

void XFoil::resetFlags() {
  invalidateWakeGeometry();
  invalidatePanelMap();
  invalidateConvergedSolution();
  analysis_state_.viscous = false;
  analysis_state_.controlByAlpha = false;
  foil.edge.sharp = false;
  boundaryLayer.setFlowRegime(FlowRegimeEnum::Laminar);
}

void XFoil::resetVariables() {
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

double XFoil::getActualMach(double cls, MachType mach_control) {
  FlowState &state = analysis_state_;
  const double cla = std::max(cls, 0.000001);
  switch (mach_control) {
  case MachType::CONSTANT: {
    state.currentMach = state.referenceMach;
    return 0.0;
  }
  case MachType::FIXED_LIFT: {
    state.currentMach = state.referenceMach / std::sqrt(cla);
    return -0.5 * state.currentMach / cla;
  }
  case MachType::FIXED_LIFT_AND_DYNAMIC_PRESSURE: {
    state.currentMach = state.referenceMach;
    return 0.0;
  }
  default:
    return 0;
  }
}

double XFoil::getActualReynolds(double cls, ReynoldsType reynolds_control) {
  FlowState &state = analysis_state_;
  const double cla = std::max(cls, 0.000001);
  switch (reynolds_control) {
  case ReynoldsType::CONSTANT: {
    state.currentRe = state.referenceRe;
    return 0.0;
  }
  case ReynoldsType::FIXED_LIFT: {
    state.currentRe = state.referenceRe / std::sqrt(cla);
    return -0.5 * state.currentRe / cla;
  }
  case ReynoldsType::FIXED_LIFT_AND_DYNAMIC_PRESSURE: {
    state.currentRe = state.referenceRe / cla;
    return -state.currentRe / cla;
  }
  default:
    return 0;
  }
}

bool XFoil::setMach() {
  state_.operatingPointCoupling.machPerLift =
      getActualMach(1.0, analysis_state_.machType);
  state_.operatingPointCoupling.reynoldsPerLift =
      getActualReynolds(1.0, analysis_state_.reynoldsType);
  const int point_count = foil.foil_shape.n;
  cpi =
      InviscidSolver::cpcalc(point_count,
                             state_.inviscid.qinvMatrix.row(0).transpose(),
                             analysis_state_.qinf, analysis_state_.currentMach);
  if (analysis_state_.viscous) {
    cpv = InviscidSolver::cpcalc(point_count + foil.wake_shape.n,
                                 state_.viscous.qvis,
                                 analysis_state_.qinf,
                                 analysis_state_.currentMach);
  }
  const auto cl_result = clcalc(cmref);
  applyClComputation(cl_result);
  aero_coeffs_.cd = cdcalc();
  invalidateConvergedSolution();
  return true;
}
