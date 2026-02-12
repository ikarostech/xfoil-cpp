#include "XFoil.h"
#include <algorithm>
#include <numbers>
#include <unordered_map>
#include "simulation/InviscidSolver.hpp"

// Initialization and global state related member functions split from XFoil.cpp

namespace {
struct InitState {
  double amax = 0.0;
  std::vector<double> qf0;
  std::vector<double> qf1;
  std::vector<double> qf2;
  std::vector<double> qf3;
};

using InitStateRegistry = std::unordered_map<const XFoil*, InitState>;

InitStateRegistry& initStateRegistry() {
  static InitStateRegistry state;
  return state;
}

InitState& ensureInitState(const XFoil* xfoil) {
  return initStateRegistry()[xfoil];
}
}  // namespace

void ClearInitState(const XFoil& xfoil) {
  initStateRegistry().erase(&xfoil);
}

bool XFoil::initialize() {

  // allocate arrays and clear containers
  initializeDataStructures();

  // reset numerical and physical variables
  resetVariables();

  //---- drop tolerance for bl system solver
  setVAccel(0.01);

  //---- set minf, reinf, based on current cl-dependence
  minf_cl = getActualMach(1.0, analysis_state_.machType);
  reinf_cl = getActualReynolds(1.0, analysis_state_.reynoldsType);

  return true;
}

void XFoil::initializeDataStructures() {
  const int point_count = foil.foil_shape.n;
  const int wake_nodes = foil.wake_shape.n;
  const int total_nodes_with_wake = point_count + wake_nodes;
  const int surface_buffer_nodes = point_count + 6;
  const int bl_node_count = point_count + wake_nodes;
  const int bl_system_size = 2 * bl_node_count + 2;
  auto& cache = ensureInitState(this);

  boundaryLayerWorkflow.lattice.top = BoundaryLayerLattice(bl_node_count);
  boundaryLayerWorkflow.lattice.bottom = BoundaryLayerLattice(bl_node_count);

  aerodynamicCache.bij = MatrixXd::Zero(point_count + 1, total_nodes_with_wake);
  aerodynamicCache.dij =
      MatrixXd::Zero(total_nodes_with_wake, total_nodes_with_wake);
  aerodynamicCache.gamu = Matrix2Xd::Zero(2, point_count + 1);
  surface_vortex = Matrix2Xd::Zero(2, point_count);
  cache.qf0.assign(surface_buffer_nodes + 1, 0.0);
  cache.qf1.assign(surface_buffer_nodes + 1, 0.0);
  cache.qf2.assign(surface_buffer_nodes + 1, 0.0);
  cache.qf3.assign(surface_buffer_nodes + 1, 0.0);
  qinv_matrix = Matrix2Xd::Zero(2, total_nodes_with_wake);
  aerodynamicCache.qinvu = Matrix2Xd::Zero(2, total_nodes_with_wake);
  qvis = VectorXd::Zero(total_nodes_with_wake);

  boundaryLayerWorkflow.wgap = VectorXd::Zero(wake_nodes);
  /*
  bl_newton_system.va.resize(bl_system_size + 1, Matrix3x2d::Zero());
  bl_newton_system.vb.resize(bl_system_size + 1, Matrix3x2d::Zero());
  bl_newton_system.vdel.resize(bl_system_size + 1, Matrix3x2d::Zero());
  bl_newton_system.vm.resize(bl_system_size + 1);
  */
  boundaryLayerWorkflow.blc.clear();
  /*
  for (auto& row : bl_newton_system.vz) {
    row.fill(0.0);
  }
  */

  qgamm = VectorXd::Zero(point_count);
}

void XFoil::resetFlags() {
  lwake = lblini = lipan = false;
  ladij = lwdij = lvconv = false;
  analysis_state_.viscous = false;
  analysis_state_.controlByAlpha = false;
  foil.edge.sharp = false;
  boundaryLayerWorkflow.flowRegime = FlowRegimeEnum::Laminar;
}

void XFoil::resetVariables() {
  boundaryLayerWorkflow.state.station1 = blData{};
  boundaryLayerWorkflow.state.station2 = blData{};
  boundaryLayerWorkflow.lattice.top.profiles.transitionIndex = 0;
  boundaryLayerWorkflow.lattice.bottom.profiles.transitionIndex = 0;
  boundaryLayerWorkflow.lattice.top.stationCount = 0;
  boundaryLayerWorkflow.lattice.bottom.stationCount = 0;
  boundaryLayerWorkflow.lattice.top.trailingEdgeIndex = 0;
  boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex = 0;
  
  analysis_state_.qinf = 1.0;
  aero_coeffs_.cl = 0.0;
  aero_coeffs_.cm = 0.0;
  aero_coeffs_.cd = 0.0;
  aero_coeffs_.xcp = 0.0;
  sigte = gamte = 0.0;
  avisc = 0.0;
  resetFlags();
  boundaryLayerWorkflow.stagnationIndex = 0;
  boundaryLayerWorkflow.stagnationSst = 0.0;
  boundaryLayerWorkflow.blCompressibility.qinfbl =
      boundaryLayerWorkflow.blCompressibility.tkbl =
          boundaryLayerWorkflow.blCompressibility.tkbl_ms = 0.0;
  boundaryLayerWorkflow.blCompressibility.rstbl =
      boundaryLayerWorkflow.blCompressibility.rstbl_ms = 0.0;
  boundaryLayerWorkflow.blCompressibility.hstinv =
      boundaryLayerWorkflow.blCompressibility.hstinv_ms = 0.0;
  boundaryLayerWorkflow.blReynolds.reybl =
      boundaryLayerWorkflow.blReynolds.reybl_ms =
          boundaryLayerWorkflow.blReynolds.reybl_re = 0.0;
  boundaryLayerWorkflow.blCompressibility.gm1bl = 0.0;
  boundaryLayerWorkflow.blTransition.xiforc = 0.0;
  boundaryLayerWorkflow.blTransition.amcrit = 0.0;
  auto& cache = ensureInitState(this);
  cache.amax = 0.0;
  analysis_state_.alpha = 0.0;
  analysis_state_.clspec = 0.0;
  foil.edge.ante = 0.0;
  foil.edge.aste = 0.0;
  foil.edge.dste = 0.0;
  analysis_state_.currentMach = 0.0;
  analysis_state_.currentRe = 0.0;
  minf_cl = reinf_cl = 0.0;
  aero_coeffs_.cl_alf = 0.0;
  aero_coeffs_.cl_msq = 0.0;
  stagnation = StagnationResult{};
}

double XFoil::getActualMach(double cls, MachType mach_control) {
  FlowState& state = analysis_state_;
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
  FlowState& state = analysis_state_;
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
  minf_cl = getActualMach(1.0, analysis_state_.machType);
  reinf_cl = getActualReynolds(1.0, analysis_state_.reynoldsType);
  const int point_count = foil.foil_shape.n;
  cpi = InviscidSolver::cpcalc(point_count, qinv_matrix.row(0).transpose(),
               analysis_state_.qinf, analysis_state_.currentMach);
  if (analysis_state_.viscous) {
    cpv = InviscidSolver::cpcalc(point_count + foil.wake_shape.n, qvis, analysis_state_.qinf, analysis_state_.currentMach);
  }
  const auto cl_result = clcalc(cmref);
  applyClComputation(cl_result);
  aero_coeffs_.cd = cdcalc();
  lvconv = false;
  return true;
}
