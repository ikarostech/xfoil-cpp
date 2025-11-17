#include "XFoil.h"
#include <algorithm>
#include <array>
#include <cstring>
#include <numbers>
#include <ranges>
#include <unordered_map>

// Initialization and global state related member functions split from XFoil.cpp

namespace {
struct InitState {
  double amax = 0.0;
  std::array<double, IQX + 1> qf0{};
  std::array<double, IQX + 1> qf1{};
  std::array<double, IQX + 1> qf2{};
  std::array<double, IQX + 1> qf3{};
  std::array<blData, 3> blsav{};
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
  dtor = std::numbers::pi / 180.0;

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
  const int total_nodes_with_wake = point_count + nw;
  apanel = VectorXd::Zero(total_nodes_with_wake);
  auto& cache = ensureInitState(this);
  cache.blsav.fill(blData{});

  boundaryLayerWorkflow.lattice.clear();
  boundaryLayerWorkflow.lattice.resizeSide(static_cast<int>(SideType::TOP), IVX);
  boundaryLayerWorkflow.lattice.resizeSide(static_cast<int>(SideType::BOTTOM), IVX);
  boundaryLayerWorkflow.lattice.stationToPanel.top.setZero();
  boundaryLayerWorkflow.lattice.stationToPanel.bottom.setZero();
  boundaryLayerWorkflow.lattice.stationToSystem.top.setZero();
  boundaryLayerWorkflow.lattice.stationToSystem.bottom.setZero();
  boundaryLayerWorkflow.lattice.transitionIndex.top = 0;
  boundaryLayerWorkflow.lattice.transitionIndex.bottom = 0;
  boundaryLayerWorkflow.lattice.ctau.top.setZero();
  boundaryLayerWorkflow.lattice.ctau.bottom.setZero();
  boundaryLayerWorkflow.lattice.ctq.top.setZero();
  boundaryLayerWorkflow.lattice.ctq.bottom.setZero();
  boundaryLayerWorkflow.lattice.dstr.top.setZero();
  boundaryLayerWorkflow.lattice.dstr.bottom.setZero();
  boundaryLayerWorkflow.lattice.thet.top.setZero();
  boundaryLayerWorkflow.lattice.thet.bottom.setZero();
  boundaryLayerWorkflow.lattice.uedg.top.setZero();
  boundaryLayerWorkflow.lattice.uedg.bottom.setZero();
  boundaryLayerWorkflow.lattice.xssi.top.setZero();
  boundaryLayerWorkflow.lattice.xssi.bottom.setZero();
  boundaryLayerWorkflow.lattice.uinv.top.setZero();
  boundaryLayerWorkflow.lattice.uinv.bottom.setZero();
  boundaryLayerWorkflow.lattice.uinv_a.top.setZero();
  boundaryLayerWorkflow.lattice.uinv_a.bottom.setZero();
  boundaryLayerWorkflow.lattice.mass.top.setZero();
  boundaryLayerWorkflow.lattice.mass.bottom.setZero();
  boundaryLayerWorkflow.lattice.vti.top.setZero();
  boundaryLayerWorkflow.lattice.vti.bottom.setZero();

  bij = MatrixXd::Zero(IQX, IZX);
  dij = MatrixXd::Zero(IZX, IZX);
  cpi = VectorXd::Zero(total_nodes_with_wake);
  cpv = VectorXd::Zero(point_count);
  gamu = Matrix2Xd::Zero(2, point_count + 1);
  surface_vortex = Matrix2Xd::Zero(2, point_count);
  std::ranges::fill(cache.qf0, 0.0);
  std::ranges::fill(cache.qf1, 0.0);
  std::ranges::fill(cache.qf2, 0.0);
  std::ranges::fill(cache.qf3, 0.0);
  qinv = VectorXd::Zero(total_nodes_with_wake);
  qinv_a = VectorXd::Zero(total_nodes_with_wake);
  qinvu = Matrix2Xd::Zero(2, total_nodes_with_wake);
  qvis = VectorXd::Zero(total_nodes_with_wake);

  memset(wgap, 0, sizeof(wgap));
  va.resize(IVX, Matrix3x2d::Zero());
  vb.resize(IVX, Matrix3x2d::Zero());
  vdel.resize(IVX, Matrix3x2d::Zero());
  memset(vm, 0, sizeof(vm));
  boundaryLayerWorkflow.blc.clear();
  memset(vz, 0, sizeof(vz));

  memset(qgamm, 0, sizeof(qgamm));
}

void XFoil::resetFlags() {
  lgamu = lwake = lblini = lipan = false;
  lqaij = ladij = lwdij = lvconv = false;
  analysis_state_.viscous = false;
  analysis_state_.controlByAlpha = false;
  foil.edge.sharp = false;
  trforc = simi = tran = turb = wake = trfree = false;
}

void XFoil::resetVariables() {
  boundaryLayerWorkflow.state.station1 = blData{};
  boundaryLayerWorkflow.state.station2 = blData{};
  boundaryLayerWorkflow.lattice.transitionIndex.top = 0;
  boundaryLayerWorkflow.lattice.transitionIndex.bottom = 0;
  boundaryLayerWorkflow.lattice.stationCount.top = 0;
  boundaryLayerWorkflow.lattice.stationCount.bottom = 0;
  boundaryLayerWorkflow.lattice.trailingEdgeIndex.top = 0;
  boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom = 0;

  gamma = 1.4;
  gamm1 = gamma - 1.0;
  analysis_state_.qinf = 1.0;
  cl = cm = cd = 0.0;
  sigte = gamte = 0.0;
  avisc = 0.0;
  resetFlags();
  cmref = Vector2d{0.25, 0.0};
  i_stagnation = 0;
  qinfbl = tkbl = tkbl_ms = 0.0;
  rstbl = rstbl_ms = 0.0;
  hstinv = hstinv_ms = 0.0;
  reybl = reybl_ms = reybl_re = 0.0;
  gm1bl = 0.0;
  xiforc = 0.0;
  amcrit = 0.0;
  auto& cache = ensureInitState(this);
  cache.amax = 0.0;
  analysis_state_.alpha = 0.0;
  analysis_state_.clspec = 0.0;
  rmxbl = rmsbl = rlx = 0.0;
  foil.edge.ante = 0.0;
  foil.edge.aste = 0.0;
  foil.edge.dste = 0.0;
  analysis_state_.currentMach = 0.0;
  analysis_state_.currentRe = 0.0;
  minf_cl = reinf_cl = 0.0;
  cl_alf = cl_msq = 0.0;
  tklam = tkl_msq = 0.0;
  sst = sst_go = sst_gp = 0.0;
  xt = xt_a1 = xt_ms = xt_re = xt_xf = 0.0;
  xt_x1 = xt_t1 = xt_d1 = xt_u1 = 0.0;
  xt_x2 = xt_t2 = xt_d2 = xt_u2 = 0.0;
}

bool XFoil::comset() {
  //---- set karman-tsien parameter tklam
  double beta, beta_msq;
  const double current_mach = analysis_state_.currentMach;
  beta = sqrt(1.0 - current_mach * current_mach);
  beta_msq = -0.5 / beta;

  tklam = MathUtil::pow(current_mach / (1.0 + beta), 2);
  tkl_msq = 1.0 / MathUtil::pow(1.0 + beta, 2) -
            2.0 * tklam / (1.0 + beta) * beta_msq;
  return true;
}

void XFoil::writeString(std::string str) { *m_pOutStream << str; }

double XFoil::getActualMach(double cls, MachType mach_control) {
  AnalysisState& state = analysis_state_;
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
  AnalysisState& state = analysis_state_;
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

bool XFoil::restoreblData(int icom) {
  auto& cache = ensureInitState(this);
  if (icom == 1) {
    boundaryLayerWorkflow.state.station1 = cache.blsav[icom];
  } else if (icom == 2) {
    boundaryLayerWorkflow.state.station2 = cache.blsav[icom];
  }
  return true;
}

bool XFoil::saveblData(int icom) {
  auto& cache = ensureInitState(this);
  if (icom == 1) {
    cache.blsav[icom] = boundaryLayerWorkflow.state.station1;
  } else {
    cache.blsav[icom] = boundaryLayerWorkflow.state.station2;
  }
  return true;
}

bool XFoil::setMach() {
  minf_cl = getActualMach(1.0, analysis_state_.machType);
  reinf_cl = getActualReynolds(1.0, analysis_state_.reynoldsType);
  comset();
  const int point_count = foil.foil_shape.n;
  cpi = cpcalc(point_count, qinv, analysis_state_.qinf, analysis_state_.currentMach);
  if (analysis_state_.viscous) {
    cpv = cpcalc(point_count + nw, qvis, analysis_state_.qinf, analysis_state_.currentMach);
  }
  const auto cl_result = clcalc(cmref);
  applyClComputation(cl_result);
  cd = cdcalc();
  lvconv = false;
  return true;
}
