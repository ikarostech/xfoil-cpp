#include "XFoil.h"
#include <numbers>
#include <cstring>

// Initialization and global state related member functions split from XFoil.cpp

bool XFoil::initialize() {
  dtor = std::numbers::pi / 180.0;

  // allocate arrays and clear containers
  initializeDataStructures();

  // reset numerical and physical variables
  resetVariables();

  //---- drop tolerance for bl system solver
  vaccel = 0.01;

  //---- set minf, reinf, based on current cl-dependence
  minf_cl = getActualMach(1.0, mach_type);
  reinf_cl = getActualReynolds(1.0, reynolds_type);

  return true;
}

void XFoil::initializeDataStructures() {
  apanel = VectorXd::Zero(n + nw);
  memset(blsav, 0, sizeof(blsav));

  bij = MatrixXd::Zero(IQX, IZX);
  dij = MatrixXd::Zero(IZX, IZX);
  cpi = VectorXd::Zero(n + nw);
  cpv = VectorXd::Zero(n);
  ctau.top = VectorXd::Zero(IVX);
  ctau.bottom = VectorXd::Zero(IVX);
  ctq.top = VectorXd::Zero(IVX);
  ctq.bottom = VectorXd::Zero(IVX);
  dstr.top = VectorXd::Zero(IVX);
  dstr.bottom = VectorXd::Zero(IVX);

  ipan.top = VectorXi::Zero(IVX);
  ipan.bottom = VectorXi::Zero(IVX);
  isys.top = VectorXi::Zero(IVX);
  isys.bottom = VectorXi::Zero(IVX);
  itran.top = 0;
  itran.bottom = 0;
  mass.top = VectorXd::Zero(IVX);
  mass.bottom = VectorXd::Zero(IVX);
  normal_vectors = Matrix2Xd::Zero(2, n + nw);
  gamu = Matrix2Xd::Zero(2, n + 1);
  surface_vortex = Matrix2Xd::Zero(2, n);
  memset(qf0, 0, sizeof(qf0));
  memset(qf1, 0, sizeof(qf1));
  memset(qf2, 0, sizeof(qf2));
  memset(qf3, 0, sizeof(qf3));
  qinv = VectorXd::Zero(n + nw);
  qinv_a = VectorXd::Zero(n + nw);
  qinvu = Matrix2Xd::Zero(2, n + nw);
  qvis = VectorXd::Zero(n + nw);
  spline_length.resize(n + nw + 1);
  snew = VectorXd::Zero(4 * IBX);
  thet.top = VectorXd::Zero(IVX);
  thet.bottom = VectorXd::Zero(IVX);
  uedg.top = VectorXd::Zero(IVX);
  uedg.bottom = VectorXd::Zero(IVX);
  uinv.top = VectorXd::Zero(IVX);
  uinv.bottom = VectorXd::Zero(IVX);
  uinv_a.top = VectorXd::Zero(IVX);
  uinv_a.bottom = VectorXd::Zero(IVX);
  vti.top = VectorXd::Zero(IVX);
  vti.bottom = VectorXd::Zero(IVX);
  dpoints_ds.resize(2, n);

  xssi.top = VectorXd::Zero(IVX);
  xssi.bottom = VectorXd::Zero(IVX);

  memset(wgap, 0, sizeof(wgap));
  va.resize(IVX, Matrix<double, 3, 2>::Zero());
  vb.resize(IVX, Matrix<double, 3, 2>::Zero());
  vdel.resize(IVX, Matrix<double, 3, 2>::Zero());
  memset(vm, 0, sizeof(vm));
  blc.clear();
  memset(vz, 0, sizeof(vz));

  memset(qgamm, 0, sizeof(qgamm));
}

void XFoil::resetFlags() {
  lgamu = lvisc = lwake = lblini = lipan = false;
  lqaij = ladij = lwdij = lvconv = false;
  sharp = lalfa = false;
  trforc = simi = tran = turb = wake = trfree = false;
}

void XFoil::resetVariables() {
  gamma = 1.4;
  gamm1 = gamma - 1.0;
  qinf = 1.0;
  cl = cm = cd = 0.0;
  sigte = gamte = 0.0;
  awake = avisc = 0.0;
  resetFlags();
  cmref = Vector2d{0.25, 0.0};
  waklen = 1.0;
  i_stagnation = 0;
  qinfbl = tkbl = tkbl_ms = 0.0;
  rstbl = rstbl_ms = 0.0;
  hstinv = hstinv_ms = 0.0;
  reybl = reybl_ms = reybl_re = 0.0;
  gm1bl = 0.0;
  xiforc = 0.0;
  amcrit = 0.0;
  alfa = amax = rmxbl = rmsbl = rlx = ante = clspec = 0.0;
  minf = reinf = 0.0;
  minf_cl = reinf_cl = 0.0;
  sle = chord = 0.0;
  cl_alf = cl_msq = 0.0;
  tklam = tkl_msq = 0.0;
  sst = sst_go = sst_gp = 0.0;
  dste = aste = 0.0;
  cfm = cfm_ms = cfm_re = 0.0;
  cfm_u1 = cfm_t1 = cfm_d1 = 0.0;
  cfm_u2 = cfm_t2 = cfm_d2 = 0.0;
  xt = xt_a1 = xt_ms = xt_re = xt_xf = 0.0;
  xt_x1 = xt_t1 = xt_d1 = xt_u1 = 0.0;
  xt_x2 = xt_t2 = xt_d2 = xt_u2 = 0.0;
}

bool XFoil::comset() {
  //---- set karman-tsien parameter tklam
  double beta, beta_msq;
  beta = sqrt(1.0 - minf * minf);
  beta_msq = -0.5 / beta;

  tklam = MathUtil::pow(minf / (1.0 + beta), 2);
  tkl_msq = 1.0 / MathUtil::pow(1.0 + beta, 2) -
            2.0 * tklam / (1.0 + beta) * beta_msq;
  return true;
}

void XFoil::writeString(std::string str) { *m_pOutStream << str; }

double XFoil::getActualMach(double cls, MachType mach_type) {
  const double cla = std::max(cls, 0.000001);
  switch (mach_type) {
  case MachType::CONSTANT: {
    minf = minf1;
    return 0.0;
  }
  case MachType::FIXED_LIFT: {
    minf = minf1 / sqrt(cla);
    return -0.5 * minf / cla;
  }
  case MachType::FIXED_LIFT_AND_DYNAMIC_PRESSURE: {
    minf = minf1;
    return 0.0;
  }
  default:
    return 0;
  }
}

double XFoil::getActualReynolds(double cls, ReynoldsType reynolds_type) {
  const double cla = std::max(cls, 0.000001);
  switch (reynolds_type) {
  case ReynoldsType::CONSTANT: {
    reinf = reinf1;
    return 0.0;
  }
  case ReynoldsType::FIXED_LIFT: {
    reinf = reinf1 / sqrt(cla);
    return -0.5 * reinf / cla;
  }
  case ReynoldsType::FIXED_LIFT_AND_DYNAMIC_PRESSURE: {
    reinf = reinf1 / cla;
    return -reinf / cla;
  }
  default:
    return 0;
  }
}

bool XFoil::restoreblData(int icom) {
  if (icom == 1) {
    blData1 = blsav[icom];
  } else if (icom == 2) {
    blData2 = blsav[icom];
  }
  return true;
}

bool XFoil::saveblData(int icom) {
  if (icom == 1) {
    blsav[icom] = blData1;
  } else {
    blsav[icom] = blData2;
  }
  return true;
}

bool XFoil::setMach() {
  minf_cl = getActualMach(1.0, mach_type);
  reinf_cl = getActualReynolds(1.0, reynolds_type);
  comset();
  cpi = cpcalc(n, qinv, qinf, minf);
  if (lvisc) {
    cpv = cpcalc(n + nw, qvis, qinf, minf);
  }
  clcalc(cmref);
  cdcalc();
  lvconv = false;
  return true;
}

