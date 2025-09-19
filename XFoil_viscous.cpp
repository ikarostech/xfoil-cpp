#include "XFoil.h"
#include <algorithm>
#include <cstring>
#include <cmath>
#include <numbers>
#include <sstream>
#include <iomanip>
using namespace Eigen;

namespace {
constexpr double kAngleTolerance = 40.0;
}  // namespace


/** Loads the Foil's geometry in XFoil,
 *  calculates the normal vectors,
 *  and sets the results in current foil */
bool XFoil::initXFoilGeometry(int fn, const double *fx, const double *fy) {

  Matrix2Xd buffer_points = Matrix2Xd::Zero(2, fn + 1);
  for (int i = 0; i < fn; i++) {
    buffer_points.col(i + 1).x() = fx[i];
    buffer_points.col(i + 1).y() = fy[i];
  }

  if (!isValidFoilPointSize(buffer_points) ||
      !isValidFoilAngles(buffer_points)) {
    writeString("Unrecognized foil format");
    return false;
  }

  abcopy(buffer_points);
  return true;
}


bool XFoil::initXFoilAnalysis(double Re, double alpha, double Mach,
                              double NCrit, double XtrTop, double XtrBot,
                              ReynoldsType reType, MachType maType,
                              bool bViscous, std::stringstream &outStream) {
  // Sets Analysis parameters in XFoil
  m_pOutStream = &outStream;

  lblini = false;
  lipan = false;

  reinf1 = Re;
  alfa = alpha * std::numbers::pi / 180.0;

  minf1 = Mach;
  reynolds_type = reType;
  mach_type = maType;
  lalfa = true;
  qinf = 1.0;

  lvisc = bViscous;

  acrit = NCrit;
  xstrip.top = XtrTop;
  xstrip.bottom = XtrBot;

  if (Mach > 0.000001) {
    if (!setMach()) {
      writeString(
          "... Invalid Analysis Settings\nCpCalc: local speed too large\n "
          "Compressibility corrections invalid ");
      return false;
    }
  }

  return true;
}


/** -------------------------------------------
 *      sets actual mach, reynolds numbers
 *      from unit-cl values and specified cls
 *      depending on matyp,retyp flags.
 * -------------------------------------------- */
// moved to XFoil_init.cpp: getActualMach()/getActualReynolds()

// moved to XFoil_geometry.cpp: ncalc()

/** --------------------------------------------------------------------
 *	   Calculates current streamfunction psi and tangential velocity
 *	   qtan at panel node or wake node i due to freestream and wake
 *	   sources sig.  also calculates sensitivity vectors dpsi/dsig
 *	   (dzdm) and dqtan/dsig (dqdm).
 *
 *			airfoil:  1   < i < n
 *			wake:	  n+1 < i < n+nw
 *-------------------------------------------------------------------- */
/** ------------------------------------------------------
 *         calculates source panel influence coefficient
 * 	   matrix for current airfoil and wake geometry.
 * ------------------------------------------------------ */
bool XFoil::qdcalc() {
  // TRACE("calculating source influence matrix ...\n");
  writeString("   Calculating source influence matrix ...\n");

  if (!ladij) {
    //----- calculate source influence matrix for airfoil surface if it doesn't
    // exist
    bij.block(0, 0, n + 1, n) =
        psi_gamma_lu.solve(bij.block(0, 0, n + 1, n)).eval();

    //------- store resulting dgam/dsig = dqtan/dsig vector
    dij.block(0, 0, n, n) = bij.block(0, 0, n, n);

    ladij = true;
  }

  //---- set up coefficient matrix of dpsi/dm on airfoil surface
  for (int i = 0; i < n; i++) {
    PsiResult psi_result = pswlin(points, i, points.col(i),
                                  normal_vectors.col(i));
    bij.row(i).segment(n, nw) = -psi_result.dzdm.segment(n, nw).transpose();
  }

  //---- set up kutta condition (no direct source influence)

  bij.row(n).segment(n, nw).setZero();

  //---- multiply by inverse of factored dpsi/dgam matrix
  bij.block(0, n, n + 1, nw) =
      psi_gamma_lu.solve(bij.block(0, n, n + 1, nw)).eval();
  //---- set the source influence matrix for the wake sources
  dij.block(0, n, n, nw) = bij.block(0, n, n, nw);

  //**** now we need to calculate the influence of sources on the wake
  // velocities

  //---- calculate dqtan/dgam and dqtan/dsig at the wake points
  MatrixXd cij = MatrixXd::Zero(nw, n);
  for (int i = n; i < n + nw; i++) {
    int iw = i - n;
    //------ airfoil contribution at wake panel node
    PsiResult psi_result =
        psilin(points, i, points.col(i),
               normal_vectors.col(i), true);
    cij.row(iw) = psi_result.dqdg.head(n).transpose();
    dij.row(i).head(n) = psi_result.dqdm.head(n).transpose();
    //------ wake contribution
    psi_result = pswlin(points, i, points.col(i), normal_vectors.col(i));
    dij.row(i).segment(n, nw) = psi_result.dqdm.segment(n, nw).transpose();
  }

  //---- add on effect of all sources on airfoil vorticity which effects wake
  // qtan
  dij.block(n, 0, nw, n) += cij * dij.topLeftCorner(n, n);

  dij.block(n, n, nw, nw) += cij * bij.block(0, n, n, nw);

  //---- make sure first wake point has same velocity as trailing edge
  dij.row(n) = dij.row(n - 1);

  lwdij = true;
  return true;
}


/** -------------------------------------------------------
 *     sets inviscid panel tangential velocity for
 *      current alpha.
 * -------------------------------------------------------- */
bool XFoil::qiset() {
  Matrix2d rotateMatrix =
      Matrix2d{{cos(alfa), sin(alfa)}, {-sin(alfa), cos(alfa)}};

  for (int i = 0; i < n + nw; i++) {
    qinv[i] = rotateMatrix.row(0).dot(qinvu.col(i));
    qinv_a[i] = rotateMatrix.row(1).dot(qinvu.col(i));
  }

  return true;
}


/** -------------------------------------------------------------
 *     sets panel viscous tangential velocity from viscous ue
 * -------------------------------------------------------------- */
bool XFoil::qvfue() {
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      int i = ipan.get(is)[ibl];
      qvis[i] = vti.get(is)[ibl] * uedg.get(is)[ibl];
    }
  }

  return true;
}


/** ---------------------------------------------------------------
 *      sets inviscid tangential velocity for alpha = 0, 90
 *      on wake due to freestream and airfoil surface vorticity.
 * --------------------------------------------------------------- */
bool XFoil::qwcalc() {

  //---- first wake point (same as te)
  qinvu.col(n) = qinvu.col(n - 1);

  //---- rest of wake
  for (int i = n + 1; i < n + nw; i++) {
    qinvu.col(i) =
        psilin(points, i, points.col(i),
               normal_vectors.col(i), false)
            .qtan;
  }

  return true;
}


// moved to XFoil_init.cpp: restoreblData()/saveblData()

void XFoil::swapEdgeVelocities(SidePair<VectorXd> &usav) {
  for (int is = 1; is <= 2; ++is) {
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      double temp = usav.get(is)[ibl];
      double cur = uedg.get(is)[ibl];
      usav.get(is)[ibl] = cur;
      uedg.get(is)[ibl] = temp;
    }
  }
}


void XFoil::computeLeTeSensitivities(int ile1, int ile2, int ite1, int ite2,
                                     VectorXd &ule1_m, VectorXd &ule2_m,
                                     VectorXd &ute1_m, VectorXd &ute2_m) {
  for (int js = 1; js <= 2; ++js) {
    for (int jbl = 0; jbl < nbl.get(js) - 1; ++jbl) {
      int j = ipan.get(js)[jbl];
      int jv = isys.get(js)[jbl];
      ule1_m[jv] = -vti.top[0] * vti.get(js)[jbl] *
                   dij(ile1, j);
      ule2_m[jv] = -vti.bottom[0] * vti.get(js)[jbl] *
                   dij(ile2, j);
      ute1_m[jv] = -vti.top[iblte.top] * vti.get(js)[jbl] *
                   dij(ite1, j);
      ute2_m[jv] = -vti.bottom[iblte.bottom] * vti.get(js)[jbl] *
                   dij(ite2, j);
    }
  }
}


void XFoil::clearDerivativeVectors(VectorXd &u_m, VectorXd &d_m) {
  for (int js = 1; js <= 2; ++js) {
    for (int jbl = 0; jbl < nbl.get(js) - 1; ++jbl) {
      int jv = isys.get(js)[jbl];
      u_m[jv] = 0.0;
      d_m[jv] = 0.0;
    }
  }
}


/**
 * @brief Compute new edge velocities and their sensitivities.
 */
void XFoil::computeNewUeDistribution(SidePair<VectorXd> &unew,
                                     SidePair<VectorXd> &u_ac) {
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      int i = ipan.get(is)[ibl];
      double dui = 0.0;
      double dui_ac = 0.0;
      for (int js = 1; js <= 2; js++) {
        for (int jbl = 0; jbl < nbl.get(js) - 1; ++jbl) {
          int j = ipan.get(js)[jbl];
          int jv = isys.get(js)[jbl];
          double ue_m = -vti.get(is)[ibl] * vti.get(js)[jbl] *
                        dij(i, j);
          dui += ue_m * (mass.get(js)[jbl] + vdel[jv](2, 0));
          dui_ac += ue_m * (-vdel[jv](2, 1));
        }
      }

      double uinv_ac = lalfa ? 0.0 : uinv_a.get(is)[ibl];
      // Store unew/u_ac at 0-based station index
      unew.get(is)[ibl] = uinv.get(is)[ibl] + dui;
      u_ac.get(is)[ibl] = uinv_ac + dui_ac;
    }
  }
}


/**
 * @brief Convert edge velocities to tangential velocities.
 */
void XFoil::computeQtan(const SidePair<VectorXd> &unew,
                        const SidePair<VectorXd> &u_ac, VectorXd &qnew,
                        VectorXd &q_ac) {
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 0; ibl < iblte.get(is); ++ibl) {
      int i = ipan.get(is)[ibl];
      const VectorXd &unew_vec = (is == 1) ? unew.top : unew.bottom;
      const VectorXd &uac_vec = (is == 1) ? u_ac.top : u_ac.bottom;
      qnew[i] = vti.get(is)[ibl] * unew_vec[ibl];
      q_ac[i] = vti.get(is)[ibl] * uac_vec[ibl];
    }
  }
}


/**
 * @brief Calculate lift coefficient contributions from tangential velocity.
 */
void XFoil::computeClFromQtan(const VectorXd &qnew, const VectorXd &q_ac,
                              double &clnew, double &cl_a, double &cl_ms,
                              double &cl_ac) {
  double beta = sqrt(1.0 - minf * minf);
  double beta_msq = -0.5 / beta;

  double bfac = 0.5 * minf * minf / (1.0 + beta);
  double bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;

  clnew = 0.0;
  cl_a = 0.0;
  cl_ms = 0.0;
  cl_ac = 0.0;

  double cginc = 1.0 - (qnew[0] / qinf) * (qnew[0] / qinf);
  double cpg1 = cginc / (beta + bfac * cginc);
  double cpg1_ms =
      -cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

  double cpi_q = -2.0 * qnew[0] / qinf / qinf;
  double cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
  double cpg1_ac = cpc_cpi * cpi_q * q_ac[0];

  for (int i = 0; i < n; i++) {
    int ip = (i + 1) % n;
    cginc = 1.0 - (qnew[ip] / qinf) * (qnew[ip] / qinf);
    double cpg2 = cginc / (beta + bfac * cginc);
    double cpg2_ms =
        -cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

    cpi_q = -2.0 * qnew[ip] / qinf / qinf;
    cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
    double cpg2_ac = cpc_cpi * cpi_q * q_ac[ip];

    Matrix2d rotateMatrix =
        Matrix2d{{cos(alfa), sin(alfa)}, {-sin(alfa), cos(alfa)}};
    Vector2d dpoint = rotateMatrix * (points.col(ip) - points.col(i));

    const double ag = 0.5 * (cpg2 + cpg1);
    const double ag_ms = 0.5 * (cpg2_ms + cpg1_ms);
    const double ag_ac = 0.5 * (cpg2_ac + cpg1_ac);

    clnew += dpoint.x() * ag;
    cl_a += dpoint.y() * ag;
    cl_ms += dpoint.x() * ag_ms;
    cl_ac += dpoint.x() * ag_ac;

    cpg1 = cpg2;
    cpg1_ms = cpg2_ms;
    cpg1_ac = cpg2_ac;
  }
}


static void applyRelaxationLimit(const VectorXd &dn, double dhi, double dlo,
                                 double &rlx) {
  double max_pos = 0.0;
  double min_neg = 0.0;
  for (int i = 0; i < dn.size(); ++i) {
    const double value = dn[i];
    if (value > max_pos)
      max_pos = value;
    if (value < min_neg)
      min_neg = value;
  }
  if (max_pos > 0.0)
    rlx = std::min(rlx, dhi / max_pos);
  if (min_neg < 0.0)
    rlx = std::min(rlx, dlo / min_neg);
}


double XFoil::computeAcChange(double clnew, double cl_current,
                              double cl_target, double cl_ac, double cl_a,
                              double cl_ms) const {
  if (lalfa) {
    return (clnew - cl_current) /
           (1.0 - cl_ac - cl_ms * 2.0 * minf * minf_cl);
  }
  return (clnew - cl_target) / (0.0 - cl_ac - cl_a);
}


double XFoil::clampRelaxationForGlobalChange(double rlx, double dac,
                                             double lower,
                                             double upper) const {
  if (dac == 0.0)
    return rlx;
  if (rlx * dac > upper)
    rlx = upper / dac;
  if (rlx * dac < lower)
    rlx = lower / dac;
  return rlx;
}


XFoil::BoundaryLayerDelta XFoil::buildBoundaryLayerDelta(
    int side, const VectorXd &unew_side, const VectorXd &u_ac_side,
    double dac) const {
  BoundaryLayerDelta delta;
  const int len = nbl.get(side) - 1;
  if (len <= 0)
    return delta;

  delta.dctau = VectorXd(len);
  delta.dthet = VectorXd(len);
  delta.ddstr = VectorXd(len);
  delta.duedg = VectorXd(len);

  const auto iv = isys.get(side).segment(0, len);
  VectorXd dmass(len);
  for (int j = 0; j < len; ++j) {
    const int idx = iv[j];
    delta.dctau[j] = vdel[idx](0, 0) - dac * vdel[idx](0, 1);
    delta.dthet[j] = vdel[idx](1, 0) - dac * vdel[idx](1, 1);
    dmass[j] = vdel[idx](2, 0) - dac * vdel[idx](2, 1);
  }

  const VectorXd uedg_segment = uedg.get(side).head(len);
  const VectorXd dstr_segment = dstr.get(side).head(len);
  const VectorXd unew_segment = unew_side.head(len);
  const VectorXd uac_segment = u_ac_side.head(len);

  delta.duedg = unew_segment + dac * uac_segment - uedg_segment;
  delta.ddstr = (dmass - dstr_segment.cwiseProduct(delta.duedg))
                    .cwiseQuotient(uedg_segment);

  return delta;
}


XFoil::BoundaryLayerMetrics XFoil::evaluateSegmentRelaxation(
    int side, const BoundaryLayerDelta &delta, double dhi, double dlo,
    double &rlx) const {
  BoundaryLayerMetrics metrics;
  const int len = delta.dctau.size();
  if (len <= 0)
    return metrics;

  const VectorXd ctau_segment = ctau.get(side).head(len);
  const VectorXd thet_segment = thet.get(side).head(len);
  const VectorXd dstr_segment = dstr.get(side).head(len);

  VectorXd dn1(len);
  const int transition_index = itran.get(side);
  for (int idx = 0; idx < len; ++idx) {
    dn1[idx] =
        (idx < transition_index) ? delta.dctau[idx] / 10.0
                                 : delta.dctau[idx] / ctau_segment[idx];
  }
  const VectorXd dn2 = delta.dthet.cwiseQuotient(thet_segment);
  const VectorXd dn3 = delta.ddstr.cwiseQuotient(dstr_segment);
  const VectorXd dn4 = delta.duedg.array().abs() / 0.25;

  applyRelaxationLimit(dn1, dhi, dlo, rlx);
  applyRelaxationLimit(dn2, dhi, dlo, rlx);
  applyRelaxationLimit(dn3, dhi, dlo, rlx);
  applyRelaxationLimit(dn4, dhi, dlo, rlx);

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


void XFoil::applyBoundaryLayerDelta(int side,
                                    const BoundaryLayerDelta &delta,
                                    double rlx) {
  const int len = delta.dctau.size();
  if (len <= 0)
    return;

  auto ctau_segment = ctau.get(side).head(len);
  auto thet_segment = thet.get(side).head(len);
  auto dstr_segment = dstr.get(side).head(len);
  auto uedg_segment = uedg.get(side).head(len);

  ctau_segment += rlx * delta.dctau;
  thet_segment += rlx * delta.dthet;
  dstr_segment += rlx * delta.ddstr;
  uedg_segment += rlx * delta.duedg;

  const int transition_index = std::max(0, itran.get(side));
  for (int idx = transition_index; idx < len; ++idx) {
    ctau_segment[idx] = std::min(ctau_segment[idx], 0.25);
  }

  for (int ibl = 0; ibl < len; ++ibl) {
    double dswaki = 0.0;
    if (ibl > iblte.get(side)) {
      const int wake_index = ibl - (iblte.get(side) + 1);
      dswaki = wgap[wake_index];
    }

    const double hklim = (ibl <= iblte.get(side)) ? 1.02 : 1.00005;
    const double uedg_val = uedg.get(side)[ibl];
    const double uedg_sq = uedg_val * uedg_val;
    const double denom = 1.0 - 0.5 * uedg_sq * hstinv;
    const double msq = uedg_sq * hstinv / (gamm1 * denom);
    double dsw = dstr.get(side)[ibl] - dswaki;
    dslim(dsw, thet.get(side)[ibl], msq, hklim);
    dstr.get(side)[ibl] = dsw + dswaki;
    mass.get(side)[ibl] = dstr.get(side)[ibl] * uedg.get(side)[ibl];
  }
}


bool XFoil::update() {
  //------------------------------------------------------------------
  //      adds on newton deltas to boundary layer variables.
  //      checks for excessive changes and underrelaxes if necessary.
  //      calculates max and rms changes.
  //      also calculates the change in the global variable "ac".
  //        if lalfa=true , "ac" is cl
  //        if lalfa=false, "ac" is alpha
  //------------------------------------------------------------------

  SidePair<VectorXd> unew, u_ac;
  unew.top = VectorXd::Zero(IVX);
  unew.bottom = VectorXd::Zero(IVX);
  u_ac.top = VectorXd::Zero(IVX);
  u_ac.bottom = VectorXd::Zero(IVX);

  VectorXd qnew = VectorXd::Zero(IQX);
  VectorXd q_ac = VectorXd::Zero(IQX);
  double clnew = 0.0, cl_a = 0.0, cl_ms = 0.0, cl_ac = 0.0;

  //---- max allowable alpha changes per iteration
  const double dalmax = 0.5 * dtor;
  const double dalmin = -0.5 * dtor;
  //---- max allowable cl change per iteration
  double dclmax = 0.5;
  double dclmin = -0.5;
  if (mach_type != MachType::CONSTANT)
    dclmin = std::max(-0.5, -0.9 * cl);
  hstinv =
      gamm1 * (minf / qinf) * (minf / qinf) / (1.0 + 0.5 * gamm1 * minf * minf);

  //--- calculate new ue distribution and tangential velocities
  computeNewUeDistribution(unew, u_ac);
  computeQtan(unew, u_ac, qnew, q_ac);
  computeClFromQtan(qnew, q_ac, clnew, cl_a, cl_ms, cl_ac);

  //--- initialize under-relaxation factor
  rlx = 1.0;
  const double cl_target = lalfa ? cl : clspec;
  double dac = computeAcChange(clnew, cl, cl_target, cl_ac, cl_a, cl_ms);

  if (lalfa)
    rlx = clampRelaxationForGlobalChange(rlx, dac, dclmin, dclmax);
  else
    rlx = clampRelaxationForGlobalChange(rlx, dac, dalmin, dalmax);

  rmsbl = 0.0;
  rmxbl = 0.0;
  const double dhi = 1.5;
  const double dlo = -0.5;

  SidePair<BoundaryLayerDelta> deltas;
  SidePair<BoundaryLayerMetrics> metrics;
  for (int side = 1; side <= 2; ++side) {
    deltas.get(side) =
        buildBoundaryLayerDelta(side, unew.get(side), u_ac.get(side), dac);
    metrics.get(side) =
        evaluateSegmentRelaxation(side, deltas.get(side), dhi, dlo, rlx);
    rmsbl += metrics.get(side).rmsContribution;
    rmxbl = std::max(rmxbl, metrics.get(side).maxChange);
  }

  rmsbl = sqrt(rmsbl / (4.0 * double(nbl.top + nbl.bottom)));

  if (lalfa)
    cl = cl + rlx * dac;
  else
    alfa = alfa + rlx * dac;

  for (int side = 1; side <= 2; ++side)
    applyBoundaryLayerDelta(side, deltas.get(side), rlx);

  //--- equate upper wake arrays to lower wake arrays
  for (int kbl = 1; kbl <= nbl.bottom - (iblte.bottom + 1); kbl++) {
    ctau.top[iblte.top + kbl] = ctau.bottom[iblte.bottom + kbl];
    thet.top[iblte.top + kbl] = thet.bottom[iblte.bottom + kbl];
    dstr.top[iblte.top + kbl] = dstr.bottom[iblte.bottom + kbl];
    uedg.top[iblte.top + kbl] =
                     uedg.bottom[iblte.bottom + kbl];
    ctq.top[iblte.top + kbl] = ctq.bottom[iblte.bottom + kbl];
  }

  return true;
}


bool XFoil::viscal() {
  ////--------------------------------------
  //     converges viscous operating point
  ////--------------------------------------

  //---- calculate wake trajectory from current inviscid solution if necessary
  if (!lwake)
    xyWake();

  //	---- set velocities on wake from airfoil vorticity for alpha=0, 90
  qwcalc();

  //	---- set velocities on airfoil and wake for initial alpha
  qiset();

  if (!lipan) {
    if (lblini)
      gamqv();

    //	----- locate stagnation point arc length position and panel index
    stfind();

    //	----- set  bl position -> panel position  pointers
    iblpan();

    //	----- calculate surface arc length array for current stagnation point
    // location
    xicalc();

    //	----- set  bl position -> system line  pointers
    iblsys();
  }

  //	---- set inviscid bl edge velocity uinv from qinv
  uicalc();

  if (!lblini) {
    //	----- set initial ue from inviscid ue
    for (int ibl = 0; ibl < nbl.top - 1; ibl++) {
      uedg.top[ibl] = uinv.top[ibl];
    }
    for (int ibl = 0; ibl < nbl.bottom - 1; ibl++) {
      uedg.bottom[ibl] = uinv.bottom[ibl];
    }
  }

  if (lvconv) {
    //	----- set correct cl if converged point exists
    qvfue();

    if (lvisc) {
      cpv = cpcalc(n + nw, qvis, qinf, minf);
      cpi = cpcalc(n + nw, qinv, qinf, minf);
    } else
      cpi = cpcalc(n, qinv, qinf, minf);

    gamqv();
    clcalc(cmref);
    cdcalc();
  }

  //	---- set up source influence matrix if it doesn't exist
  if (!lwdij || !ladij)
    qdcalc();

  return true;
}


bool XFoil::ViscalEnd() {

  cpi = cpcalc(n + nw, qinv, qinf, minf);
  cpv = cpcalc(n + nw, qvis, qinf, minf);

  return true;
}


bool XFoil::ViscousIter() {
  //	Performs one iteration
  std::stringstream ss;
  double eps1 = 0.0001;

  setbl(); //	------ fill newton system for bl variables

  blsolve(); //	------ solve newton system with custom solver

  update(); //	------ update bl variables

  if (lalfa) { //	------- set new freestream mach, re from new cl
    minf_cl = getActualMach(cl, mach_type);
    reinf_cl = getActualReynolds(cl, reynolds_type);
    comset();
  } else { //	------- set new inviscid speeds qinv and uinv for new alpha
    qiset();
    uicalc();
  }

  qvfue();  //	------ calculate edge velocities qvis(.) from uedg(..)
  gamqv();  //	------ set gam distribution from qvis
  stmove(); //	------ relocate stagnation point

  //	------ set updated cl,cd
  clcalc(cmref);
  cdcalc();

  if (rmsbl < eps1) {
    lvconv = true;
    avisc = alfa;
    mvisc = minf;
    writeString("----------CONVERGED----------\n\n");
  }

  return true;
}


bool XFoil::setexp(double spline_length[], double ds1, double smax, int nn) {
  //........................................................
  //     sets geometriy stretched array s:
  //
  //       s(i+1) - s(i)  =  r * [s(i) - s(i-1)]
  //
  //       s     (output)  array to be set
  //       ds1   (input)   first s increment:  spline_length[2] -
  //       spline_length[1] smax  (input)   final s value:      s(nn) nn (input)
  //       number of points
  //........................................................
  int nex;
  double sigma, rnex, rni, aaa, bbb, ccc;
  double disc, ratio, sigman, res;
  double dresdr, dratio, ds;

  sigma = smax / ds1;
  nex = nn - 1;
  rnex = (double)nex;
  rni = 1.0 / rnex;

  //-- solve quadratic for initial geometric ratio guess
  aaa = rnex * (rnex - 1.0) * (rnex - 2.0) / 6.0;
  bbb = rnex * (rnex - 1.0) / 2.0;
  ccc = rnex - sigma;

  disc = bbb * bbb - 4.0 * aaa * ccc;
  disc = std::max(0.0, disc);

  if (nex <= 1) {
    writeString("setexp: cannot fill array.  n too small\n");
    return false;
  } else {
    if (nex == 2)
      ratio = -ccc / bbb + 1.0;
    else
      ratio = (-bbb + sqrt(disc)) / (2.0 * aaa) + 1.0;
  }

  //-- newton iteration for actual geometric ratio
  for (int iter = 1; iter <= 100; iter++) {
    sigman = (pow(ratio, (double)nex) - 1.0) / (ratio - 1.0);
    res = pow(sigman, rni) - pow(sigma, rni);
    dresdr = rni * pow(sigman, rni) *
             (rnex * pow(ratio, (double)(nex - 1)) - sigman) /
             (pow(ratio, (double)nex) - 1.0);

    dratio = -res / dresdr;
    ratio = ratio + dratio;

    if (fabs(dratio) < 1.0e-5) {
      break;
    }
  }

  // Fill 0-based: length[0]=0, then cumulative with geometric step
  spline_length[0] = 0.0;
  ds = ds1;
  for (int i = 1; i < nn; i++) {
    spline_length[i] = spline_length[i - 1] + ds;
    ds = ds * ratio;
  }
  return true;
}


// moved to XFoil_geometry.cpp: xyWake()

bool XFoil::isValidFoilAngles(Matrix2Xd points) {

  double max_angle = cang(points);
  return max_angle <= kAngleTolerance;
}


bool XFoil::isValidFoilPointSize(Matrix2Xd points) {
  return points.cols() >= 3;
}
