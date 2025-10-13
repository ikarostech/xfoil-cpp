#include "XFoil.h"
#include <algorithm>
#include <cstring>
#include <cmath>
#include <numbers>
#include <sstream>
#include <iomanip>
#include <utility>
using namespace Eigen;

namespace {
constexpr double kAngleTolerance = 40.0;
}  // namespace


/** Loads the Foil's geometry in XFoil,
 *  calculates the normal vectors,
 *  and sets the results in current foil */
bool XFoil::initXFoilGeometry(int fn, const double *fx, const double *fy) {

  Matrix2Xd buffer_points = Matrix2Xd::Zero(2, fn);
  for (int i = 0; i < fn; i++) {
    buffer_points.col(i).x() = fx[i];
    buffer_points.col(i).y() = fy[i];
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
  const int point_count = foil.foil_shape.n;

  if (!ladij) {
    //----- calculate source influence matrix for airfoil surface if it doesn't
    // exist
    bij.block(0, 0, point_count + 1, point_count) =
        psi_gamma_lu.solve(bij.block(0, 0, point_count + 1, point_count)).eval();

    //------- store resulting dgam/dsig = dqtan/dsig vector
    dij.block(0, 0, point_count, point_count) = bij.block(0, 0, point_count, point_count);

    ladij = true;
  }

  //---- set up coefficient matrix of dpsi/dm on airfoil surface
  for (int i = 0; i < point_count; i++) {
    PsiResult psi_result =
        pswlin(foil, i, foil.foil_shape.points.col(i),
               foil.foil_shape.normal_vector.col(i), point_count, nw, apanel);
    bij.row(i).segment(point_count, nw) = -psi_result.dzdm.segment(point_count, nw).transpose();
  }

  //---- set up kutta condition (no direct source influence)

  bij.row(point_count).segment(point_count, nw).setZero();

  //---- multiply by inverse of factored dpsi/dgam matrix
  bij.block(0, point_count, point_count + 1, nw) =
      psi_gamma_lu.solve(bij.block(0, point_count, point_count + 1, nw)).eval();
  //---- set the source influence matrix for the wake sources
  dij.block(0, point_count, point_count, nw) = bij.block(0, point_count, point_count, nw);

  //**** now we need to calculate the influence of sources on the wake
  // velocities

  //---- calculate dqtan/dgam and dqtan/dsig at the wake points
  MatrixXd cij = MatrixXd::Zero(nw, point_count);
  for (int i = point_count; i < point_count + nw; i++) {
    int iw = i - point_count;
    //------ airfoil contribution at wake panel node
    PsiResult psi_result =
        psilin(foil, i, foil.wake_shape.points.col(i),
               foil.wake_shape.normal_vector.col(i), true, point_count, gamu,
               surface_vortex, alfa, qinf, apanel, foil.edge.sharp,
               foil.edge.ante, foil.edge.dste, foil.edge.aste);
    cij.row(iw) = psi_result.dqdg.head(point_count).transpose();
    dij.row(i).head(point_count) = psi_result.dqdm.head(point_count).transpose();
    //------ wake contribution
    psi_result =
        pswlin(foil, i, foil.wake_shape.points.col(i),
               foil.wake_shape.normal_vector.col(i), point_count, nw, apanel);
    dij.row(i).segment(point_count, nw) = psi_result.dqdm.segment(point_count, nw).transpose();
  }

  //---- add on effect of all sources on airfoil vorticity which effects wake
  // qtan
  dij.block(point_count, 0, nw, point_count) += cij * dij.topLeftCorner(point_count, point_count);

  dij.block(point_count, point_count, nw, nw) += cij * bij.block(0, point_count, point_count, nw);

  //---- make sure first wake point has same velocity as trailing edge
  dij.row(point_count) = dij.row(point_count - 1);

  lwdij = true;
  return true;
}


/** -------------------------------------------------------
 *     sets inviscid panel tangential velocity for
 *      current alpha.
 * -------------------------------------------------------- */
XFoil::TangentialVelocityResult XFoil::qiset() const {
  Matrix2d rotateMatrix =
      Matrix2d{{cos(alfa), sin(alfa)}, {-sin(alfa), cos(alfa)}};
  const int point_count = foil.foil_shape.n;
  const int total_nodes_with_wake = point_count + nw;

  TangentialVelocityResult result;
  result.qinv = VectorXd::Zero(total_nodes_with_wake);
  result.qinv_a = VectorXd::Zero(total_nodes_with_wake);

  for (int i = 0; i < total_nodes_with_wake; i++) {
    result.qinv[i] = rotateMatrix.row(0).dot(qinvu.col(i));
    result.qinv_a[i] = rotateMatrix.row(1).dot(qinvu.col(i));
  }

  return result;
}


/** -------------------------------------------------------------
 *     sets panel viscous tangential velocity from viscous ue
 * -------------------------------------------------------------- */
VectorXd XFoil::qvfue() const {
  VectorXd updated_qvis = qvis;
  for (int is = 1; is <= 2; is++) {
    const auto& vti_side = vti.get(is);
    const auto& uedg_side = uedg.get(is);
    const int limit = nbl.get(is) - 1;
    for (int ibl = 0; ibl < limit; ++ibl) {
      int i = ipan.get(is)[ibl];
      updated_qvis[i] = vti_side[ibl] * uedg_side[ibl];
    }
  }

  return updated_qvis;
}


/** ---------------------------------------------------------------
 *      sets inviscid tangential velocity for alpha = 0, 90
 *      on wake due to freestream and airfoil surface vorticity.
 * --------------------------------------------------------------- */
Matrix2Xd XFoil::qwcalc() {
  const int point_count = foil.foil_shape.n;
  Matrix2Xd updated_qinvu = qinvu;

  if (point_count >= 1 && point_count < updated_qinvu.cols()) {
    updated_qinvu.col(point_count) = updated_qinvu.col(point_count - 1);
  }

  for (int i = point_count + 1; i < point_count + nw; i++) {
    updated_qinvu.col(i) =
        psilin(foil, i, foil.wake_shape.points.col(i),
               foil.wake_shape.normal_vector.col(i), false, point_count, gamu,
               surface_vortex, alfa, qinf, apanel, foil.edge.sharp,
               foil.edge.ante, foil.edge.dste, foil.edge.aste)
            .qtan;
  }

  return updated_qinvu;
}

XFoil::EdgeVelocitySwapResult XFoil::swapEdgeVelocities(
    const SidePair<VectorXd> &usav) const {
  EdgeVelocitySwapResult result;
  result.swappedUsav = usav;
  result.restoredUedg = uedg;
  for (int is = 1; is <= 2; ++is) {
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      result.swappedUsav.get(is)[ibl] = uedg.get(is)[ibl];
      result.restoredUedg.get(is)[ibl] = usav.get(is)[ibl];
    }
  }
  return result;
}


XFoil::LeTeSensitivities XFoil::computeLeTeSensitivities(int ile1, int ile2,
                                                          int ite1,
                                                          int ite2) const {
  LeTeSensitivities sensitivities;
  sensitivities.ule1_m = VectorXd::Zero(2 * IVX + 1);
  sensitivities.ule2_m = VectorXd::Zero(2 * IVX + 1);
  sensitivities.ute1_m = VectorXd::Zero(2 * IVX + 1);
  sensitivities.ute2_m = VectorXd::Zero(2 * IVX + 1);
  for (int js = 1; js <= 2; ++js) {
    for (int jbl = 0; jbl < nbl.get(js) - 1; ++jbl) {
      const int j = ipan.get(js)[jbl];
      const int jv = isys.get(js)[jbl];
      const double vti_js = vti.get(js)[jbl];
      sensitivities.ule1_m[jv] = -vti.top[0] * vti_js * dij(ile1, j);
      sensitivities.ule2_m[jv] = -vti.bottom[0] * vti_js * dij(ile2, j);
      sensitivities.ute1_m[jv] = -vti.top[iblte.top] * vti_js * dij(ite1, j);
      sensitivities.ute2_m[jv] = -vti.bottom[iblte.bottom] * vti_js * dij(ite2, j);
    }
  }
  return sensitivities;
}


XFoil::DerivativeVectors XFoil::clearDerivativeVectors(const VectorXd &u_m,
                                                        const VectorXd &d_m) const {
  DerivativeVectors result{u_m, d_m};
  for (int js = 1; js <= 2; ++js) {
    for (int jbl = 0; jbl < nbl.get(js) - 1; ++jbl) {
      const int jv = isys.get(js)[jbl];
      result.u[jv] = 0.0;
      result.d[jv] = 0.0;
    }
  }
  return result;
}


/**
 * @brief Compute new edge velocities and their sensitivities.
 */
XFoil::EdgeVelocityDistribution XFoil::computeNewUeDistribution() const {
  EdgeVelocityDistribution distribution;
  distribution.unew.top = VectorXd::Zero(IVX);
  distribution.unew.bottom = VectorXd::Zero(IVX);
  distribution.u_ac.top = VectorXd::Zero(IVX);
  distribution.u_ac.bottom = VectorXd::Zero(IVX);
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      const int i = ipan.get(is)[ibl];
      double dui = 0.0;
      double dui_ac = 0.0;
      for (int js = 1; js <= 2; js++) {
        for (int jbl = 0; jbl < nbl.get(js) - 1; ++jbl) {
          const int j = ipan.get(js)[jbl];
          const int jv = isys.get(js)[jbl];
          const double ue_m = -vti.get(is)[ibl] * vti.get(js)[jbl] *
                              dij(i, j);
          dui += ue_m * (mass.get(js)[jbl] + vdel[jv](2, 0));
          dui_ac += ue_m * (-vdel[jv](2, 1));
        }
      }

      const double uinv_ac = lalfa ? 0.0 : uinv_a.get(is)[ibl];
      // Store unew/u_ac at 0-based station index
      distribution.unew.get(is)[ibl] = uinv.get(is)[ibl] + dui;
      distribution.u_ac.get(is)[ibl] = uinv_ac + dui_ac;
    }
  }
  return distribution;
}


/**
 * @brief Convert edge velocities to tangential velocities.
 */
XFoil::QtanResult XFoil::computeQtan(const SidePair<VectorXd> &unew,
                                      const SidePair<VectorXd> &u_ac) const {
  QtanResult result;
  result.qnew = VectorXd::Zero(IQX);
  result.q_ac = VectorXd::Zero(IQX);
  for (int is = 1; is <= 2; is++) {
    const VectorXd &unew_vec = (is == 1) ? unew.top : unew.bottom;
    const VectorXd &uac_vec = (is == 1) ? u_ac.top : u_ac.bottom;
    for (int ibl = 0; ibl < iblte.get(is); ++ibl) {
      const int i = ipan.get(is)[ibl];
      result.qnew[i] = vti.get(is)[ibl] * unew_vec[ibl];
      result.q_ac[i] = vti.get(is)[ibl] * uac_vec[ibl];
    }
  }
  return result;
}


/**
 * @brief Calculate lift coefficient contributions from tangential velocity.
 */
XFoil::ClContributions XFoil::computeClFromQtan(const VectorXd &qnew,
                                                const VectorXd &q_ac) const {
  ClContributions contributions;
  const auto compressibility = buildCompressibilityParams();
  const Matrix2d rotateMatrix = buildBodyToFreestreamRotation();
  const int point_count = foil.foil_shape.n;
  if (point_count == 0) {
    return contributions;
  }

  const PressureCoefficientResult cp_first = computePressureCoefficient(
      qnew[0], q_ac[0], compressibility);

  double cpg1 = cp_first.cp;
  double cpg1_ms = cp_first.cp_msq;
  double cpg1_ac = cp_first.cp_velocity_derivative;

  for (int i = 0; i < point_count; i++) {
    const int ip = (i + 1) % point_count;
    const PressureCoefficientResult cp_next = computePressureCoefficient(
        qnew[ip], q_ac[ip], compressibility);

    const double cpg2 = cp_next.cp;
    const double cpg2_ms = cp_next.cp_msq;
    const double cpg2_ac = cp_next.cp_velocity_derivative;

    const Vector2d dpoint =
        rotateMatrix * (foil.foil_shape.points.col(ip) - foil.foil_shape.points.col(i));

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


static void applyRelaxationLimit(const VectorXd &dn, double dhi, double dlo,
                                 double &relaxation) {
  double max_pos = 0.0;
  double min_neg = 0.0;
  for (const double value : dn) {
    if (value > max_pos)
      max_pos = value;
    if (value < min_neg)
      min_neg = value;
  }
  if (max_pos > 0.0)
    relaxation = std::min(relaxation, dhi / max_pos);
  if (min_neg < 0.0)
    relaxation = std::min(relaxation, dlo / min_neg);
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


double XFoil::clampRelaxationForGlobalChange(double relaxation, double dac,
                                             double lower,
                                             double upper) const {
  if (dac == 0.0)
    return relaxation;
  if (relaxation * dac > upper)
    relaxation = upper / dac;
  if (relaxation * dac < lower)
    relaxation = lower / dac;
  return relaxation;
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
    double &relaxation) const {
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


XFoil::BoundaryLayerSideState XFoil::applyBoundaryLayerDelta(
    int side, const BoundaryLayerDelta &delta, double relaxation) {
  BoundaryLayerSideState state;
  state.ctau = ctau.get(side);
  state.thet = thet.get(side);
  state.dstr = dstr.get(side);
  state.uedg = uedg.get(side);
  state.mass = mass.get(side);

  const int len = delta.dctau.size();
  if (len <= 0)
    return state;

  state.ctau.head(len) += relaxation * delta.dctau;
  state.thet.head(len) += relaxation * delta.dthet;
  state.dstr.head(len) += relaxation * delta.ddstr;
  state.uedg.head(len) += relaxation * delta.duedg;

  const int transition_index = std::max(0, itran.get(side));
  for (int idx = transition_index; idx < len; ++idx) {
    state.ctau[idx] = std::min(state.ctau[idx], 0.25);
  }

  for (int ibl = 0; ibl < len; ++ibl) {
    double dswaki = 0.0;
    if (ibl > iblte.get(side)) {
      const int wake_index = ibl - (iblte.get(side) + 1);
      dswaki = wgap[wake_index];
    }

    const double hklim = (ibl <= iblte.get(side)) ? 1.02 : 1.00005;
    const double uedg_val = state.uedg[ibl];
    const double uedg_sq = uedg_val * uedg_val;
    const double denom = 1.0 - 0.5 * uedg_sq * hstinv;
    const double msq = uedg_sq * hstinv / (gamm1 * denom);
    double dsw = state.dstr[ibl] - dswaki;
    dslim(dsw, state.thet[ibl], msq, hklim);
    state.dstr[ibl] = dsw + dswaki;
    state.mass[ibl] = state.dstr[ibl] * state.uedg[ibl];
  }

  return state;
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
  const auto ue_distribution = computeNewUeDistribution();
  unew = ue_distribution.unew;
  u_ac = ue_distribution.u_ac;
  const auto qtan = computeQtan(unew, u_ac);
  const auto cl_contributions = computeClFromQtan(qtan.qnew, qtan.q_ac);

  //--- initialize under-relaxation factor
  rlx = 1.0;
  const double cl_target = lalfa ? cl : clspec;
  double dac = computeAcChange(cl_contributions.cl, cl, cl_target,
                               cl_contributions.cl_ac, cl_contributions.cl_a,
                               cl_contributions.cl_ms);

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
  SidePair<BoundaryLayerSideState> updated_boundary_layer;
  for (int side = 1; side <= 2; ++side) {
    deltas.get(side) =
        buildBoundaryLayerDelta(side, unew.get(side), u_ac.get(side), dac);
    metrics.get(side) =
        evaluateSegmentRelaxation(side, deltas.get(side), dhi, dlo, rlx);
    rmsbl += metrics.get(side).rmsContribution;
    rmxbl = std::max(rmxbl, metrics.get(side).maxChange);
    updated_boundary_layer.get(side) =
        applyBoundaryLayerDelta(side, deltas.get(side), rlx);
  }

  rmsbl = sqrt(rmsbl / (4.0 * double(nbl.top + nbl.bottom)));

  if (lalfa)
    cl = cl + rlx * dac;
  else
    alfa = alfa + rlx * dac;

  for (int side = 1; side <= 2; ++side) {
    ctau.get(side) = updated_boundary_layer.get(side).ctau;
    thet.get(side) = updated_boundary_layer.get(side).thet;
    dstr.get(side) = updated_boundary_layer.get(side).dstr;
    uedg.get(side) = updated_boundary_layer.get(side).uedg;
    mass.get(side) = updated_boundary_layer.get(side).mass;
  }

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
  const int point_count = foil.foil_shape.n;
  const int total_nodes_with_wake = point_count + nw;

  //---- calculate wake trajectory from current inviscid solution if necessary
  if (!lwake)
    xyWake();

  //	---- set velocities on wake from airfoil vorticity for alpha=0, 90
  qinvu = qwcalc();

  //	---- set velocities on airfoil and wake for initial alpha
  {
    auto qiset_result = qiset();
    qinv = std::move(qiset_result.qinv);
    qinv_a = std::move(qiset_result.qinv_a);
  }

  if (!lipan) {
    if (lblini)
      surface_vortex = gamqv();

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
    qvis = qvfue();

    if (lvisc) {
      cpv = cpcalc(total_nodes_with_wake, qvis, qinf, minf);
      cpi = cpcalc(total_nodes_with_wake, qinv, qinf, minf);
    } else
      cpi = cpcalc(point_count, qinv, qinf, minf);

    surface_vortex = gamqv();
    const auto cl_result = clcalc(cmref);
    applyClComputation(cl_result);
    cd = cdcalc();
  }

  //	---- set up source influence matrix if it doesn't exist
  if (!lwdij || !ladij)
    qdcalc();

  return true;
}


XFoil::ViscalEndResult XFoil::ViscalEnd() {
  ViscalEndResult result;
  const int total_nodes_with_wake = foil.foil_shape.n + nw;
  result.inviscidCp = cpcalc(total_nodes_with_wake, qinv, qinf, minf);
  result.viscousCp = cpcalc(total_nodes_with_wake, qvis, qinf, minf);
  return result;
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
    auto qiset_result = qiset();
    qinv = std::move(qiset_result.qinv);
    qinv_a = std::move(qiset_result.qinv_a);
    uicalc();
  }

  qvis = qvfue();  //	------ calculate edge velocities qvis(.) from uedg(..)
  surface_vortex = gamqv();  //	------ set gam distribution from qvis
  stmove(); //	------ relocate stagnation point

  //	------ set updated cl,cd
  const auto cl_result = clcalc(cmref);
  applyClComputation(cl_result);
  cd = cdcalc();

  if (rmsbl < eps1) {
    lvconv = true;
    avisc = alfa;
    mvisc = minf;
    writeString("----------CONVERGED----------\n\n");
  }

  return true;
}

bool XFoil::isValidFoilAngles(Matrix2Xd points) {

  double max_angle = cang(points);
  return max_angle <= kAngleTolerance;
}


bool XFoil::isValidFoilPointSize(Matrix2Xd points) {
  return points.cols() >= 3;
}
