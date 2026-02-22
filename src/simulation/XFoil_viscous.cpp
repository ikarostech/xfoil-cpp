#include "XFoil.h"
#include "simulation/InviscidSolver.hpp"
#include "simulation/Blsolve.hpp"
#include "simulation/BoundaryLayer_march.hpp"
#include "simulation/psi.hpp"
#include "domain/boundary_layer/boundary_layer_builder.hpp"
#include "core/math_util.hpp"
#include "infrastructure/logger.hpp"
#include <algorithm>
#include <cstring>
#include <cmath>
#include <numbers>
#include <sstream>
#include <iomanip>
#include <utility>
using Eigen::Matrix2Xd;
using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::VectorXd;

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
    Logger::instance().write("Unrecognized foil format");
    return false;
  }

  abcopy(buffer_points);
  aerodynamicCache = ggcalc();
  return true;
}


bool XFoil::initXFoilAnalysis(double Re, double alpha, double Mach,
                              double NCrit, double XtrTop, double XtrBot,
                              ReynoldsType reType, MachType maType,
                              bool bViscous) {

  setBLInitialized(false);
  invalidatePanelMap();
  invalidateWakeGeometry();
  invalidateConvergedSolution();

  FlowState& state = analysis_state_;
  state.referenceRe = Re;
  state.alpha = alpha * std::numbers::pi / 180.0;
  state.referenceMach = Mach;
  state.reynoldsType = reType;
  state.machType = maType;
  state.controlByAlpha = true;
  state.qinf = 1.0;
  state.viscous = bViscous;

  acrit = NCrit;
  boundaryLayerWorkflow.lattice.top.transitionLocation = XtrTop;
  boundaryLayerWorkflow.lattice.bottom.transitionLocation = XtrBot;

  if (Mach > 0.000001) {
    setMach();
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
 *			wake:	  n+1 < i < n+foil.wake_shape.n
 *-------------------------------------------------------------------- */
/** ------------------------------------------------------
 *         calculates source panel influence coefficient
 * 	   matrix for current airfoil and wake geometry.
 * ------------------------------------------------------ */
bool XFoil::qdcalc() {
  // TRACE("calculating source influence matrix ...\n");
  Logger::instance().write("   Calculating source influence matrix ...\n");
  const int point_count = foil.foil_shape.n;

  if (!hasAirfoilInfluenceMatrix()) {
    //----- calculate source influence matrix for airfoil surface
    aerodynamicCache.bij.block(0, 0, point_count + 1, point_count) =
        aerodynamicCache.psi_gamma_lu.solve(aerodynamicCache.bij.block(0, 0, point_count + 1, point_count)).eval();

    //------- store resulting dgam/dsig = dqtan/dsig vector
    aerodynamicCache.dij.block(0, 0, point_count, point_count) =
        aerodynamicCache.bij.block(0, 0, point_count, point_count);
  }

  //---- set up coefficient matrix of dpsi/dm on airfoil surface
  for (int i = 0; i < point_count; i++) {
    PsiResult psi_result =
        pswlin(foil, i, foil.foil_shape.points.col(i),
               foil.foil_shape.normal_vector.col(i), foil.wake_shape.angle_panel);
    aerodynamicCache.bij.row(i).segment(point_count, foil.wake_shape.n) = -psi_result.dzdm.segment(point_count, foil.wake_shape.n).transpose();
  }

  //---- set up kutta condition (no direct source influence)

  aerodynamicCache.bij.row(point_count).segment(point_count, foil.wake_shape.n).setZero();

  //---- multiply by inverse of factored dpsi/dgam matrix
  aerodynamicCache.bij.block(0, point_count, point_count + 1, foil.wake_shape.n) =
      aerodynamicCache.psi_gamma_lu.solve(aerodynamicCache.bij.block(0, point_count, point_count + 1, foil.wake_shape.n)).eval();
  //---- set the source influence matrix for the wake sources
  aerodynamicCache.dij.block(0, point_count, point_count, foil.wake_shape.n) = aerodynamicCache.bij.block(0, point_count, point_count, foil.wake_shape.n);

  //**** now we need to calculate the influence of sources on the wake
  // velocities

  //---- calculate dqtan/dgam and dqtan/dsig at the wake points
  MatrixXd cij = MatrixXd::Zero(foil.wake_shape.n, point_count);
  for (int i = point_count; i < point_count + foil.wake_shape.n; i++) {
    int iw = i - point_count;
    //------ airfoil contribution at wake panel node
    PsiResult psi_result =
        psilin(foil, i, foil.wake_shape.points.col(i),
               foil.wake_shape.normal_vector.col(i), true, aerodynamicCache.gamu,
               surface_vortex, analysis_state_.alpha, analysis_state_.qinf,
               foil.wake_shape.angle_panel);
    cij.row(iw) = psi_result.dqdg.head(point_count).transpose();
    aerodynamicCache.dij.row(i).head(point_count) = psi_result.dqdm.head(point_count).transpose();
    //------ wake contribution
    psi_result =
        pswlin(foil, i, foil.wake_shape.points.col(i),
               foil.wake_shape.normal_vector.col(i), foil.wake_shape.angle_panel);
    aerodynamicCache.dij.row(i).segment(point_count, foil.wake_shape.n) = psi_result.dqdm.segment(point_count, foil.wake_shape.n).transpose();
  }

  //---- add on effect of all sources on airfoil vorticity which effects wake
  // qtan
  aerodynamicCache.dij.block(point_count, 0, foil.wake_shape.n, point_count) += cij * aerodynamicCache.dij.topLeftCorner(point_count, point_count);

  aerodynamicCache.dij.block(point_count, point_count, foil.wake_shape.n, foil.wake_shape.n) += cij * aerodynamicCache.bij.block(0, point_count, point_count, foil.wake_shape.n);

  //---- make sure first wake point has same velocity as trailing edge
  aerodynamicCache.dij.row(point_count) = aerodynamicCache.dij.row(point_count - 1);

  return true;
}


/** -------------------------------------------------------------
 *     sets panel viscous tangential velocity from viscous ue
 * -------------------------------------------------------------- */
VectorXd XFoil::qvfue(const VectorXd& base_qvis,
                      const SidePair<BoundaryLayerLattice>& lattice) const {
  VectorXd updated_qvis = base_qvis;
  for (int is = 1; is <= 2; is++) {
    const auto& panelInfluenceFactor_side = lattice.get(is).panelInfluenceFactor;
    const auto& edgeVelocity_side = lattice.get(is).profiles.edgeVelocity;
    const int limit = lattice.get(is).stationCount - 1;
    for (int ibl = 0; ibl < limit; ++ibl) {
      int i = lattice.get(is).stationToPanel[ibl];
      updated_qvis[i] = panelInfluenceFactor_side[ibl] * edgeVelocity_side[ibl];
    }
  }

  return updated_qvis;
}


/** ---------------------------------------------------------------
 *      sets inviscid tangential velocity for alpha = 0, 90
 *      on wake due to freestream and airfoil surface vorticity.
 * --------------------------------------------------------------- */
Matrix2Xd XFoil::qwcalc(const Foil& foil, const Matrix2Xd& base_qinvu,
                        const Matrix2Xd& gamu,
                        const Matrix2Xd& surface_vortex, double alpha,
                        double qinf) const {
  const int point_count = foil.foil_shape.n;
  Matrix2Xd updated_qinvu = base_qinvu;

  if (point_count >= 1 && point_count < updated_qinvu.cols()) {
    updated_qinvu.col(point_count) = updated_qinvu.col(point_count - 1);
  }

  for (int i = point_count + 1; i < point_count + foil.wake_shape.n; i++) {
    updated_qinvu.col(i) =
        psilin(foil, i, foil.wake_shape.points.col(i),
               foil.wake_shape.normal_vector.col(i), false, gamu,
               surface_vortex, alpha, qinf, foil.wake_shape.angle_panel)
            .qtan;
  }

  return updated_qinvu;
}

XFoil::EdgeVelocitySwapResult XFoil::swapEdgeVelocities(
    const SidePair<VectorXd> &usav) const {
  EdgeVelocitySwapResult result;
  result.swappedUsav = usav;
  result.restoredUedg.top = boundaryLayerWorkflow.lattice.top.profiles.edgeVelocity;
  result.restoredUedg.bottom = boundaryLayerWorkflow.lattice.bottom.profiles.edgeVelocity;
  for (int is = 1; is <= 2; ++is) {
    for (int ibl = 0; ibl < boundaryLayerWorkflow.lattice.get(is).stationCount - 1; ++ibl) {
      result.swappedUsav.get(is)[ibl] = boundaryLayerWorkflow.lattice.get(is).profiles.edgeVelocity[ibl];
      result.restoredUedg.get(is)[ibl] = usav.get(is)[ibl];
    }
  }
  return result;
}

double XFoil::computeAcChange(double clnew, double cl_current,
                              double cl_target, double cl_ac, double cl_a,
                              double cl_ms) const {
  if (analysis_state_.controlByAlpha) {
    return (clnew - cl_current) /
           (1.0 - cl_ac -
            cl_ms * 2.0 * analysis_state_.currentMach * minf_cl);
  }
  return (clnew - cl_target) / (0.0 - cl_ac - cl_a);
}

double XFoil::rlxCalc(double dac) const {
  //---- max allowable alpha changes per iteration
  const double dtor = std::numbers::pi / 180.0;
  double dalmax = 0.5 * dtor;
  double dalmin = -0.5 * dtor;
  //---- max allowable cl change per iteration
  double dclmax = 0.5;
  double dclmin = -0.5;

  if (analysis_state_.machType != MachType::CONSTANT)
    dclmin = std::max(-0.5, -0.9 * aero_coeffs_.cl);

  auto clampRelaxationForGlobalChange = [&](double relaxation, double dac,
                                         double lower, double upper) {
    if (dac == 0.0)
      return relaxation;
    if (relaxation * dac > upper)
      relaxation = upper / dac;
    if (relaxation * dac < lower)
      relaxation = lower / dac;
    return relaxation;
  };

  if (analysis_state_.controlByAlpha)
    return clampRelaxationForGlobalChange(1.0, dac, dclmin, dclmax);

  return clampRelaxationForGlobalChange(1.0, dac, dalmin, dalmax);
}

XFoil::UpdateResult XFoil::update(const XFoil::Matrix3x2dVector& vdel) const {
  //------------------------------------------------------------------
  //      adds on newton deltas to boundary layer variables.
  //      checks for excessive changes and underrelaxes if necessary.
  //      calculates max and rms changes.
  //      also calculates the change in the global variable "ac".
  //        if controlByAlpha=true , "ac" is cl
  //        if controlByAlpha=false, "ac" is alpha
  //------------------------------------------------------------------

  UpdateResult result;
  result.analysis_state = analysis_state_;
  result.aero_coeffs = aero_coeffs_;
  result.profiles.top.skinFrictionCoeffHistory =
      boundaryLayerWorkflow.lattice.top.profiles.skinFrictionCoeffHistory;
  result.profiles.bottom.skinFrictionCoeffHistory =
      boundaryLayerWorkflow.lattice.bottom.profiles.skinFrictionCoeffHistory;
  const double gamma = 1.4;

  BoundaryLayerMarcher marcher;
  result.hstinv =
      (gamma - 1) * MathUtil::pow(analysis_state_.currentMach / analysis_state_.qinf, 2) /
      (1.0 + 0.5 * (gamma - 1) * analysis_state_.currentMach * analysis_state_.currentMach);

  //--- calculate new ue distribution and tangential velocities
  const auto ue_distribution =
      marcher.computeNewUeDistribution(boundaryLayerWorkflow, *this, vdel);
  const auto cl_contributions =
      marcher.computeClFromEdgeVelocityDistribution(boundaryLayerWorkflow,
                                                    *this, ue_distribution);

  const double cl_target =
      analysis_state_.controlByAlpha ? aero_coeffs_.cl : analysis_state_.clspec;
  double dac = computeAcChange(cl_contributions.cl, aero_coeffs_.cl, cl_target,
                               cl_contributions.cl_ac, cl_contributions.cl_a,
                               cl_contributions.cl_ms);

  double rlx = rlxCalc(dac);

  double rmsbl = 0.0;
  double rmxbl = 0.0;
  const double dhi = 1.5;
  const double dlo = -0.5;

  SidePair<BoundaryLayerDelta> deltas;
  SidePair<BoundaryLayerMetrics> metrics;
  for (int side = 1; side <= 2; ++side) {
    deltas.get(side) =
        marcher.buildBoundaryLayerDelta(
            boundaryLayerWorkflow, side, ue_distribution.unew.get(side),
            ue_distribution.u_ac.get(side), dac, *this, vdel);
    metrics.get(side) =
        marcher.evaluateSegmentRelaxation(
            boundaryLayerWorkflow, side, deltas.get(side), dhi, dlo, rlx);
    rmsbl += metrics.get(side).rmsContribution;
    rmxbl = std::max(rmxbl, metrics.get(side).maxChange);
    result.profiles.get(side) =
        marcher.applyBoundaryLayerDelta(
            boundaryLayerWorkflow, side, deltas.get(side), rlx, result.hstinv,
            gamma - 1);
  }

  rmsbl = sqrt(rmsbl / (4.0 * double(boundaryLayerWorkflow.lattice.top.stationCount + boundaryLayerWorkflow.lattice.bottom.stationCount)));

  if (analysis_state_.controlByAlpha)
    result.aero_coeffs.cl = aero_coeffs_.cl + rlx * dac;
  else
    result.analysis_state.alpha =
        analysis_state_.alpha + rlx * dac;

  //--- equate upper wake arrays to lower wake arrays
  for (int kbl = 1; kbl <= boundaryLayerWorkflow.lattice.bottom.stationCount - (boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex + 1); kbl++) {
    const int top_index = boundaryLayerWorkflow.lattice.top.trailingEdgeIndex + kbl;
    const int bottom_index =
        boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex + kbl;
    result.profiles.top.skinFrictionCoeff[top_index] =
        result.profiles.bottom.skinFrictionCoeff[bottom_index];
    result.profiles.top.momentumThickness[top_index] =
        result.profiles.bottom.momentumThickness[bottom_index];
    result.profiles.top.displacementThickness[top_index] =
        result.profiles.bottom.displacementThickness[bottom_index];
    result.profiles.top.edgeVelocity[top_index] =
        result.profiles.bottom.edgeVelocity[bottom_index];
    result.profiles.top.skinFrictionCoeffHistory[top_index] =
        result.profiles.bottom.skinFrictionCoeffHistory[bottom_index];
  }

  result.rlx = rlx;
  result.rmsbl = rmsbl;
  result.rmxbl = rmxbl;
  result.dac = dac;
  return result;
}

void XFoil::applyUpdateResult(UpdateResult result) {
  boundaryLayerWorkflow.blCompressibility.hstinv = result.hstinv;
  analysis_state_ = std::move(result.analysis_state);
  aero_coeffs_ = std::move(result.aero_coeffs);
  boundaryLayerWorkflow.lattice.top.profiles =
      std::move(result.profiles.top);
  boundaryLayerWorkflow.lattice.bottom.profiles =
      std::move(result.profiles.bottom);
}


bool XFoil::viscal() {
  ////--------------------------------------
  //     converges viscous operating point
  ////--------------------------------------
  const int point_count = foil.foil_shape.n;
  const int total_nodes_with_wake = point_count + foil.wake_shape.n;

  //---- calculate wake trajectory from current inviscid solution if necessary
  xyWake();

  //	---- set velocities on wake from airfoil vorticity for alpha=0, 90
  aerodynamicCache.qinvu =
      qwcalc(foil, aerodynamicCache.qinvu, aerodynamicCache.gamu,
             surface_vortex, analysis_state_.alpha, analysis_state_.qinf);

  //	---- set velocities on airfoil and wake for initial alpha
  qinv_matrix =
      InviscidSolver::qiset(analysis_state_.alpha, aerodynamicCache.qinvu);
  

  if (!hasPanelMap()) {
    if (isBLInitialized())
      surface_vortex = gamqv();

    //	----- locate stagnation point arc length position and panel index
    const auto stagnation = boundaryLayerWorkflow.geometry.stfind(
        surface_vortex, foil.foil_shape.spline_length);
    if (!stagnation.found) {
      Logger::instance().write("stfind: Stagnation point not found. Continuing ...\n");
    }
    boundaryLayerWorkflow.stagnationIndex = stagnation.stagnationIndex;
    this->stagnation = stagnation;
    boundaryLayerWorkflow.stagnationSst = stagnation.sst;

    //	----- set  bl position -> panel position  pointers
    boundaryLayerWorkflow.geometry.iblpan(point_count, foil.wake_shape.n);

    //	----- calculate surface arc length array for current stagnation point
    // location
    boundaryLayerWorkflow.geometry.xicalc(foil);

    //	----- set  bl position -> system line  pointers
    boundaryLayerWorkflow.geometry.iblsys(boundaryLayerWorkflow.nsys);
  }

  //	---- set inviscid bl edge velocity inviscidEdgeVelocityMatrix from qinv_matrix
  {
    const auto inviscid_edge_velocity = boundaryLayerWorkflow.geometry.uicalc(qinv_matrix);
    boundaryLayerWorkflow.lattice.top.inviscidEdgeVelocityMatrix = inviscid_edge_velocity.top;
    boundaryLayerWorkflow.lattice.bottom.inviscidEdgeVelocityMatrix = inviscid_edge_velocity.bottom;
  }

  if (!isBLInitialized()) {
    //	----- set initial ue from inviscid ue
    for (int ibl = 0; ibl < boundaryLayerWorkflow.lattice.top.stationCount - 1; ibl++) {
      boundaryLayerWorkflow.lattice.top.profiles.edgeVelocity[ibl] =
          boundaryLayerWorkflow.lattice.top.inviscidEdgeVelocityMatrix(0, ibl);
    }
    for (int ibl = 0; ibl < boundaryLayerWorkflow.lattice.bottom.stationCount - 1; ibl++) {
      boundaryLayerWorkflow.lattice.bottom.profiles.edgeVelocity[ibl] =
          boundaryLayerWorkflow.lattice.bottom.inviscidEdgeVelocityMatrix(0, ibl);
    }
  }

  if (hasConvergedSolution()) {
    //	----- set correct cl if converged point exists
    qvis = qvfue(qvis, boundaryLayerWorkflow.lattice);

    if (analysis_state_.viscous) {
      cpv = InviscidSolver::cpcalc(total_nodes_with_wake, qvis, analysis_state_.qinf,
                   analysis_state_.currentMach);
      cpi = InviscidSolver::cpcalc(total_nodes_with_wake, qinv_matrix.row(0).transpose(),
                   analysis_state_.qinf,
                   analysis_state_.currentMach);
    } else {
      cpi = InviscidSolver::cpcalc(point_count, qinv_matrix.row(0).transpose(),
                   analysis_state_.qinf,
                   analysis_state_.currentMach);
    }

    const auto cl_result = clcalc(cmref);
    applyClComputation(cl_result);
    aero_coeffs_.cd = cdcalc();
  }

  //	---- set up source influence matrix if it doesn't exist
  if (!hasAirfoilInfluenceMatrix() || !hasWakeInfluenceMatrix())
    qdcalc();

  return true;
}


XFoil::ViscalEndResult XFoil::ViscalEnd() {
  ViscalEndResult result;
  const int total_nodes_with_wake = foil.foil_shape.n + foil.wake_shape.n;
  result.inviscidCp = InviscidSolver::cpcalc(total_nodes_with_wake, qinv_matrix.row(0).transpose(),
                             analysis_state_.qinf,
                             analysis_state_.currentMach);
  result.viscousCp = InviscidSolver::cpcalc(total_nodes_with_wake, qvis, analysis_state_.qinf,
                             analysis_state_.currentMach);
  return result;
}


bool XFoil::ViscousIter() {
  //	Performs one iteration
  std::stringstream ss;
  double eps1 = 0.0001;

  auto setbl_output = setbl(SidePairRef<const BoundaryLayerSideProfiles>{
      boundaryLayerWorkflow.lattice.top.profiles,
      boundaryLayerWorkflow.lattice.bottom.profiles});
  setbl_output.applyToXFoil(
      *this); //	------ fill newton system for bl variables

  Blsolve solver;
  SidePair<int> ivte{
      boundaryLayerWorkflow.lattice.top.stationToSystem[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex],
      boundaryLayerWorkflow.lattice.bottom.stationToSystem[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex]
  };
  auto result =
      solver.solve(boundaryLayerWorkflow.nsys, ivte, VAccel(),
                   setbl_output.bl_newton_system);

  const auto update_result = update(result.vdel);
  applyUpdateResult(update_result); //	------ update bl variables
  const double rmsbl = update_result.rmsbl;

  if (analysis_state_.controlByAlpha) { //	------- set new freestream mach, re from new cl
    minf_cl = getActualMach(aero_coeffs_.cl, analysis_state_.machType);
    reinf_cl = getActualReynolds(aero_coeffs_.cl, analysis_state_.reynoldsType);
  } else { //	------- set new inviscid speeds qinv_matrix and inviscidEdgeVelocityMatrix for new alpha
    qinv_matrix =
        InviscidSolver::qiset(analysis_state_.alpha, aerodynamicCache.qinvu);
    const auto inviscid_edge_velocity = boundaryLayerWorkflow.geometry.uicalc(qinv_matrix);
    boundaryLayerWorkflow.lattice.top.inviscidEdgeVelocityMatrix = inviscid_edge_velocity.top;
    boundaryLayerWorkflow.lattice.bottom.inviscidEdgeVelocityMatrix = inviscid_edge_velocity.bottom;
  }

  qvis = qvfue(qvis, boundaryLayerWorkflow.lattice);  //	------ calculate edge velocities qvis(.) from edgeVelocity(..)
  surface_vortex = gamqv();  //	------ set gam distribution from qvis
  boundaryLayerWorkflow.geometry.stmove(surface_vortex, foil.foil_shape.spline_length,
                               foil, qinv_matrix, stagnation,
                               boundaryLayerWorkflow.nsys); //	------ relocate stagnation point

  //	------ set updated cl,cd
  const auto cl_result = clcalc(cmref);
  applyClComputation(cl_result);
  aero_coeffs_.cd = cdcalc();

  if (rmsbl < eps1) {
    avisc = analysis_state_.alpha;
    mvisc = analysis_state_.currentMach;
    Logger::instance().write("----------CONVERGED----------\n\n");
  }

  return true;
}

bool XFoil::isValidFoilAngles(Matrix2Xd points) {
  auto cang = [&](Matrix2Xd points) {
    double max_angle = 0;
    //---- go over each point, calculating corner angle
    for (int i = 1; i < points.cols() - 1; i++) {
      Vector2d delta_former = points.col(i) - points.col(i - 1);
      Vector2d delta_later = points.col(i) - points.col(i + 1);
      double sin = MathUtil::cross2(delta_later, delta_former) / delta_former.norm() / delta_later.norm();
      double delta_angle = asin(sin) * 180.0 / std::numbers::pi;
      max_angle = std::max(fabs(delta_angle), max_angle);
    }
    return max_angle;
  };
  double max_angle = cang(points);
  return max_angle <= kAngleTolerance;
}


bool XFoil::isValidFoilPointSize(Matrix2Xd points) {
  return points.cols() >= 3;
}
