#include "XFoil.h"
#include "simulation/InviscidSolver.hpp"
#include "simulation/psi.hpp"
#include "core/math_util.hpp"
#include "infrastructure/logger.hpp"
#include <algorithm>
#include <cmath>
using Eigen::FullPivLU;
using Eigen::Matrix2d;
using Eigen::Matrix2Xd;
using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::VectorXd;

double XFoil::cdcalc() const {
  if (!(analysis_state_.viscous && isBLInitialized())) {
    return 0.0;
  }

  const double beta =
      std::sqrt(std::max(0.0, 1.0 - analysis_state_.currentMach * analysis_state_.currentMach));
  const double tklam_local =
      MathUtil::pow(analysis_state_.currentMach / (1.0 + beta), 2);

  const double thwake = boundaryLayerWorkflow.lattice.get(2).profiles.momentumThickness[boundaryLayerWorkflow.lattice.bottom.stationCount - 2];
  const double edgeVelocity_bottom = boundaryLayerWorkflow.lattice.bottom.profiles.edgeVelocity[boundaryLayerWorkflow.lattice.bottom.stationCount - 2];
  const double urat = edgeVelocity_bottom / analysis_state_.qinf;
  const double uewake =
      edgeVelocity_bottom * (1.0 - tklam_local) /
      (1.0 - tklam_local * urat * urat);
  const double shwake =
      boundaryLayerWorkflow.lattice.get(2).profiles.displacementThickness[boundaryLayerWorkflow.lattice.bottom.stationCount - 2] /
      boundaryLayerWorkflow.lattice.get(2).profiles.momentumThickness[boundaryLayerWorkflow.lattice.bottom.stationCount - 2];

  const double exponent = 0.5 * (5.0 + shwake);
  const double wake_ratio = uewake / analysis_state_.qinf;
  const double wake_term = std::pow(wake_ratio, exponent);
  return 2.0 * thwake * wake_term;
}

double XFoil::getXcp() const {
  return aero_coeffs_.xcp;
}

XFoil::ClComputation XFoil::clcalc(Vector2d ref) const {

  //-----------------------------------------------------------
  //	   integrates surface pressures to get cl and cm.
  //	   integrates skin friction to get cdf.
  //	   calculates dcl/dalpha for prescribed-cl routines.
  //-----------------------------------------------------------

  ClComputation result;
  double xcp_accumulator = 0.0;

  const auto compressibility = buildCompressibilityParams();
  const Matrix2d rotateMatrix = MathUtil::getRotateMatrix(analysis_state_.alpha);
  const int point_count = foil.foil_shape.n;

  const PressureCoefficientResult cp_first = computePressureCoefficient(
      surface_vortex(0, 0), surface_vortex(1, 0), compressibility);

  double cpg1 = cp_first.cp;
  double cpg1_msq = cp_first.cp_msq;
  double cpg1_alf = cp_first.cp_velocity_derivative;

  for (int i = 0; i < point_count; i++) {
    const int ip = (i + 1) % point_count;
    const PressureCoefficientResult cp_next = computePressureCoefficient(
        surface_vortex(0, ip), surface_vortex(1, ip), compressibility);

    const double cpg2 = cp_next.cp;
    const double cpg2_msq = cp_next.cp_msq;
    const double cpg2_alf = cp_next.cp_velocity_derivative;

    const Vector2d delta = foil.foil_shape.points.col(ip) - foil.foil_shape.points.col(i);
    const Vector2d dpoint = rotateMatrix * delta;
    const double dg = cpg2 - cpg1;

    const Vector2d apoint =
        rotateMatrix * ((foil.foil_shape.points.col(ip) + foil.foil_shape.points.col(i)) / 2 - ref);
    const double ag = 0.5 * (cpg2 + cpg1);

    const double dx_alf = MathUtil::cross2(delta, rotateMatrix.row(0).transpose());
    const double ag_alf = 0.5 * (cpg2_alf + cpg1_alf);
    const double ag_msq = 0.5 * (cpg2_msq + cpg1_msq);

    result.cl += dpoint.x() * ag;
    result.cm -= dpoint.dot(ag * apoint + dg * dpoint / 12.0);

    xcp_accumulator += dpoint.x() * ag *
                       (foil.foil_shape.points.col(ip).x() +
                        foil.foil_shape.points.col(i).x()) /
                       2.0;

    result.cl_alf += dpoint.x() * ag_alf + ag * dx_alf;
    result.cl_msq += dpoint.x() * ag_msq;

    cpg1 = cpg2;
    cpg1_alf = cpg2_alf;
    cpg1_msq = cpg2_msq;
  }

  if (fabs(result.cl) > 0.0)
    result.xcp = xcp_accumulator / result.cl;
  else
    result.xcp = 0.0;

  return result;
}

void XFoil::applyClComputation(const ClComputation &result) {
  aero_coeffs_.cl = result.cl;
  aero_coeffs_.cm = result.cm;
  aero_coeffs_.cl_alf = result.cl_alf;
  aero_coeffs_.cl_msq = result.cl_msq;
  aero_coeffs_.xcp = result.xcp;
}

Matrix2Xd XFoil::gamqv() const {
  const int point_count = foil.foil_shape.n;
  Matrix2Xd updated_surface_vortex(2, point_count);
  for (int i = 0; i < point_count; i++) {
    updated_surface_vortex(0, i) = qvis[i];
    updated_surface_vortex(1, i) = qinv_matrix(1, i);
  }
  return updated_surface_vortex;
}


/** --------------------------------------------------------------
 *     Calculates two surface vorticity (gamma) distributions
 *     for alpha = 0, 90  degrees.  These are superimposed
 *     in specal or speccl for specified alpha or cl.
 *-------------------------------------------------------------- */
FoilAerodynamicCache XFoil::ggcalc() {
  FoilAerodynamicCache cache = aerodynamicCache;

  Logger::instance().write("   Calculating unit vorticity distributions ...\n");

  MatrixXd dpsi_dgam = MatrixXd::Zero(foil.foil_shape.n + 1, foil.foil_shape.n + 1);

  Matrix2Xd psi = Matrix2Xd::Zero(2, foil.foil_shape.n + 1);

  //---- set up matrix system for  psi = psio  on airfoil surface.
  //-    the unknowns are (dgamma)i and dpsio.
  for (int i = 0; i < foil.foil_shape.n; i++) {
    //------ calculate psi and dpsi/dgamma array for current node
    PsiResult psi_result =
        psilin(foil, i, foil.foil_shape.points.col(i),
               foil.foil_shape.normal_vector.col(i), true, cache.gamu,
               surface_vortex, analysis_state_.alpha, analysis_state_.qinf,
               foil.foil_shape.angle_panel);
    //------ dres/dgamma
    dpsi_dgam.row(i).head(foil.foil_shape.n) = psi_result.dzdg;
    cache.bij.row(i).head(foil.foil_shape.n) = -psi_result.dzdm;
  }
  psi.leftCols(foil.foil_shape.n) = -analysis_state_.qinf * Matrix2d{{0.0, 1.0}, {-1.0, 0.0}} * foil.foil_shape.points;
  //------ dres/dpsio
  dpsi_dgam.col(foil.foil_shape.n) = - VectorXd::Ones(foil.foil_shape.n);
  dpsi_dgam.row(foil.foil_shape.n) = VectorXd::Zero(foil.foil_shape.n + 1);
  dpsi_dgam(foil.foil_shape.n, 0) = 1;
  dpsi_dgam(foil.foil_shape.n, foil.foil_shape.n - 1) = 1;

  //---- set up Kutta condition (no direct source influence)
  cache.bij.row(foil.foil_shape.n).head(foil.foil_shape.n) = VectorXd::Zero(foil.foil_shape.n);

  if (foil.edge.sharp) {
    //----- set zero internal velocity in TE corner

    //----- set TE bisector angle
    const double ag1 = atan2(-foil.foil_shape.dpoints_ds.col(0).y(), -foil.foil_shape.dpoints_ds.col(0).x());
    const double ag2 =
        atanc(foil.foil_shape.dpoints_ds.col(foil.foil_shape.n - 1).y(), foil.foil_shape.dpoints_ds.col(foil.foil_shape.n - 1).x(), ag1);
    const double abis = 0.5 * (ag1 + ag2);

    Vector2d bis_vector{cos(abis), sin(abis)};

    //----- minimum panel length adjacent to TE
    const double dsmin = std::min((foil.foil_shape.points.col(0) - foil.foil_shape.points.col(1)).norm(),
                                  (foil.foil_shape.points.col(foil.foil_shape.n - 1) - foil.foil_shape.points.col(foil.foil_shape.n - 2)).norm());

    //---- distance of internal control point ahead of sharp TE
    //-    (fraction of smaller panel length adjacent to TE)
    const double bwt = 0.1;

    //----- control point on bisector just ahead of TE point
    const Vector2d bis = foil.edge.point_te - bwt * dsmin * bis_vector;
    const Vector2d normal_bis{-bis_vector.y(), bis_vector.x()};

    //----- set velocity component along bisector line
    PsiResult psi_result =
        psilin(foil, -1, bis, normal_bis, true, cache.gamu,
               surface_vortex, analysis_state_.alpha, analysis_state_.qinf,
               foil.foil_shape.angle_panel);

    //----- dres/dgamma
    dpsi_dgam.row(foil.foil_shape.n - 1).head(foil.foil_shape.n) = psi_result.dzdg;
    //----- -dres/dmass
    cache.bij.row(foil.foil_shape.n - 1).head(foil.foil_shape.n) = -psi_result.dzdm;

    //----- dres/dpsio

    //----- -dres/duinf -dres/dvinf
    psi.col(foil.foil_shape.n - 1) = -bis_vector;
  }

  //---- lu-factor coefficient matrix aij
  cache.psi_gamma_lu = FullPivLU<MatrixXd>(dpsi_dgam);
  cache.gamu.resize(2, foil.foil_shape.n + 1);
  cache.gamu.row(0) = cache.psi_gamma_lu.solve(psi.row(0).transpose()).transpose();
  cache.gamu.row(1) = cache.psi_gamma_lu.solve(psi.row(1).transpose()).transpose();
  
  //---- set inviscid alpha=0,90 surface speeds for this geometry
  cache.qinvu.leftCols(foil.foil_shape.n + 1) = cache.gamu;

  return cache;
}

bool XFoil::specal() {
  return InviscidSolver::specal(*this);
}

bool XFoil::speccl() {
  return InviscidSolver::speccl(*this);
}
