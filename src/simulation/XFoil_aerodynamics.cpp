#include "XFoil.h"
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <utility>
using namespace Eigen;

namespace {
struct AerodynamicsState {
  double xcp = 0.0;
};

std::unordered_map<const XFoil*, AerodynamicsState> g_aerodynamics_state;

AerodynamicsState& ensureAerodynamicsState(const XFoil* foil) {
  return g_aerodynamics_state[foil];
}
}  // namespace

double cross2(const Eigen::Vector2d &a, const Eigen::Vector2d &b);


/**
 * @brief find index and angle of max foil curv diff
 *
 * @param x foil x parameters
 * @param y foil y parameters
 * @param n foil plot size
 * @return PairIndex index: index of max angle diff, value: angle of max angle
 * diff[degree]
 */
// moved to XFoil_geometry.cpp: cang()

double XFoil::cdcalc() const {
  if (!(lvisc && lblini)) {
    return 0.0;
  }

  const double beta = std::sqrt(std::max(0.0, 1.0 - minf * minf));
  const double tklam_local = MathUtil::pow(minf / (1.0 + beta), 2);

  const double thwake = thet.get(2)[nbl.bottom - 2];
  const double uedg_bottom = uedg.bottom[nbl.bottom - 2];
  const double urat = uedg_bottom / qinf;
  const double uewake =
      uedg_bottom * (1.0 - tklam_local) /
      (1.0 - tklam_local * urat * urat);
  const double shwake =
      dstr.get(2)[nbl.bottom - 2] / thet.get(2)[nbl.bottom - 2];

  const double exponent = 0.5 * (5.0 + shwake);
  const double wake_ratio = uewake / qinf;
  const double wake_term = std::pow(wake_ratio, exponent);
  return 2.0 * thwake * wake_term;
}

double XFoil::getXcp() const {
  auto it = g_aerodynamics_state.find(this);
  if (it == g_aerodynamics_state.end()) {
    return 0.0;
  }
  return it->second.xcp;
}

void ClearAerodynamicsState(const XFoil& foil) {
  g_aerodynamics_state.erase(&foil);
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
  const Matrix2d rotateMatrix = buildBodyToFreestreamRotation();

  const PressureCoefficientResult cp_first = computePressureCoefficient(
      surface_vortex(0, 0), surface_vortex(1, 0), compressibility);

  double cpg1 = cp_first.cp;
  double cpg1_msq = cp_first.cp_msq;
  double cpg1_alf = cp_first.cp_velocity_derivative;

  for (int i = 0; i < n; i++) {
    const int ip = (i + 1) % n;
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

    const double dx_alf = cross2(delta, rotateMatrix.row(0));
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
  cl = result.cl;
  cm = result.cm;
  cl_alf = result.cl_alf;
  cl_msq = result.cl_msq;

  auto &aero = ensureAerodynamicsState(this);
  aero.xcp = result.xcp;
}


// moved to XFoil_init.cpp: comset()

/** ---------------------------------------------
 *      sets compressible cp from speed.
 * ---------------------------------------------- */
VectorXd XFoil::cpcalc(int n, VectorXd q, double qinf, double minf) {
  VectorXd cp = VectorXd::Zero(n);
  bool denneg = false;
  double beta, bfac;

  beta = sqrt(1.0 - MathUtil::pow(minf, 2));
  bfac = 0.5 * MathUtil::pow(minf, 2) / (1.0 + beta);

  for (int i = 0; i < n; i++) {
    const double cpinc = 1.0 - (q[i] / qinf) * (q[i] / qinf);
    const double den = beta + bfac * cpinc;
    cp[i] = cpinc / den;
    if (den <= 0.0)
      denneg = true;
  }

  if (denneg) {
    writeString("CpCalc: local speed too larger\n Compressibility corrections "
                "invalid\n");
  }

  return cp;
}


/** Laminar dissipation function  ( 2 cd/h* )     (from Falkner-Skan)*/
// moved to XFoil_boundary.cpp: dslim()

Matrix2Xd XFoil::gamqv() const {
  Matrix2Xd updated_surface_vortex(2, n);
  for (int i = 0; i < n; i++) {
    updated_surface_vortex(0, i) = qvis[i];
    updated_surface_vortex(1, i) = qinv_a[i];
  }
  return updated_surface_vortex;
}


/** --------------------------------------------------------------
 *     Calculates two surface vorticity (gamma) distributions
 *     for alpha = 0, 90  degrees.  These are superimposed
 *     in specal or speccl for specified alpha or cl.
 *-------------------------------------------------------------- */
bool XFoil::ggcalc() {

  double res;

  writeString("   Calculating unit vorticity distributions ...\n");

  MatrixXd dpsi_dgam = MatrixXd::Zero(n + 1, n + 1);

  Matrix2Xd psi = Matrix2Xd::Zero(2, n + 1);

  //---- set up matrix system for  psi = psio  on airfoil surface.
  //-    the unknowns are (dgamma)i and dpsio.
  for (int i = 0; i < n; i++) {
    //------ calculate psi and dpsi/dgamma array for current node
    PsiResult psi_result = psilin(
        foil.foil_shape.points, i, foil.foil_shape.points.col(i),
        foil.foil_shape.normal_vector.col(i), true, spline_length, n, gamu,
        surface_vortex, alfa, qinf);

    const Vector2d res = qinf * Vector2d{foil.foil_shape.points.col(i).y(),
                                         -foil.foil_shape.points.col(i).x()};

    //------ dres/dgamma
    dpsi_dgam.row(i).head(n) = psi_result.dzdg.head(n);
    bij.row(i).head(n) = -psi_result.dzdm.head(n);

    //------ dres/dpsio
    dpsi_dgam(i, n) = -1.0;

    psi.col(i) = -res;
  }

  //---- set Kutta condition
  //-    res = gam(1) + gam[n]
  res = 0.0;

  dpsi_dgam.row(n).head(n + 1) = VectorXd::Zero(n + 1);
  dpsi_dgam(n, 0) = 1;
  dpsi_dgam(n, n - 1) = 1;

  psi.col(n).x() = -res;
  psi.col(n).y() = -res;

  //---- set up Kutta condition (no direct source influence)
  bij.row(n).head(n) = VectorXd::Zero(n);

  if (sharp) {
    //----- set zero internal velocity in TE corner

    //----- set TE bisector angle
    const double ag1 = atan2(-foil.foil_shape.dpoints_ds.col(0).y(), -foil.foil_shape.dpoints_ds.col(0).x());
    const double ag2 =
        atanc(foil.foil_shape.dpoints_ds.col(n - 1).y(), foil.foil_shape.dpoints_ds.col(n - 1).x(), ag1);
    const double abis = 0.5 * (ag1 + ag2);

    Vector2d bis_vector{cos(abis), sin(abis)};

    //----- minimum panel length adjacent to TE
    const double dsmin = std::min((foil.foil_shape.points.col(0) - foil.foil_shape.points.col(1)).norm(),
                                  (foil.foil_shape.points.col(foil.foil_shape.n - 1) - foil.foil_shape.points.col(foil.foil_shape.n - 2)).norm());

    //---- distance of internal control point ahead of sharp TE
    //-    (fraction of smaller panel length adjacent to TE)
    const double bwt = 0.1;

    //----- control point on bisector just ahead of TE point
    const Vector2d bis = foil.edge_data.point_te - bwt * dsmin * bis_vector;
    const Vector2d normal_bis{-bis_vector.y(), bis_vector.x()};

    //----- set velocity component along bisector line
    PsiResult psi_result = psilin(foil.foil_shape.points, -1, bis, normal_bis,
                                  true, spline_length, n, gamu,
                                  surface_vortex, alfa, qinf);

    //----- dres/dgamma
    dpsi_dgam.row(n - 1).head(n) = psi_result.dzdg.head(n);
    //----- -dres/dmass
    bij.row(n - 1).head(n) = -psi_result.dzdm.head(n);

    //----- dres/dpsio
    dpsi_dgam(n - 1, n);

    //----- -dres/duinf -dres/dvinf
    psi.col(n - 1) = -bis_vector;
  }

  //---- lu-factor coefficient matrix aij
  psi_gamma_lu = FullPivLU<MatrixXd>(dpsi_dgam);
  lqaij = true;
  VectorXd gamu_temp(n + 1);
  //---- solve system for the two vorticity distributions

  gamu_temp = psi_gamma_lu.solve(psi.row(0).transpose());

  for (int iu = 0; iu <= n; iu++) {
    gamu.col(iu).x() = gamu_temp[iu];
  }

  gamu_temp = psi_gamma_lu.solve(psi.row(1).transpose());
  for (int iu = 0; iu <= n; iu++) {
    gamu.col(iu).y() = gamu_temp[iu];
  }

  //---- set inviscid alpha=0,90 surface speeds for this geometry
  for (int i = 0; i <= n; i++) {
    qinvu.col(i) = gamu.col(i);
  }

  lgamu = true;

  return true;
}


// moved to XFoil_init.cpp: setMach()

/** returns the absolute value of "a" x sign(b) */
// moved to XFoil_geometry.cpp: sign()

/**
 *      Converges to specified alpha.
 */
namespace {
enum class SpecTarget { AngleOfAttack, LiftCoefficient };

void updateSurfaceVortexFromGamu(XFoil &foil) {
  Matrix2d rotateMatrix =
      Matrix2d{{cos(foil.alfa), sin(foil.alfa)}, {-sin(foil.alfa), cos(foil.alfa)}};

  for (int i = 0; i < foil.n; i++) {
    foil.surface_vortex(0, i) = rotateMatrix.row(0).dot(foil.gamu.col(i));
    foil.surface_vortex(1, i) = rotateMatrix.row(1).dot(foil.gamu.col(i));
  }
}

bool specConverge(XFoil &foil, SpecTarget target) {
  // Ensure unit vorticity distributions are available.
  if (!foil.lgamu || !foil.lqaij)
    foil.ggcalc();

  updateSurfaceVortexFromGamu(foil);

  auto applyQiset = [&foil]() {
    auto qiset_result = foil.qiset();
    foil.qinv = std::move(qiset_result.qinv);
    foil.qinv_a = std::move(qiset_result.qinv_a);
  };

  if (target == SpecTarget::AngleOfAttack) {
    foil.tecalc();
    applyQiset();
  } else {
    foil.minf_cl = foil.getActualMach(foil.clspec, foil.mach_type);
    foil.reinf_cl = foil.getActualReynolds(foil.clspec, foil.reynolds_type);
    foil.comset();
  }

  foil.applyClComputation(foil.clcalc(foil.cmref));

  bool bConv = false;

  if (target == SpecTarget::AngleOfAttack) {
    double clm = 1.0;
    double minf_clm = foil.getActualMach(clm, foil.mach_type);

    for (int itcl = 1; itcl <= 20; itcl++) {
      const double msq_clm = 2.0 * foil.minf * minf_clm;
      const double dclm = (foil.cl - clm) / (1.0 - foil.cl_msq * msq_clm);

      const double clm1 = clm;
      foil.rlx = 1.0;

      //------ under-relaxation loop to avoid driving m(cl) above 1
      for (int irlx = 1; irlx <= 12; irlx++) {
        clm = clm1 + foil.rlx * dclm;

        //-------- set new freestream mach m(clm)
        minf_clm = foil.getActualMach(clm, foil.mach_type);

        //-------- if mach is ok, go do next newton iteration
        // FIXME double型の==比較
        if (foil.mach_type == XFoil::MachType::CONSTANT || foil.minf == 0.0 ||
            minf_clm != 0.0)
          break;

        foil.rlx = 0.5 * foil.rlx;
      }

      //------ set new cl(m)
      foil.applyClComputation(foil.clcalc(foil.cmref));

      if (fabs(dclm) <= 1.0e-6) {
        bConv = true;
        break;
      }
    }

    if (!bConv) {
      foil.writeString("Specal:  MInf convergence failed\n");
      return false;
    }

    //---- set final mach, cl, cp distributions, and hinge moment
    foil.minf_cl = foil.getActualMach(foil.cl, foil.mach_type);
    foil.reinf_cl = foil.getActualReynolds(foil.cl, foil.reynolds_type);
    foil.comset();
    applyQiset();
    foil.applyClComputation(foil.clcalc(foil.cmref));

    foil.cpi = foil.cpcalc(foil.n, foil.qinv, foil.qinf, foil.minf);
    if (foil.lvisc) {
      foil.cpv = foil.cpcalc(foil.n + foil.nw, foil.qvis, foil.qinf, foil.minf);
      foil.cpi = foil.cpcalc(foil.n + foil.nw, foil.qinv, foil.qinf, foil.minf);
    } else
      foil.cpi = foil.cpcalc(foil.n, foil.qinv, foil.qinf, foil.minf);

    for (int i = 0; i < foil.n; i++) {
      foil.qgamm[i] = foil.surface_vortex(0, i);
    }

    return true;
  }

  for (int ital = 1; ital <= 20; ital++) {
    const double dalfa = (foil.clspec - foil.cl) / foil.cl_alf;
    foil.rlx = 1.0;

    foil.alfa = foil.alfa + foil.rlx * dalfa;

    updateSurfaceVortexFromGamu(foil);

    //------ set new cl(alpha)
    foil.applyClComputation(foil.clcalc(foil.cmref));

    if (fabs(dalfa) <= 1.0e-6) {
      bConv = true;
      break;
    }
  }
  if (!bConv) {
    foil.writeString("Speccl:  cl convergence failed");
    return false;
  }

  //---- set final surface speed and cp distributions
  foil.tecalc();
  applyQiset();

  if (foil.lvisc) {
    foil.cpv = foil.cpcalc(foil.n + foil.nw, foil.qvis, foil.qinf, foil.minf);
    foil.cpi = foil.cpcalc(foil.n + foil.nw, foil.qinv, foil.qinf, foil.minf);

  } else {
    foil.cpi = foil.cpcalc(foil.n, foil.qinv, foil.qinf, foil.minf);
  }

  return true;
}
}  // namespace

bool XFoil::specal() {
  return specConverge(*this, SpecTarget::AngleOfAttack);
}

bool XFoil::speccl() {
  return specConverge(*this, SpecTarget::LiftCoefficient);
}
