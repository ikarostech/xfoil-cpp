#include "XFoil.h"
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <utility>
using Eigen::FullPivLU;
using Eigen::Matrix2d;
using Eigen::Matrix2Xd;
using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::VectorXd;

namespace {
struct AerodynamicsState {
  double xcp = 0.0;
};

using AerodynamicsStateRegistry =
    std::unordered_map<const XFoil*, AerodynamicsState>;

AerodynamicsStateRegistry& aerodynamicsStateRegistry() {
  static AerodynamicsStateRegistry state;
  return state;
}

AerodynamicsState& ensureAerodynamicsState(const XFoil* xfoil) {
  return aerodynamicsStateRegistry()[xfoil];
}
}  // namespace

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
  const auto& state = aerodynamicsStateRegistry();
  auto it = state.find(this);
  if (it == state.end()) {
    return 0.0;
  }
  return it->second.xcp;
}

void ClearAerodynamicsState(const XFoil& xfoil) {
  aerodynamicsStateRegistry().erase(&xfoil);
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
  cl = result.cl;
  cm = result.cm;
  cl_alf = result.cl_alf;
  cl_msq = result.cl_msq;

  auto &aero = ensureAerodynamicsState(this);
  aero.xcp = result.xcp;
}

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

Matrix2Xd XFoil::gamqv() const {
  const int point_count = foil.foil_shape.n;
  Matrix2Xd updated_surface_vortex(2, point_count);
  for (int i = 0; i < point_count; i++) {
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
  const int point_count = foil.foil_shape.n;

  MatrixXd dpsi_dgam = MatrixXd::Zero(point_count + 1, point_count + 1);

  Matrix2Xd psi = Matrix2Xd::Zero(2, point_count + 1);

  //---- set up matrix system for  psi = psio  on airfoil surface.
  //-    the unknowns are (dgamma)i and dpsio.
  for (int i = 0; i < point_count; i++) {
    //------ calculate psi and dpsi/dgamma array for current node
    PsiResult psi_result =
        psilin(foil, i, foil.foil_shape.points.col(i),
               foil.foil_shape.normal_vector.col(i), true, point_count, gamu,
               surface_vortex, alfa, qinf, apanel, foil.edge.sharp,
               foil.edge.ante, foil.edge.dste, foil.edge.aste);

    const Vector2d res = qinf * Vector2d{foil.foil_shape.points.col(i).y(),
                                         -foil.foil_shape.points.col(i).x()};

    //------ dres/dgamma
    dpsi_dgam.row(i).head(point_count) = psi_result.dzdg.head(point_count);
    bij.row(i).head(point_count) = -psi_result.dzdm.head(point_count);

    //------ dres/dpsio
    dpsi_dgam(i, point_count) = -1.0;

    psi.col(i) = -res;
  }

  //---- set Kutta condition
  //-    res = gam(1) + gam[n]
  res = 0.0;

  dpsi_dgam.row(point_count).head(point_count + 1) = VectorXd::Zero(point_count + 1);
  dpsi_dgam(point_count, 0) = 1;
  dpsi_dgam(point_count, point_count - 1) = 1;

  psi.col(point_count).x() = -res;
  psi.col(point_count).y() = -res;

  //---- set up Kutta condition (no direct source influence)
  bij.row(point_count).head(point_count) = VectorXd::Zero(point_count);

  if (foil.edge.sharp) {
    //----- set zero internal velocity in TE corner

    //----- set TE bisector angle
    const double ag1 = atan2(-foil.foil_shape.dpoints_ds.col(0).y(), -foil.foil_shape.dpoints_ds.col(0).x());
    const double ag2 =
        atanc(foil.foil_shape.dpoints_ds.col(point_count - 1).y(), foil.foil_shape.dpoints_ds.col(point_count - 1).x(), ag1);
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
    PsiResult psi_result = psilin(foil, -1, bis, normal_bis, true, point_count,
                                  gamu, surface_vortex, alfa, qinf, apanel,
                                  foil.edge.sharp, foil.edge.ante, foil.edge.dste,
                                  foil.edge.aste);

    //----- dres/dgamma
    dpsi_dgam.row(point_count - 1).head(point_count) = psi_result.dzdg.head(point_count);
    //----- -dres/dmass
    bij.row(point_count - 1).head(point_count) = -psi_result.dzdm.head(point_count);

    //----- dres/dpsio
    dpsi_dgam(point_count - 1, point_count);

    //----- -dres/duinf -dres/dvinf
    psi.col(point_count - 1) = -bis_vector;
  }

  //---- lu-factor coefficient matrix aij
  psi_gamma_lu = FullPivLU<MatrixXd>(dpsi_dgam);
  lqaij = true;
  VectorXd gamu_temp(point_count + 1);
  //---- solve system for the two vorticity distributions

  gamu_temp = psi_gamma_lu.solve(psi.row(0).transpose());

  for (int iu = 0; iu <= point_count; iu++) {
    gamu.col(iu).x() = gamu_temp[iu];
  }

  gamu_temp = psi_gamma_lu.solve(psi.row(1).transpose());
  for (int iu = 0; iu <= point_count; iu++) {
    gamu.col(iu).y() = gamu_temp[iu];
  }

  //---- set inviscid alpha=0,90 surface speeds for this geometry
  for (int i = 0; i <= point_count; i++) {
    qinvu.col(i) = gamu.col(i);
  }

  lgamu = true;

  return true;
}

/**
 *      Converges to specified alpha.
 */
namespace {
enum class SpecTarget { AngleOfAttack, LiftCoefficient };

void updateSurfaceVortexFromGamu(XFoil &xfoil) {
  Matrix2d rotateMatrix =
      Matrix2d{{cos(xfoil.alfa), sin(xfoil.alfa)}, {-sin(xfoil.alfa), cos(xfoil.alfa)}};

  for (int i = 0; i < xfoil.foil.foil_shape.n; i++) {
    xfoil.surface_vortex(0, i) = rotateMatrix.row(0).dot(xfoil.gamu.col(i));
    xfoil.surface_vortex(1, i) = rotateMatrix.row(1).dot(xfoil.gamu.col(i));
  }
}

bool specConverge(XFoil &xfoil, SpecTarget target) {
  // Ensure unit vorticity distributions are available.
  if (!xfoil.lgamu || !xfoil.lqaij)
    xfoil.ggcalc();

  updateSurfaceVortexFromGamu(xfoil);

  auto applyQiset = [&xfoil]() {
    auto qiset_result = xfoil.qiset();
    xfoil.qinv = std::move(qiset_result.qinv);
    xfoil.qinv_a = std::move(qiset_result.qinv_a);
  };

  if (target == SpecTarget::AngleOfAttack) {
    xfoil.updateTrailingEdgeState();
    applyQiset();
  } else {
    xfoil.minf_cl = xfoil.getActualMach(xfoil.clspec, xfoil.mach_type);
    xfoil.reinf_cl = xfoil.getActualReynolds(xfoil.clspec, xfoil.reynolds_type);
    xfoil.comset();
  }

  xfoil.applyClComputation(xfoil.clcalc(xfoil.cmref));

  bool bConv = false;

  if (target == SpecTarget::AngleOfAttack) {
    double clm = 1.0;
    double minf_clm = xfoil.getActualMach(clm, xfoil.mach_type);

    for (int itcl = 1; itcl <= 20; itcl++) {
      const double msq_clm = 2.0 * xfoil.minf * minf_clm;
      const double dclm = (xfoil.cl - clm) / (1.0 - xfoil.cl_msq * msq_clm);

      const double clm1 = clm;
      xfoil.rlx = 1.0;

      //------ under-relaxation loop to avoid driving m(cl) above 1
      for (int irlx = 1; irlx <= 12; irlx++) {
        clm = clm1 + xfoil.rlx * dclm;

        //-------- set new freestream mach m(clm)
        minf_clm = xfoil.getActualMach(clm, xfoil.mach_type);

        //-------- if mach is ok, go do next newton iteration
        // FIXME double型の==比較
        if (xfoil.mach_type == XFoil::MachType::CONSTANT || xfoil.minf == 0.0 ||
            minf_clm != 0.0)
          break;

        xfoil.rlx = 0.5 * xfoil.rlx;
      }

      //------ set new cl(m)
      xfoil.applyClComputation(xfoil.clcalc(xfoil.cmref));

      if (fabs(dclm) <= 1.0e-6) {
        bConv = true;
        break;
      }
    }

    if (!bConv) {
      xfoil.writeString("Specal:  MInf convergence failed\n");
      return false;
    }

    //---- set final mach, cl, cp distributions, and hinge moment
    applyQiset();
    xfoil.applyClComputation(xfoil.clcalc(xfoil.cmref));

    xfoil.cpi = xfoil.cpcalc(xfoil.foil.foil_shape.n, xfoil.qinv, xfoil.qinf, xfoil.minf);
    if (xfoil.lvisc) {
      xfoil.cpv = xfoil.cpcalc(xfoil.foil.foil_shape.n + xfoil.nw, xfoil.qvis, xfoil.qinf, xfoil.minf);
      xfoil.cpi = xfoil.cpcalc(xfoil.foil.foil_shape.n + xfoil.nw, xfoil.qinv, xfoil.qinf, xfoil.minf);
    } else
      xfoil.cpi = xfoil.cpcalc(xfoil.foil.foil_shape.n, xfoil.qinv, xfoil.qinf, xfoil.minf);

    for (int i = 0; i < xfoil.foil.foil_shape.n; i++) {
      xfoil.qgamm[i] = xfoil.surface_vortex(0, i);
    }

    return true;
  }

  for (int ital = 1; ital <= 20; ital++) {
    const double dalfa = (xfoil.clspec - xfoil.cl) / xfoil.cl_alf;
    xfoil.rlx = 1.0;

    xfoil.alfa = xfoil.alfa + xfoil.rlx * dalfa;

    updateSurfaceVortexFromGamu(xfoil);

    //------ set new cl(alpha)
    xfoil.applyClComputation(xfoil.clcalc(xfoil.cmref));

    if (fabs(dalfa) <= 1.0e-6) {
      bConv = true;
      break;
    }
  }
  if (!bConv) {
    xfoil.writeString("Speccl:  cl convergence failed");
    return false;
  }

  //---- set final surface speed and cp distributions
  xfoil.updateTrailingEdgeState();
  applyQiset();

  if (xfoil.lvisc) {
    xfoil.cpv = xfoil.cpcalc(xfoil.foil.foil_shape.n + xfoil.nw, xfoil.qvis, xfoil.qinf, xfoil.minf);
    xfoil.cpi = xfoil.cpcalc(xfoil.foil.foil_shape.n + xfoil.nw, xfoil.qinv, xfoil.qinf, xfoil.minf);

  } else {
    xfoil.cpi = xfoil.cpcalc(xfoil.foil.foil_shape.n, xfoil.qinv, xfoil.qinf, xfoil.minf);
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
