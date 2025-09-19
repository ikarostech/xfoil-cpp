#include "XFoil.h"
#include <algorithm>
#include <cmath>
#include <unordered_map>
using namespace Eigen;

namespace {
struct AerodynamicsState {
  double xcp = 0.0;
};

std::unordered_map<const XFoil*, AerodynamicsState> g_aerodynamics_state;

AerodynamicsState& ensureAerodynamicsState(XFoil* foil) {
  return g_aerodynamics_state[foil];
}
}  // namespace

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

bool XFoil::cdcalc() {
  // Ensure compressibility parameters reflect the current Mach number
  comset();

  if (!(lvisc && lblini)) {
    cd = 0.0;
    return true;
  }

  //---- set variables at the end of the wake
  double thwake = thet.get(2)[nbl.bottom - 2];
  double urat = uedg.bottom[nbl.bottom - 2] / qinf;
  double uewake =
      uedg.bottom[nbl.bottom - 2] * (1.0 - tklam) /
      (1.0 - tklam * urat * urat);
  double shwake = dstr.get(2)[nbl.bottom - 2] /
                   thet.get(2)[nbl.bottom - 2];

  //---- extrapolate wake to downstream infinity using squire-young relation
  //      (reduces errors of the wake not being long enough)
  cd = 2.0 * thwake * pow((uewake / qinf), (0.5 * (5.0 + shwake)));

  return true;
}

bool XFoil::clcalc(Vector2d ref) {

  //-----------------------------------------------------------
  //	   integrates surface pressures to get cl and cm.
  //	   integrates skin friction to get cdf.
  //	   calculates dcl/dalpha for prescribed-cl routines.
  //-----------------------------------------------------------

  auto& aero = ensureAerodynamicsState(this);
  aero.xcp = 0.0;

  const double beta = sqrt(1.0 - minf * minf);
  const double beta_msq = -0.5 / beta;

  const double bfac = 0.5 * minf * minf / (1.0 + beta);
  const double bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;

  cl = 0.0;
  cm = 0.0;

  cl_alf = 0.0;
  cl_msq = 0.0;

  double cginc = 1.0 - MathUtil::pow((surface_vortex(0, 0) / qinf), 2);
  double cpg1 = cginc / (beta + bfac * cginc);
  double cpg1_msq =
      -cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

  double cpi_gam = -2.0 * surface_vortex(0, 0) / qinf / qinf;
  double cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
  double cpg1_alf = cpc_cpi * cpi_gam * surface_vortex(1, 0);

  for (int i = 0; i < n; i++) {
    int ip = (i + 1) % n;

    cginc = 1.0 - MathUtil::pow((surface_vortex(0, ip) / qinf), 2);
    double cpg2 = cginc / (beta + bfac * cginc);
    double cpg2_msq =
        -cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

    cpi_gam = -2.0 * surface_vortex(0, ip) / qinf / qinf;
    cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
    double cpg2_alf = cpc_cpi * cpi_gam * surface_vortex(1, ip);

    Matrix2d rotateMatrix =
        Matrix2d{{cos(alfa), sin(alfa)}, {-sin(alfa), cos(alfa)}};
    const Vector2d dpoint = rotateMatrix * (points.col(ip) -
                                            points.col(i));
    const double dg = cpg2 - cpg1;

    const Vector2d apoint = rotateMatrix * ((points.col(ip) +
                                             points.col(i)) /
                                                2 -
                                            ref);
    const double ag = 0.5 * (cpg2 + cpg1);

    const double dx_alf = cross2(points.col(ip) -
                                     points.col(i),
                                 rotateMatrix.row(0));
    const double ag_alf = 0.5 * (cpg2_alf + cpg1_alf);
    const double ag_msq = 0.5 * (cpg2_msq + cpg1_msq);

    cl = cl + dpoint.x() * ag;
    cm = cm - dpoint.dot(ag * apoint + dg * dpoint / 12.0);

    aero.xcp += dpoint.x() * ag *
                (points.col(ip).x() +
                 points.col(i).x()) /
                2.0;

    cl_alf = cl_alf + dpoint.x() * ag_alf + ag * dx_alf;
    cl_msq = cl_msq + dpoint.x() * ag_msq;

    cpg1 = cpg2;
    cpg1_alf = cpg2_alf;
    cpg1_msq = cpg2_msq;
  }

  if (fabs(cl) > 0.0)
    aero.xcp /= cl;
  else
    aero.xcp = 0.0;

  return true;
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

bool XFoil::gamqv() {
  for (int i = 0; i < n; i++) {
    surface_vortex(0, i) = qvis[i];
    surface_vortex(1, i) = qinv_a[i];
  }

  return true;
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
    PsiResult psi_result =
        psilin(points, i, points.col(i),
               normal_vectors.col(i), true);

    const Vector2d res = qinf * Vector2d{points.col(i).y(),
                                         -points.col(i).x()};

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
    const double ag1 = atan2(-dpoints_ds.col(0).y(), -dpoints_ds.col(0).x());
    const double ag2 =
        atanc(dpoints_ds.col(n - 1).y(), dpoints_ds.col(n - 1).x(), ag1);
    const double abis = 0.5 * (ag1 + ag2);

    Vector2d bis_vector{cos(abis), sin(abis)};

    //----- minimum panel length adjacent to TE
    const double dsmin = std::min((points.col(0) - points.col(1)).norm(),
                                  (points.col(n - 1) - points.col(n - 2)).norm());

    //---- distance of internal control point ahead of sharp TE
    //-    (fraction of smaller panel length adjacent to TE)
    const double bwt = 0.1;

    //----- control point on bisector just ahead of TE point
    const Vector2d bis = point_te - bwt * dsmin * bis_vector;
    const Vector2d normal_bis{-bis_vector.y(), bis_vector.x()};

    //----- set velocity component along bisector line
    PsiResult psi_result = psilin(points, -1, bis, normal_bis, true);

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
bool XFoil::specal() {
  double minf_clm;
  double clm;

  //---- calculate surface vorticity distributions for alpha = 0, 90 degrees
  if (!lgamu || !lqaij)
    ggcalc();

  Matrix2d rotateMatrix =
      Matrix2d{{cos(alfa), sin(alfa)}, {-sin(alfa), cos(alfa)}};

  //---- superimpose suitably weighted  alpha = 0, 90  distributions
  for (int i = 0; i < n; i++) {
    surface_vortex(0, i) = rotateMatrix.row(0).dot(gamu.col(i));
    surface_vortex(1, i) = rotateMatrix.row(1).dot(gamu.col(i));
  }

  tecalc();
  qiset();

  //---- set initial guess for the newton variable clm
  clm = 1.0;

  //---- set corresponding  m(clm), re(clm)
  minf_clm = getActualMach(clm, mach_type);

  //---- set corresponding cl(m)
  clcalc(cmref);
  //---- iterate on clm
  bool bConv = false;
  for (int itcl = 1; itcl <= 20; itcl++) {
    const double msq_clm = 2.0 * minf * minf_clm;
    const double dclm = (cl - clm) / (1.0 - cl_msq * msq_clm);

    const double clm1 = clm;
    rlx = 1.0;

    //------ under-relaxation loop to avoid driving m(cl) above 1
    for (int irlx = 1; irlx <= 12; irlx++) {
      clm = clm1 + rlx * dclm;

      //-------- set new freestream mach m(clm)
      minf_clm = getActualMach(clm, mach_type);

      //-------- if mach is ok, go do next newton iteration
      // FIXME double型の==比較
      if (mach_type == MachType::CONSTANT || minf == 0.0 || minf_clm != 0.0)
        break;

      rlx = 0.5 * rlx;
    }

    //------ set new cl(m)
    clcalc(cmref);

    if (fabs(dclm) <= 1.0e-6) {
      bConv = true;
      break;
    }
  }
  if (!bConv) {
    writeString("Specal:  MInf convergence failed\n");
    return false;
  }

  //---- set final mach, cl, cp distributions, and hinge moment
  minf_cl = getActualMach(cl, mach_type);
  reinf_cl = getActualReynolds(cl, reynolds_type);
  comset();
  clcalc(cmref);

  cpi = cpcalc(n, qinv, qinf, minf);
  if (lvisc) {
    cpv = cpcalc(n + nw, qvis, qinf, minf);
    cpi = cpcalc(n + nw, qinv, qinf, minf);
  } else
    cpi = cpcalc(n, qinv, qinf, minf);

  for (int i = 0; i < n; i++) {
    qgamm[i + INDEX_START_WITH] = surface_vortex(0, i);
  }

  return true;
}


bool XFoil::speccl() {
  //-----------------------------------------
  //     converges to specified inviscid cl.
  //-----------------------------------------

  //---- calculate surface vorticity distributions for alpha = 0, 90 degrees
  if (!lgamu || !lqaij)
    ggcalc();

  //---- set freestream mach from specified cl -- mach will be held fixed
  minf_cl = getActualMach(clspec, mach_type);
  reinf_cl = getActualReynolds(clspec, reynolds_type);
  comset();

  Matrix2d rotateMatrix =
      Matrix2d{{cos(alfa), sin(alfa)}, {-sin(alfa), cos(alfa)}};

  //---- superimpose suitably weighted  alpha = 0, 90  distributions
  for (int i = 0; i < n; i++) {
    surface_vortex(0, i) = rotateMatrix.row(0).dot(gamu.col(i));
    surface_vortex(1, i) = rotateMatrix.row(1).dot(gamu.col(i));
  }

  //---- get corresponding cl, cl_alpha, cl_mach
  clcalc(cmref);

  //---- newton loop for alpha to get specified inviscid cl
  bool bConv = false;
  for (int ital = 1; ital <= 20; ital++) {
    const double dalfa = (clspec - cl) / cl_alf;
    rlx = 1.0;

    alfa = alfa + rlx * dalfa;

    //------ set new surface speed distribution
    Matrix2d rotateMatrix =
        Matrix2d{{cos(alfa), sin(alfa)}, {-sin(alfa), cos(alfa)}};
    for (int i = 0; i < n; i++) {
      surface_vortex(0, i) = rotateMatrix.row(0).dot(gamu.col(i));
      surface_vortex(1, i) = rotateMatrix.row(1).dot(gamu.col(i));
    }

    //------ set new cl(alpha)
    clcalc(cmref);

    if (fabs(dalfa) <= 1.0e-6) {
      bConv = true;
      break;
    }
  }
  if (!bConv) {
    writeString("Speccl:  cl convergence failed");
    return false;
  }

  //---- set final surface speed and cp distributions
  tecalc();
  qiset();

  if (lvisc) {
    cpv = cpcalc(n + nw, qvis, qinf, minf);
    cpi = cpcalc(n + nw, qinv, qinf, minf);

  } else {
    cpi = cpcalc(n, qinv, qinf, minf);
  }

  return true;
}
