/****************************************************************************

        XFoil Class
        Copyright (C) 2000 Mark Drela
        Copyright (C) 2003 Andre Deperrois techwinder@gmail.com

        This program is free software; you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation; either version 2 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program; if not, write to the Free Software
        Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

*****************************************************************************/

#include "XFoil.h"
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/StdVector"
#include "model/coefficient/skin_friction.hpp"
#include "model/coefficient/dissipation.hpp"
#include <algorithm>
#include <cstring>
#include <numbers>
using namespace Eigen;

// determinant
double cross2(const Eigen::Vector2d &a, const Eigen::Vector2d &b) {
  return a[0] * b[1] - a[1] * b[0];
}

bool XFoil::s_bCancel = false;
double XFoil::vaccel = 0.01;
const int INDEX_START_WITH = 1;

XFoil::XFoil() {
  m_pOutStream = nullptr;

  // fortran seems to initializes variables to 0
  mvisc = 0.0;

  // initialize transition parameters until user changes them
  acrit = 9.0;
  xstrip.top = 1.0;
  xstrip.bottom = 1.0;

  //---- initialize freestream mach number to zero
  mach_type = MachType::CONSTANT;
  minf1 = 0.0;

  //---- drop tolerance for bl system solver
  vaccel = 0.01;
  //---- default viscous parameters
  reynolds_type = ReynoldsType::CONSTANT;
  reinf1 = 0.0;
}

XFoil::~XFoil() {}

/** ---------------------------------------------------
 *      variable initialization/default routine.
 * --------------------------------------------------- */
// moved to XFoil_init.cpp: initialize()

/** -------------------------------------------------------
 * @brief Allocate and zero out large data containers.
 * ------------------------------------------------------- */
// moved to XFoil_init.cpp: initializeDataStructures()

/** -------------------------------------------------------
 * @brief Reset boolean state flags to defaults.
 * ------------------------------------------------------- */
// moved to XFoil_init.cpp: resetFlags()

/** -------------------------------------------------------
 * @brief Reset scalar state variables to default values.
 * ------------------------------------------------------- */
// moved to XFoil_init.cpp: resetVariables()

// moved to XFoil_boundary.cpp: computeShapeParameters()

/**
 * @brief Compute shear and skin friction coefficients.
 *
 * Determines the equilibrium shear coefficient (cq) and the
 * skin friction coefficient (cf) for the current flow regime.
 * Sensitivities with respect to the primary variables are also
 * accumulated in @p ref.
 *
 * @param ref           Boundary layer data container updated with the results.
 * @param flowRegimeType Flow regime type (1=laminar, 2=turbulent, 3=wake).
 */
// moved to XFoil_boundary.cpp: computeCoefficients()

/**
 * @brief Compute dissipation and boundary-layer thickness.
 *
 * Calculates the dissipation coefficient and boundary-layer thickness
 * based on the current flow regime.  This routine also accounts for
 * wake modifications and adds turbulent outer-layer contributions when
 * applicable.  Resulting derivatives are stored in @p ref.
 *
 * @param ref           Boundary layer data container updated with the results.
 * @param flowRegimeType Flow regime type (1=laminar, 2=turbulent, 3=wake).
 */
// moved to XFoil_boundary.cpp: computeDissipationAndThickness()

// moved to XFoil_geometry.cpp: abcopy()

// moved to XFoil_geometry.cpp: apcalc()

/** -------------------------------------------------------------
 *     atan2 function with branch cut checking.
 *
 *     increments position angle of point x,y from some previous
 *     value thold due to a change in position, ensuring that the
 *     position change does not cross the atan2 branch cut
 *     (which is in the -x direction).  for example:
 *
 *       atanc( -1.0 , -1.0 , 0.75*PI )  returns  1.25*PI , whereas
 *       atan2( -1.0 , -1.0 )            returns  -.75*PI .
 *
 *     typically, atanc is used to fill an array of angles:
 *
 *        theta(1) = atan2( y(1) , x(1) )
 *        do i=2, n
 *          theta[i] = atanc( y[i] , x[i] , theta(i-1) )
 *        end do
 *
 *     this will prevent the angle array theta(i) from jumping by
 *     +/- 2 pi when the path x(i),y(i) crosses the negative x axis.
 *
 *     input:
 *       x,y     point position coordinates
 *       thold   position angle of nearby point
 *
 *     output:
 *       atanc   position angle of x,y
 * -------------------------------------------------------------- */
// moved to XFoil_geometry.cpp: atanc()

/** ----------------------------------------------------------
 *      returns average amplification ax over interval 1..2
 * ----------------------------------------------------------- */
XFoil::AxResult XFoil::axset(double hk1, double t1, double rt1, double a1,
                             double hk2, double t2, double rt2, double a2,
                             double acrit) {
  AxResult result;
  //
  //==========================
  //---- 2nd-order
  double axsq = 0.0;
  double axa = 0.0, axa_ax1 = 0.0, axa_ax2 = 0.0;
  double exn = 0.0, exn_a1 = 0.0, exn_a2 = 0.0, dax = 0.0, dax_a1 = 0.0,
         dax_a2 = 0.0, dax_t1 = 0.0, dax_t2 = 0.0;
  double f_arg = 0.0; // ex arg

  EnvEnResult envEnResult1 = dampl(hk1, t1, rt1);
  EnvEnResult envEnResult2 = dampl(hk2, t2, rt2);

  //---- rms-average version (seems a little better on coarse grids)
  axsq = 0.5 * (envEnResult1.ax * envEnResult1.ax +
                envEnResult2.ax * envEnResult2.ax);
  if (axsq <= 0.0) {
    axa = 0.0;
    axa_ax1 = 0.0;
    axa_ax2 = 0.0;
  } else {
    axa = sqrt(axsq);
    axa_ax1 = 0.5 * envEnResult1.ax / axa;
    axa_ax2 = 0.5 * envEnResult2.ax / axa;
  }

  //----- small additional term to ensure  dn/dx > 0  near  n = ncrit
  f_arg = std::min(20.0 * (acrit - 0.5 * (a1 + a2)), 20.0);
  if (f_arg <= 0.0) {
    exn = 1.0;
    exn_a1 = 0.0;
    exn_a2 = 0.0;
  } else {
    exn = exp(-f_arg);
    exn_a1 = 20.0 * 0.5 * exn;
    exn_a2 = 20.0 * 0.5 * exn;
  }

  dax = exn * 0.002 / (t1 + t2);
  dax_a1 = exn_a1 * 0.002 / (t1 + t2);
  dax_a2 = exn_a2 * 0.002 / (t1 + t2);
  dax_t1 = -dax / (t1 + t2);
  dax_t2 = -dax / (t1 + t2);

  //==========================

  result.ax = axa + dax;
  result.ax_hk1 = axa_ax1 * envEnResult1.ax_hk;
  result.ax_t1 = axa_ax1 * envEnResult1.ax_th + dax_t1;
  result.ax_rt1 = axa_ax1 * envEnResult1.ax_rt;
  result.ax_a1 = dax_a1;

  result.ax_hk2 = axa_ax2 * envEnResult2.ax_hk;
  result.ax_t2 = axa_ax2 * envEnResult2.ax_th + dax_t2;
  result.ax_rt2 = axa_ax2 * envEnResult2.ax_rt;
  result.ax_a2 = dax_a2;

  return result;
}

/** -----------------------------------------------------------
 *     sets up the newton system coefficients and residuals
 *
 *         flowRegimeType = 0 :  similarity station
 *         flowRegimeType = 1 :  laminar interval
 *         flowRegimeType = 2 :  turbulent interval
 *         flowRegimeType = 3 :  wake interval
 *
 *      this routine knows nothing about a transition interval,
 *      which is taken care of by trdif.
 * ------------------------------------------------------------ */
// moved to XFoil_bldif.cpp: bldif()

/**
 * @brief Build laminar amplification equation coefficients.
 */
// moved to XFoil_bldif.cpp: bldifLaminar()

/**
 * @brief Build turbulent or wake shear lag equation coefficients.
 */
// moved to XFoil_bldif.cpp: bldifTurbulent()

/**
 * @brief Build momentum equation coefficients.
 */
// moved to XFoil_bldif.cpp: bldifMomentum()

/**
 * @brief Build shape parameter equation coefficients.
 */
// moved to XFoil_bldif.cpp: bldifShape()

bool XFoil::blkin() {
  //----------------------------------------------------------
  //     calculates turbulence-independent secondary "2"
  //     variables from the primary "2" variables.
  //----------------------------------------------------------
  double tr2, herat, he_u2, he_ms, v2_he;
  //---- set edge mach number ** 2
  blData2.param.mz =
      blData2.param.uz * blData2.param.uz * hstinv /
      (gm1bl * (1.0 - 0.5 * blData2.param.uz * blData2.param.uz * hstinv));
  tr2 = 1.0 + 0.5 * gm1bl * blData2.param.mz;
  blData2.param.mz_uz = 2.0 * blData2.param.mz * tr2 / blData2.param.uz;
  blData2.param.mz_ms =
      blData2.param.uz * blData2.param.uz * tr2 /
      (gm1bl * (1.0 - 0.5 * blData2.param.uz * blData2.param.uz * hstinv)) *
      hstinv_ms;

  //---- set edge density (isentropic relation)
  blData2.param.rz = rstbl * pow(tr2, (-1.0 / gm1bl));
  blData2.param.rz_uz = -blData2.param.rz / tr2 * 0.5 * blData2.param.mz_uz;
  blData2.param.rz_ms = -blData2.param.rz / tr2 * 0.5 * blData2.param.mz_ms +
                        rstbl_ms * pow(tr2, (-1.0 / gm1bl));

  //---- set shape parameter
  blData2.param.hz = blData2.param.dz / blData2.param.tz;
  blData2.param.hz_dz = 1.0 / blData2.param.tz;
  blData2.param.hz_tz = -blData2.param.hz / blData2.param.tz;

  //---- set edge static/stagnation enthalpy
  herat = 1.0 - 0.5 * blData2.param.uz * blData2.param.uz * hstinv;
  he_u2 = -blData2.param.uz * hstinv;
  he_ms = -0.5 * blData2.param.uz * blData2.param.uz * hstinv_ms;
  //---- set molecular viscosity
  v2_he = (1.5 / herat - 1.0 / (herat + hvrat));

  //---- set kinematic shape parameter
  boundary_layer::KineticShapeParameterResult hkin_result =
      boundary_layer::hkin(blData2.param.hz, blData2.param.mz);
  blData2.hkz.scalar = hkin_result.hk;

  blData2.hkz.u() = hkin_result.hk_msq * blData2.param.mz_uz;
  blData2.hkz.t() = hkin_result.hk_h * blData2.param.hz_tz;
  blData2.hkz.d() = hkin_result.hk_h * blData2.param.hz_dz;
  blData2.hkz.ms() = hkin_result.hk_msq * blData2.param.mz_ms;

  //---- set momentum thickness reynolds number
  blData2.rtz.scalar =
      blData2.param.rz * blData2.param.uz * blData2.param.tz /
      (sqrt(herat * herat * herat) * (1.0 + hvrat) / (herat + hvrat) / reybl);
  blData2.rtz.u() = blData2.rtz.scalar *
                    (1.0 / blData2.param.uz +
                     blData2.param.rz_uz / blData2.param.rz - v2_he * he_u2);
  blData2.rtz.t() = blData2.rtz.scalar / blData2.param.tz;
  blData2.rtz.ms() =
      blData2.rtz.scalar * (blData2.param.rz_ms / blData2.param.rz +
                            (1 / reybl * reybl_ms - v2_he * he_ms));
  blData2.rtz.re() = blData2.rtz.scalar * (reybl_re / reybl);

  return true;
}

bool XFoil::blmid(FlowRegimeEnum flowRegimeType) {
  //----------------------------------------------------
  //     calculates midpoint skin friction cfm
  //
  //      flowRegimeType = 1 :  laminar
  //      flowRegimeType = 2 :  turbulent
  //      flowRegimeType = 3 :  turbulent wake
  //----------------------------------------------------
  //

  //---- set similarity variables if not defined
  if (simi) {
    blData1.hkz = blData2.hkz;
    blData1.rtz = blData2.rtz;
    blData1.param.mz = blData2.param.mz;
    blData1.param.mz_uz = blData2.param.mz_uz;
    blData1.param.mz_ms = blData2.param.mz_ms;
  }

  //---- define stuff for midpoint cf
  double hka = 0.5 * (blData1.hkz.scalar + blData2.hkz.scalar);
  double rta = 0.5 * (blData1.rtz.scalar + blData2.rtz.scalar);
  double ma = 0.5 * (blData1.param.mz + blData2.param.mz);

  //---- compute midpoint skin friction coefficient
  skin_friction::C_f cf_res = skin_friction::getSkinFriction(hka, rta, ma, flowRegimeType);

  cfm = cf_res.cf;
  double cfm_hka = cf_res.hk;
  double cfm_rta = cf_res.rt;
  double cfm_ma = cf_res.msq;

  cfm_u1 = 0.5 * (cfm_hka * blData1.hkz.u() + cfm_ma * blData1.param.mz_uz +
                  cfm_rta * blData1.rtz.u());
  cfm_t1 = 0.5 * (cfm_hka * blData1.hkz.t() + cfm_rta * blData1.rtz.t());
  cfm_d1 = 0.5 * (cfm_hka * blData1.hkz.d());

  cfm_u2 = 0.5 * (cfm_hka * blData2.hkz.u() + cfm_ma * blData2.param.mz_uz +
                  cfm_rta * blData2.rtz.u());
  cfm_t2 = 0.5 * (cfm_hka * blData2.hkz.t() + cfm_rta * blData2.rtz.t());
  cfm_d2 = 0.5 * (cfm_hka * blData2.hkz.d());

  cfm_ms = 0.5 * (cfm_hka * blData1.hkz.ms() + cfm_ma * blData1.param.mz_ms +
                  cfm_rta * blData1.rtz.ms() + cfm_hka * blData2.hkz.ms() +
                  cfm_ma * blData2.param.mz_ms + cfm_rta * blData2.rtz.ms());
  cfm_re = 0.5 * (cfm_rta * blData1.rtz.re() + cfm_rta * blData2.rtz.re());

  return true;
}

/** ----------------------------------------------------------
 *     set bl primary "2" variables from parameter list
 *  ---------------------------------------------------------- */
bool XFoil::blprv(double xsi, double ami, double cti, double thi, double dsi,
                  double dswaki, double uei) {
  blData2.param.xz = xsi;
  blData2.param.amplz = ami;
  blData2.param.sz = cti;
  blData2.param.tz = thi;
  blData2.param.dz = dsi - dswaki;
  blData2.param.dwz = dswaki;

  blData2.param.uz =
      uei * (1.0 - tkbl) / (1.0 - tkbl * (uei / qinfbl) * (uei / qinfbl));
  blData2.param.uz_uei =
      (1.0 + tkbl * (2.0 * blData2.param.uz * uei / qinfbl / qinfbl - 1.0)) /
      (1.0 - tkbl * (uei / qinfbl) * (uei / qinfbl));
  blData2.param.uz_ms =
      (blData2.param.uz * (uei / qinfbl) * (uei / qinfbl) - uei) * tkbl_ms /
      (1.0 - tkbl * (uei / qinfbl) * (uei / qinfbl));
  return true;
}

/** -----------------------------------------------------------------
 *      custom solver for coupled viscous-inviscid newton system:
 *
 *        a  |  |  .  |  |  .  |    d       r       s
 *        b  a  |  .  |  |  .  |    d       r       s
 *        |  b  a  .  |  |  .  |    d       r       s
 *        .  .  .  .  |  |  .  |    d   =   r - dre s
 *        |  |  |  b  a  |  .  |    d       r       s
 *        |  z  |  |  b  a  .  |    d       r       s
 *        .  .  .  .  .  .  .  |    d       r       s
 *        |  |  |  |  |  |  b  a    d       r       s
 *
 *       a, b, z  3x3  blocks containing linearized bl equation coefficients
 *       |        3x1  vectors containing mass defect influence
 *                     coefficients on ue
 *       d        3x1  unknown vectors (newton deltas for ctau][ theta][ m)
 *       r        3x1  residual vectors
 *       s        3x1  re influence vectors
 * ------------------------------------------------------------------ */
// moved to XFoil_blsolve.cpp: plu/lu helpers and blsolve()
// moved to XFoil_blsolve.cpp: blsolve()

/** ----------------------------------------------------
 *      calculates all secondary "2" variables from
 *      the primary "2" variables x2, u2, t2, d2, s2.
 *      also calculates the sensitivities of the
 *      secondary variables wrt the primary variables.
 *
 *       flowRegimeType = 1 :  laminar
 *       flowRegimeType = 2 :  turbulent
 *       flowRegimeType = 3 :  turbulent wake
 * ---------------------------------------------------- */
bool XFoil::blvar(blData &ref, FlowRegimeEnum flowRegimeType) {
  // This routine is now decomposed into helper functions to simplify
  // the original Fortran translation.
  computeShapeParameters(ref, flowRegimeType);
  computeCoefficients(ref, flowRegimeType);
  computeDissipationAndThickness(ref, flowRegimeType);
  return true;
}

/** ------------------------------------------------------------------
 *
 *      sets up the bl newton system governing the current interval:
 *
 *      |       ||da1|     |       ||da2|       |     |
 *      |  vs1  ||dt1|  +  |  vs2  ||dt2|   =   |vsrez|
 *      |       ||dd1|     |       ||dd2|       |     |
 *               |du1|              |du2|
 *               |dx1|              |dx2|
 *
 *         3x5    5x1         3x5    5x1          3x1
 *
 *      the system as shown corresponds to a laminar station
 *      if tran, then  ds2  replaces  da2
 *      if turb, then  ds1, ds2  replace  da1, da2
 *
 * ------------------------------------------------------------------ */
bool XFoil::blsys() {

  //---- calculate secondary bl variables and their sensitivities
  if (wake) {
    blvar(blData2, FlowRegimeEnum::Wake);
    blmid(FlowRegimeEnum::Wake);
  } else {
    if (turb || tran) {
      blvar(blData2, FlowRegimeEnum::Turbulent);
      blmid(FlowRegimeEnum::Turbulent);
    } else {
      blvar(blData2, FlowRegimeEnum::Laminar);
      blmid(FlowRegimeEnum::Laminar);
    }
  }

  //---- for the similarity station, "1" and "2" variables are the same
  if (simi) {
    //		for(int icom=1;icom<= ncom;icom++) com1[icom] = com2[icom];
    stepbl();
  }

  //---- set up appropriate finite difference system for current interval
  if (tran)
    trdif();
  else if (simi)
    bldif(0);
  else if (!turb)
    bldif(1);
  else if (wake)
    bldif(3);
  else
    bldif(2);

  if (simi) {
    //----- at similarity station, "1" variables are really "2" variables
    blc.a2 += blc.a1;
    blc.a1 = Matrix<double, 4, 5>::Zero();
  }

  //---- change system over into incompressible uei and mach
  for (int k = 0; k < 4; k++) {
    //------ residual derivatives wrt compressible uec
    double res_u1 = blc.a1(k, 3);
    double res_u2 = blc.a2(k, 3);
    double res_ms = blc.d_msq[k];

    //------ combine with derivatives of compressible  u1,u2 = uec(uei m)
    blc.a1(k, 3) *= blData1.param.uz_uei;
    blc.a2(k, 3) *= blData2.param.uz_uei;
    blc.d_msq[k] =
        res_u1 * blData1.param.uz_ms + res_u2 * blData2.param.uz_ms + res_ms;
  }
  return true;
}

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

  xcp = 0.0;

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

    xcp += dpoint.x() * ag *
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
    xcp /= cl;
  else
    xcp = 0.0;

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

// moved to XFoil_init.cpp: writeString()

/** ==============================================================
 *      amplification rate routine for envelope e^n method.
 *      reference:
 *                 drela, m., giles, m.,
 *                "viscous/inviscid analysis of transonic and
 *                 low reynolds number airfoils",
 *                 aiaa journal, oct. 1987.
 *
 *      new version.   march 1991       (latest bug fix  july 93)
 *           - m(h) correlation made more accurate up to h=20
 *           - for h > 5, non-similar profiles are used
 *             instead of falkner-skan profiles.  these
 *             non-similar profiles have smaller reverse
 *             velocities, are more representative of typical
 *             separation bubble profiles.
 * --------------------------------------------------------------
 *
 *      input :   hk     kinematic shape parameter
 *                th     momentum thickness
 *                rt     momentum-thickness reynolds number
 *
 *      output:   ax     envelope spatial amplification rate
 *                ax_(.) sensitivity of ax to parameter (.)
 *
 *
 *      usage: the log of the envelope amplitude n(x) is
 *             calculated by integrating ax (= dn/dx) with
 *             respect to the streamwise distance x.
 *                       x
 *                      /
 *               n(x) = | ax(h(x),th(x),rth(x)) dx
 *                      /
 *                       0
 *             the integration can be started from the leading
 *             edge since ax will be returned as zero when rt
 *             is below the critical rtheta.  transition occurs
 *             when n(x) reaches ncrit (ncrit= 9 is "standard").
 * ============================================================== */
XFoil::EnvEnResult XFoil::dampl(double hk, double th, double rt) {
  EnvEnResult result;
  double dgr = 0.08;

  const double hmi = 1.0 / (hk - 1.0);
  const double hmi_hk = -hmi * hmi;

  //---- log10(critical rth) - h   correlation for falkner-skan profiles
  const double aa = 2.492 * pow(hmi, 0.43);
  const double aa_hk = (aa / hmi) * 0.43 * hmi_hk;
  const double bb = tanh(14.0 * hmi - 9.24);
  const double bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk;
  const double grcrit = aa + 0.7 * (bb + 1.0);
  const double grc_hk = aa_hk + 0.7 * bb_hk;
  const double gr = log10(rt);
  const double gr_rt = 1.0 / (2.3025851 * rt);
  if (gr < grcrit - dgr) {
    //----- no amplification for rtheta < rcrit
    result.ax = 0.0;
    result.ax_hk = 0.0;
    result.ax_th = 0.0;
    result.ax_rt = 0.0;
  } else {
    //----- set steep cubic ramp used to turn on ax smoothly as rtheta
    //-     exceeds rcrit (previously, this was done discontinuously).
    //-     the ramp goes between  -dgr < log10(rtheta/rcrit) < dgr

    const double rnorm = (gr - (grcrit - dgr)) / (2.0 * dgr);
    double rn_hk = -grc_hk / (2.0 * dgr);
    double rn_rt = gr_rt / (2.0 * dgr);

    double rfac, rfac_hk, rfac_rt;
    if (rnorm >= 1.0) {
      rfac = 1.0;
      rfac_hk = 0.0;
      rfac_rt = 0.0;
    } else {
      rfac = 3.0 * rnorm * rnorm - 2.0 * rnorm * rnorm * rnorm;
      const double rfac_rn = 6.0 * rnorm - 6.0 * rnorm * rnorm;

      rfac_hk = rfac_rn * rn_hk;
      rfac_rt = rfac_rn * rn_rt;
    }

    //----- amplification envelope slope correlation for falkner-skan
    const double f_arg = 3.87 * hmi - 2.52;
    const double arg_hk = 3.87 * hmi_hk;

    const double ex = exp(-f_arg * f_arg);
    const double ex_hk = ex * (-2.0 * f_arg * arg_hk);

    const double dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
    const double dadr_hk = 0.028 - 0.0345 * ex_hk;

    //----- new m(h) correlation    1 march 91
    const double af =
        -0.05 + 2.7 * hmi - 5.5 * hmi * hmi + 3.0 * hmi * hmi * hmi;
    const double af_hmi = 2.7 - 11.0 * hmi + 9.0 * hmi * hmi;
    const double af_hk = af_hmi * hmi_hk;

    result.ax = (af * dadr / th) * rfac;
    result.ax_hk = (af_hk * dadr / th + af * dadr_hk / th) * rfac +
                   (af * dadr / th) * rfac_hk;
    result.ax_th = -(result.ax) / th;
    result.ax_rt = (af * dadr / th) * rfac_rt;
  }

  return result;
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

/** -----------------------------------------------------------
 *     sets  bl location -> panel location  pointer array ipan
 * -----------------------------------------------------------*/
bool XFoil::iblpan() {
  std::stringstream ss;

  //-- top surface first
  // store ipan with 0-based BL station index, and set vti at 0-based
  for (int i = 0; i <= i_stagnation; i++) {
    ipan.top[i] = i_stagnation - i; // panel index
    vti.top[i] = 1.0;
  }
  // store TE as 0-based logical index
  iblte.top = i_stagnation;
  // nbl exclusive upper bound (0-based TE -> 1-based TE station + 1)
  nbl.top = iblte.top + 2;

  //-- bottom surface next
  // Bottom side: station 0 just after stagnation on bottom
  for (int index = 0; index <= n - i_stagnation; ++index) {
    ipan.bottom[index] = i_stagnation + INDEX_START_WITH + index;
    vti.bottom[index] = -1.0;
  }

  //-- wake
  iblte.bottom = n - (i_stagnation + INDEX_START_WITH) - 1; // logical 0-based TE

  for (int iw = 0; iw < nw; iw++) {
    int i = n + iw; // panel index in wake
    int index = iblte.bottom + iw + 2; // 1-based BL station for wake (bottom)
    ipan.bottom[index - 1] = i;        // ipan is 0-based in BL station
    vti.bottom[index - 1] = -1.0;
  }

  nbl.bottom = iblte.bottom + nw + 2;

  //-- upper wake pointers (for plotting only)
  for (int iw = 0; iw < nw; iw++) {
    // copy wake panel pointer from bottom to top (for plotting)
    ipan.top[iblte.top + iw + 1] = ipan.bottom[iblte.bottom + iw + 1]; // both sides are 0-based indices, hence -1 vs vti
    vti.top[iblte.top + iw + 1] = 1.0;
  }
  int iblmax = std::max(iblte.top, iblte.bottom) + nw + 2;
  if (iblmax > IVX) {
    ss << "iblpan :  ***  bl array overflow\n";
    ss << "Increase IVX to at least " << iblmax << "\n";
    writeString(ss.str());
    return false;
  }

  lipan = true;
  return true;
}

/** ---------------------------------------------
 *     sets the bl newton system line number
 *     corresponding to each bl station.
 * --------------------------------------------- */
bool XFoil::iblsys() {
  int iv = 0;
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      iv++;
      isys.get(is)[ibl] = iv;
    }
  }

  nsys = iv;
  if (nsys > 2 * IVX) {
    writeString("*** iblsys: bl system array overflow. ***");
    return false;
  }

  return true;
}

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

/** ------------------------------------------------------
 *     locates leading edge spline-parameter value sle
 *
 *     the defining condition is
 *
 *      (x-xte,y-yte) . (x',y') = 0     at  s = sle
 *
 *     i.e. the surface tangent is normal to the chord
 *     line connecting x(sle),y(sle) and the te point.
//------------------------------------------------------ */
// moved to XFoil_geometry.cpp: lefind()

/** ----------------------------------------------------
 *      marches the bls and wake in mixed mode using
 *      the current ue and hk.  the calculated ue
 *      and hk lie along a line quasi-normal to the
 *      natural ue-hk characteristic line of the
 *      current bl so that the goldstein or levy-lees
 *      singularity is never encountered.  continuous
 *      checking of transition onset is performed.
 * ----------------------------------------------------- */
bool XFoil::mrchdu() {
  std::stringstream ss;

  const double deps = 0.000005;
  int itrold, iw = 0; // icom

  double senswt = 0.0, thi, uei, dsi, cti, dswaki, ratlen = 0.0;
  double sens = 0.0, sennew = 0.0, msq = 0.0, thm = 0.0, dsm = 0.0, uem = 0.0;
  double xsi, hklim = 0.0, dsw = 0.0;
  double ami = 0.0, dte = 0.0, cte = 0.0, tte = 0.0, ueref = 0.0, hkref = 0.0,
         dmax = 0.0;

  //---- constant controlling how far hk is allowed to deviate
  //    from the specified value.
  senswt = 1000.0;
  sens = 0.0;
  sennew = 0.0;

  for (int is = 1; is <= 2; is++) {
    //---- set forced transition arc length position
    xiforc = xifset(is);

    //---- old transition station
    itrold = itran.get(is);

    tran = false;
    turb = false;
    itran.get(is) = iblte.get(is);
    //---- march downstream
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      //int ibl = ibl - 1;
      //int ibm = ibl - 1;

      simi = (ibl == 0);
      wake = ibl > iblte.get(is);

      //------ initialize current station to existing variables (xssi now 0-based)
      xsi = xssi.get(is)[ibl];
      uei = uedg.get(is)[ibl];
      thi = thet.get(is)[ibl];
      dsi = dstr.get(is)[ibl];

      //------ fixed bug   md 7 june 99
      if (ibl < itrold) {
        ami = ctau.get(is)[ibl]; // ami must be initialized
        cti = 0.03;
      } else {
        cti = ctau.get(is)[ibl];
        if (cti <= 0.0)
          cti = 0.03;
      }

      if (wake) {
        iw = ibl - iblte.get(is);
        dswaki = wgap[iw - 1];
      } else
        dswaki = 0.0;

      if (ibl <= iblte.get(is))
        dsi =std::max(dsi - dswaki, 1.02000 * thi) + dswaki;
      if (ibl > iblte.get(is))
        dsi = std::max(dsi - dswaki, 1.00005 * thi) + dswaki;

      //------ newton iteration loop for current station

      bool converged = false;
      for (int itbl = 1; itbl <= 25; itbl++) {
        //-------- assemble 10x3 linearized system for dctau, dth, dds, due, dxi
        //         at the previous "1" station and the current "2" station
        //         (the "1" station coefficients will be ignored)

        blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
        blkin();

        //-------- check for transition and set appropriate flags and things
        if ((!simi) && (!turb)) {
          trchek();
          ami = blData2.param.amplz;
          if (tran)
            itran.get(is) = ibl;
          if (!tran)
            itran.get(is) = ibl + 1;
        }
        if (ibl == iblte.get(is) + 1) {
          tte = thet.get(1)[iblte.top] +
                thet.get(2)[iblte.bottom];
          dte = dstr.get(1)[iblte.top] +
                dstr.get(2)[iblte.bottom] + ante;
          cte = (ctau.get(1)[iblte.top] *
                     thet.get(1)[iblte.top] +
                 ctau.get(2)[iblte.bottom] *
                     thet.get(2)[iblte.bottom]) /
                tte;
          tesys(cte, tte, dte);
        } else {
          blsys();
        }

        //-------- set stuff at first iteration...
        if (itbl == 1) {
          //--------- set "baseline" ue and hk for forming  ue(hk)  relation
          ueref = blData2.param.uz;
          hkref = blData2.hkz.scalar;

          //--------- if current point ibl was turbulent and is now laminar,
          // then...
          if (ibl < itran.get(is) && ibl >= itrold) {
            //---------- extrapolate baseline hk
            if (ibl > 0) {
              uem = uedg.get(is)[ibl - 1];
              dsm = dstr.get(is)[ibl - 1];
              thm = thet.get(is)[ibl - 1];
            } else {
              uem = uedg.get(is)[ibl];
              dsm = dstr.get(is)[ibl];
              thm = thet.get(is)[ibl];
            }
            msq =
                uem * uem * hstinv / (gm1bl * (1.0 - 0.5 * uem * uem * hstinv));
            boundary_layer::KineticShapeParameterResult hkin_result =
                boundary_layer::hkin(dsm / thm, msq);
            hkref = hkin_result.hk;
          }

          //--------- if current point ibl was laminar, then...
          if (ibl < itrold) {
            //---------- reinitialize or extrapolate ctau if it's now turbulent
            if (tran)
              ctau.get(is)[ibl] = 0.03;
            if (turb) {
              double prev = (ibl >= 1) ? ctau.get(is)[ibl - 1]
                                     : ctau.get(is)[ibl];
              ctau.get(is)[ibl] = prev;
            }
            if (tran || turb) {
              cti = ctau.get(is)[ibl - 1];
              blData2.param.sz = cti;
            }
          }
        }

        if (simi || ibl == iblte.get(is) + 1) {
          //--------- for similarity station or first wake point, prescribe ue
          blc.a2(3, 0) = 0.0;
          blc.a2(3, 1) = 0.0;
          blc.a2(3, 2) = 0.0;
          blc.a2(3, 3) = blData2.param.uz_uei;
          blc.rhs[3] = ueref - blData2.param.uz;
        } else {
          //******** calculate ue-hk characteristic slope

          //--------- set unit dhk
          blc.a2(3, 0) = 0.0;
          blc.a2(3, 1) = blData2.hkz.t();
          blc.a2(3, 2) = blData2.hkz.d();
          blc.a2(3, 3) = blData2.hkz.u() * blData2.param.uz_uei;
          blc.rhs[3] = 1.0;

          //--------- calculate due response
          double delta_sen = blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs)[3];

          //--------- set  senswt * (normalized due/dhk)
          sennew = senswt * delta_sen * hkref / ueref;
          if (itbl <= 5)
            sens = sennew;
          else if (itbl <= 15)
            sens = 0.5 * (sens + sennew);

          //--------- set prescribed ue-hk combination
          blc.a2(3, 0) = 0.0;
          blc.a2(3, 1) = blData2.hkz.t() * hkref;
          blc.a2(3, 2) = blData2.hkz.d() * hkref;
          blc.a2(3, 3) =
              (blData2.hkz.u() * hkref + sens / ueref) * blData2.param.uz_uei;
          blc.rhs[3] = -(hkref * hkref) * (blData2.hkz.scalar / hkref - 1.0) -
                     sens * (blData2.param.uz / ueref - 1.0);
        }

        //-------- solve newton system for current "2" station
        blc.rhs = blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);

        //-------- determine max changes and underrelax if necessary
        dmax = std::max(fabs(blc.rhs[1] / thi), fabs(blc.rhs[2] / dsi));
        if (ibl >= itran.get(is))
          dmax = std::max(dmax, fabs(blc.rhs[0] / (10.0 * cti)));

        rlx = 1.0;
        if (dmax > 0.3)
          rlx = 0.3 / dmax;

        //-------- update as usual
        if (ibl < itran.get(is))
          ami = ami + rlx * blc.rhs[0];
        if (ibl >= itran.get(is))
          cti = cti + rlx * blc.rhs[0];
        thi = thi + rlx * blc.rhs[1];
        dsi = dsi + rlx * blc.rhs[2];
        uei = uei + rlx * blc.rhs[3];

        //-------- eliminate absurd transients
        if (ibl >= itran.get(is)) {
          cti = std::min(cti, 0.30);
          cti = std::max(cti, 0.0000001);
        }

        if (ibl <= iblte.get(is))
          hklim = 1.02;
        else
          hklim = 1.00005;

        msq = uei * uei * hstinv / (gm1bl * (1.0 - 0.5 * uei * uei * hstinv));
        dsw = dsi - dswaki;
        dslim(dsw, thi, msq, hklim);
        dsi = dsw + dswaki;

        if (dmax <= deps) {
          converged = true;
          break;
        }
      }

      if (!converged) {
        ss << "     mrchdu: convergence failed at " << ibl << " ,  side " << is
           << ", res=" << std::setw(4) << std::fixed << std::setprecision(3)
           << dmax << "\n";
        writeString(ss.str());
        ss.str("");

        //------ the current unconverged solution might still be reasonable...
        if (dmax > 0.1) {
          //------- the current solution is garbage --> extrapolate values
          if (ibl >= 2) {
            if (ibl <= iblte.get(is)) {
              thi = thet.get(is)[ibl - 1] *
                    sqrt(xssi.get(is)[ibl] / xssi.get(is)[ibl - 1]);
              dsi = dstr.get(is)[ibl - 1] *
                    sqrt(xssi.get(is)[ibl] / xssi.get(is)[ibl - 1]);
              uei = uedg.get(is)[ibl - 1];
            } else {
              if (ibl == iblte.get(is) + 1) {
                cti = cte;
                thi = tte;
                dsi = dte;
                uei = uedg.get(is)[ibl - 1];
              } else {
                thi = thet.get(is)[ibl - 1];
                ratlen = (xssi.get(is)[ibl] - xssi.get(is)[ibl - 1]) /
                         (10.0 * dstr.get(is)[ibl - 1]);
                dsi = (dstr.get(is)[ibl - 1] +
                       thi * ratlen) /
                      (1.0 + ratlen);
                uei = uedg.get(is)[ibl - 1];
              }
            }
            if (ibl == itran.get(is))
              cti = 0.05;
            if (ibl > itran.get(is))
              cti = ctau.get(is)[ibl - 1];
          }
        }

        blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
        blkin();

        //------- check for transition and set appropriate flags and things
        if ((!simi) && (!turb)) {
          trchek();
          ami = blData2.param.amplz;
          if (tran)
            itran.get(is) = ibl;
          if (!tran)
            itran.get(is) = ibl + 2;
        }

        //------- set all other extrapolated values for current station
        if (ibl < itran.get(is))
          blvar(blData2, FlowRegimeEnum::Laminar);
        if (ibl >= itran.get(is))
          blvar(blData2, FlowRegimeEnum::Turbulent);
        if (wake)
          blvar(blData2, FlowRegimeEnum::Wake);

        if (ibl < itran.get(is))
          blmid(FlowRegimeEnum::Laminar);
        if (ibl >= itran.get(is))
          blmid(FlowRegimeEnum::Turbulent);
        if (wake)
          blmid(FlowRegimeEnum::Wake);
      }

      //------ pick up here after the newton iterations
      sens = sennew;

      //------ store primary variables
      if (ibl < itran.get(is))
        ctau.get(is)[ibl] = ami;
      else
        ctau.get(is)[ibl] = cti;
      thet.get(is)[ibl] = thi;
      dstr.get(is)[ibl] = dsi;
      uedg.get(is)[ibl] = uei;
      mass.get(is)[ibl] = dsi * uei;
      ctq.get(is)[ibl] = blData2.cqz.scalar;

      //------ set "1" variables to "2" variables for next streamwise station
      blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
      blkin();

      stepbl();

      //------ turbulent intervals will follow transition interval or te
      if (tran || ibl == iblte.get(is)) {
        turb = true;
      }

      tran = false;
      //                        qApp->processEvents();
      if (s_bCancel)
        return false;
    }
  }
  return true;
}
double XFoil::calcHtarg(int ibl, int is, bool wake) {
  if (ibl < itran.get(is)) {
    return blData1.hkz.scalar +
           0.03 * (blData2.param.xz - blData1.param.xz) / blData1.param.tz;
  } else if (ibl == itran.get(is)) {
    return blData1.hkz.scalar +
           (0.03 * (xt - blData1.param.xz) - 0.15 * (blData2.param.xz - xt)) /
               blData1.param.tz;
  } else if (wake) {
    const double cst =
        0.03 * (blData2.param.xz - blData1.param.xz) / blData1.param.tz;
    auto euler = [](double hk2, double hk1, double cst) {
      return hk2 - (hk2 + cst * pow(hk2 - 1, 3) - hk1) /
                       (1 + 3 * cst * pow(hk2 - 1, 2));
    };
    blData2.hkz.scalar = blData1.hkz.scalar;
    for (int i = 0; i < 3; i++) {
      blData2.hkz.scalar = euler(blData2.hkz.scalar, blData1.hkz.scalar, cst);
    }
    return blData2.hkz.scalar;
  } else {
    return blData1.hkz.scalar -
           0.15 * (blData2.param.xz - blData1.param.xz) / blData1.param.tz;
  }
}

/** ----------------------------------------------------
 *      marches the bls and wake in direct mode using
 *      the uedg array. if separation is encountered,
 *      a plausible value of hk extrapolated from
 *      upstream is prescribed instead.  continuous
 *      checking of transition onset is performed.
 * ----------------------------------------------------*/
bool XFoil::mrchue() {
  std::stringstream ss;
  bool direct;

  double msq, ratlen, dsw, hklim;
  double xsi, uei, ucon, tsq, thi, ami, cti, dsi;
  double dswaki;
  double htest, hktest;
  double cte, dte, tte, dmax, hmax, htarg = 0.0;
  // double cte = dte = tte = dmax = hmax = htarg = 0.0;

  //---- shape parameters for separation criteria
  const double hlmax = 3.8;
  const double htmax = 2.5;

  for (int is = 1; is <= 2; is++) { // 2000
    ss << "    Side " << is << " ...\n";
    writeString(ss.str());

    //---- set forced transition arc length position
    xiforc = xifset(is);

    //---- initialize similarity station with thwaites' formula (keep ibl=1)
    // xssi is 0-based: first BL station arc length at index 0
    xsi = xssi.get(is)[0];
    uei = uedg.get(is)[0];

    ucon = uei / xsi;
    tsq = 0.45 / (ucon * 6.0 * reybl);
    thi = sqrt(tsq);
    dsi = 2.2 * thi;
    ami = 0.0;

    //---- initialize ctau for first turbulent station
    cti = 0.03;

    tran = false;
    turb = false;
    itran.get(is) = iblte.get(is);

    //---- march downstream
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      //int ibl = ibl - 1;
      int ibm = ibl + 1;
      int iw = ibl - iblte.get(is);
      simi = (ibl == 0);
      wake = ibl > iblte.get(is);

      //------ prescribed quantities (xssi now 0-based)
      xsi = xssi.get(is)[ibl];
      uei = uedg.get(is)[ibl];

      if (wake) {
        iw = ibl - iblte.get(is);
        dswaki = wgap[iw - 1];
      } else
        dswaki = 0.0;

      direct = true;
      bool converged = false;

      //------ newton iteration loop for current station
      for (int itbl = 1; itbl <= 25; itbl++) { // 100

        //-------- assemble 10x3 linearized system for dctau, dth, dds, due, dxi
        //         at the previous "1" station and the current "2" station
        //         (the "1" station coefficients will be ignored)

        blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
        blkin();

        //-------- check for transition and set appropriate flags and things
        if ((!simi) && (!turb)) {
          trchek();
          ami = blData2.param.amplz;

          //--------- fixed bug   md 7 jun 99
          if (tran) {
            itran.get(is) = ibl;
            if (cti <= 0.0) {
              cti = 0.03;
              blData2.param.sz = cti;
            }
          } else
            itran.get(is) = ibl + 2;
        }

        if (ibl == iblte.get(is) + 1) {
          tte = thet.get(1)[iblte.top] +
                thet.get(2)[iblte.bottom];
          dte = dstr.get(1)[iblte.top] +
                dstr.get(2)[iblte.bottom] + ante;
          cte = (ctau.get(1)[iblte.top] *
                     thet.get(1)[iblte.top] +
                 ctau.get(2)[iblte.bottom] *
                     thet.get(2)[iblte.bottom]) /
                tte;
          tesys(cte, tte, dte);
        } else
          blsys();

        if (direct) {
          //--------- try direct mode (set due = 0 in currently empty 4th line)
          blc.a2(3, 0) = 0.0;
          blc.a2(3, 1) = 0.0;
          blc.a2(3, 2) = 0.0;
          blc.a2(3, 3) = 1.0;
          blc.rhs[3] = 0.0;
          //--------- solve newton system for current "2" station
          blc.rhs = blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);
          //--------- determine max changes and underrelax if necessary
          dmax = std::max(fabs(blc.rhs[1] / thi), fabs(blc.rhs[2] / dsi));
          if (ibl < itran.get(is))
            dmax = std::max(dmax, fabs(blc.rhs[0] / 10.0));
          if (ibl >= itran.get(is))
            dmax = std::max(dmax, fabs(blc.rhs[0] / cti));

          rlx = 1.0;
          if (dmax > 0.3)
            rlx = 0.3 / dmax;
          //--------- see if direct mode is not applicable
          if (ibl != iblte.get(is) + 1) {
            //---------- calculate resulting kinematic shape parameter hk
            msq =
                uei * uei * hstinv / (gm1bl * (1.0 - 0.5 * uei * uei * hstinv));
            htest = (dsi + rlx * blc.rhs[2]) / (thi + rlx * blc.rhs[1]);
            boundary_layer::KineticShapeParameterResult hkin_result =
                boundary_layer::hkin(htest, msq);
            hktest = hkin_result.hk;

            //---------- decide whether to do direct or inverse problem based on
            // hk
            if (ibl < itran.get(is))
              hmax = hlmax;
            if (ibl >= itran.get(is))
              hmax = htmax;
            direct = (hktest < hmax);
          }
          if (direct) {
            //---------- update as usual
            if (ibl >= itran.get(is))
              cti = cti + rlx * blc.rhs[0];
            thi = thi + rlx * blc.rhs[1];
            dsi = dsi + rlx * blc.rhs[2];
          } else {
            //---------- set prescribed hk for inverse calculation at the
            // current station
            htarg = calcHtarg(ibl, is, wake);

            //---------- limit specified hk to something reasonable
            if (wake)
              htarg = std::max(htarg, 1.01);
            else
              htarg = std::max(htarg, hmax);

            ss.str("");
            ss << "     mrchue: inverse mode at " << ibl
               << "    hk=" << std::fixed << setprecision(3) << htarg << "\n";
            writeString(ss.str());
            ss.str("");

            //---------- try again with prescribed hk

            continue;
          }
        } else {
          //-------- inverse mode (force hk to prescribed value htarg)
          blc.a2(3, 0) = 0.0;
          blc.a2(3, 1) = blData2.hkz.t();
          blc.a2(3, 2) = blData2.hkz.d();
          blc.a2(3, 3) = blData2.hkz.u();
          blc.rhs[3] = htarg - blData2.hkz.scalar;
          blc.rhs = blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);

          dmax = std::max(fabs(blc.rhs[1] / thi), fabs(blc.rhs[2] / dsi));
          if (ibl >= itran.get(is))
            dmax = std::max(dmax, fabs(blc.rhs[0] / cti));
          rlx = 1.0;
          if (dmax > 0.3)
            rlx = 0.3 / dmax;
          //--------- update variables
          if (ibl >= itran.get(is))
            cti = cti + rlx * blc.rhs[0];
          thi = thi + rlx * blc.rhs[1];
          dsi = dsi + rlx * blc.rhs[2];
          uei = uei + rlx * blc.rhs[3];
        }
        //-------- eliminate absurd transients

        if (ibl >= itran.get(is)) {
          cti = std::min(cti, 0.30);
          cti = std::max(cti, 0.0000001);
        }
        if (ibl <= iblte.get(is))
          hklim = 1.02;
        else
          hklim = 1.00005;
        msq = uei * uei * hstinv / (gm1bl * (1.0 - 0.5 * uei * uei * hstinv));
        dsw = dsi - dswaki;
        dslim(dsw, thi, msq, hklim);
        dsi = dsw + dswaki;
        if (dmax <= 0.00001) {
          converged = true;
          break;
        }

      } // end itbl loop

      if (!converged) {
        ss.str("");
        ss << "     mrchue: convergence failed at " << ibl << ",  side " << is
           << ", res =" << std::fixed << std::setprecision(3) << dmax << "\n";
        writeString(ss.str());

        //------ the current unconverged solution might still be reasonable...
        if (dmax > 0.1) {
          //------- the current solution is garbage --> extrapolate values
          // instead
          if (ibl >= 2) {
            if (ibl <= iblte.get(is)) {
              thi = thet.get(is)[ibl - 1] *
                    sqrt(xssi.get(is)[ibl] / xssi.get(is)[ibl - 1]);
              dsi = dstr.get(is)[ibl - 1] *
                    sqrt(xssi.get(is)[ibl] / xssi.get(is)[ibl - 1]);
            } else {
              if (ibl == iblte.get(is) + 1) {
                cti = cte;
                thi = tte;
                dsi = dte;
              } else {
                thi = thet.get(is)[ibl - 1];
                ratlen = (xssi.get(is)[ibl] - xssi.get(is)[ibl - 1]) /
                         (10.0 * dstr.get(is)[ibl - 1]);
                dsi = (dstr.get(is)[ibl - 1] +
                       thi * ratlen) /
                      (1.0 + ratlen);
              }
            }
            if (ibl == itran.get(is))
              cti = 0.05;
            if (ibl > itran.get(is))
              cti = ctau.get(is)[ibl - 1];

            uei = uedg.get(is)[ibl];

            if (ibl < nbl.get(is) - 1)
              uei = 0.5 *
                    (uedg.get(is)[ibl - 1] +
                     uedg.get(is)[ibl + 1]);
          }
        }
        // 109
        blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
        blkin();
        //------- check for transition and set appropriate flags and things
        if ((!simi) && (!turb)) {
          trchek();
          ami = blData2.param.amplz;
          if (tran)
            itran.get(is) = ibl;
          if (!tran)
            itran.get(is) = ibl + 2;
        }
        //------- set all other extrapolated values for current station
        if (ibl < itran.get(is))
          blvar(blData2, FlowRegimeEnum::Laminar);
        if (ibl >= itran.get(is))
          blvar(blData2, FlowRegimeEnum::Turbulent);
        if (wake)
          blvar(blData2, FlowRegimeEnum::Wake);
        if (ibl < itran.get(is))
          blmid(FlowRegimeEnum::Laminar);
        if (ibl >= itran.get(is))
          blmid(FlowRegimeEnum::Turbulent);
        if (wake)
          blmid(FlowRegimeEnum::Wake);
      }
      //------ store primary variables
      if (ibl < itran.get(is))
        ctau.get(is)[ibl] = ami;
      if (ibl >= itran.get(is))
        ctau.get(is)[ibl] = cti;
      thet.get(is)[ibl] = thi;
      dstr.get(is)[ibl] = dsi;
      uedg.get(is)[ibl] = uei;
      mass.get(is)[ibl] = dsi * uei;
      ctq.get(is)[ibl] = blData2.cqz.scalar;

      //------ set "1" variables to "2" variables for next streamwise station
      blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
      blkin();

      stepbl();

      //------ turbulent intervals will follow transition interval or te
      if (tran || ibl == iblte.get(is)) {
        turb = true;
      }

      tran = false;

      if (ibl == iblte.get(is)) {
        thi = thet.get(1)[iblte.top] +
              thet.get(2)[iblte.bottom];
        dsi = dstr.get(1)[iblte.top] +
              dstr.get(2)[iblte.bottom] + ante;
      }
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
    for (int jbl = 1; jbl < nbl.get(js); ++jbl) {
      int j = ipan.get(js)[jbl - 1];
      int jv = isys.get(js)[jbl - 1];
      ule1_m[jv] = -vti.get(1)[0] * vti.get(js)[jbl - 1] *
                   dij(ile1, j);
      ule2_m[jv] = -vti.get(2)[0] * vti.get(js)[jbl - 1] *
                   dij(ile2, j);
      ute1_m[jv] = -vti.get(1)[iblte.top] * vti.get(js)[jbl - 1] *
                   dij(ite1, j);
      ute2_m[jv] = -vti.get(2)[iblte.bottom] * vti.get(js)[jbl - 1] *
                   dij(ite2, j);
    }
  }
}

void XFoil::clearDerivativeVectors(VectorXd &u_m, VectorXd &d_m) {
  for (int js = 1; js <= 2; ++js) {
    for (int jbl = 1; jbl < nbl.get(js); ++jbl) {
      int jv = isys.get(js)[jbl - 1];
      u_m[jv] = 0.0;
      d_m[jv] = 0.0;
    }
  }
}

bool XFoil::setbl() {
  //-------------------------------------------------
  //	   sets up the bl newton system coefficients for the current bl
  // variables
  //     and the edge velocities received from setup. the local bl system
  //     coefficients are then incorporated into the global newton system.
  //-------------------------------------------------

  std::stringstream ss;
  int jvte1 = 0, jvte2 = 0;
  VectorXd u1_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd u2_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd d1_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd d2_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd ule1_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd ule2_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd ute1_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd ute2_m = VectorXd::Zero(2 * IVX + 1);

  double msq_clmr = 0.0, mdi;
  double herat = 0.0, herat_ms = 0.0;

  double clmr = 0.0, ma_clmr = 0.0, re_clmr = 0.0;
  double ule1_a = 0.0, ule2_a = 0.0, u2_a, due2, dds2;
  double xsi, cti = 0.0, uei, thi, dsi, dswaki;
  double d2_a, d2_m2, d2_u2, dte_mte1, dte_ute1, dte_mte2, dte_ute2;
  double tte, cte, dte, dule1 = 0.0, dule2 = 0.0;
  double xi_ule1, xi_ule2;
  double ami = 0.0, tte_tte1 = 0.0, tte_tte2 = 0.0, cte_tte1 = 0.0,
         cte_tte2 = 0.0, cte_cte1 = 0.0, cte_cte2 = 0.0;

  //---- set the cl used to define mach, reynolds numbers
  if (lalfa)
    clmr = cl;
  else
    clmr = clspec;

  cti = 0.0; // techwinder added, otherwise variable is not initialized

  //---- set current minf(cl)
  ma_clmr = getActualMach(clmr, mach_type);
  re_clmr = getActualReynolds(clmr, reynolds_type);

  msq_clmr = 2.0 * minf * ma_clmr;

  //---- set compressibility parameter tklam and derivative tk_msq
  comset();

  //---- set gas constant (= cp/cv)
  gm1bl = gamm1;

  //---- set parameters for compressibility correction
  qinfbl = qinf;
  tkbl = tklam;
  tkbl_ms = tkl_msq;

  //---- stagnation density and 1/enthalpy
  rstbl = pow((1.0 + 0.5 * gm1bl * minf * minf), (1.0 / gm1bl));
  rstbl_ms = 0.5 * rstbl / (1.0 + 0.5 * gm1bl * minf * minf);
  hstinv = gm1bl * (minf / qinfbl) * (minf / qinfbl) /
           (1.0 + 0.5 * gm1bl * minf * minf);
  hstinv_ms = gm1bl * (1.0 / qinfbl) * (1.0 / qinfbl) /
                  (1.0 + 0.5 * gm1bl * minf * minf) -
              0.5 * gm1bl * hstinv / (1.0 + 0.5 * gm1bl * minf * minf);

  //---- set reynolds number based on freestream density, velocity, viscosity
  herat = 1.0 - 0.5 * qinfbl * qinfbl * hstinv;
  herat_ms = -0.5 * qinfbl * qinfbl * hstinv_ms;

  reybl = reinf * sqrt(herat * herat * herat) * (1.0 + hvrat) / (herat + hvrat);
  reybl_re = sqrt(herat * herat * herat) * (1.0 + hvrat) / (herat + hvrat);
  reybl_ms = reybl * (1.5 / herat - 1.0 / (herat + hvrat)) * herat_ms;

  amcrit = acrit;

  if (!lblini) {
    //----- initialize bl by marching with ue (fudge at separation)
    // TRACE(" initializing bl ...\n");
    writeString("   Initializing bl ...\n");

    mrchue();
    lblini = true;
  }

  //---- march bl with current ue and ds to establish transition
  mrchdu();

  SidePair<VectorXd> usav = uedg;

  ueset();
  swapEdgeVelocities(usav);
jvte1 = isys.top[iblte.top];
jvte2 = isys.bottom[iblte.bottom];

  dule1 = uedg.top[0] - usav.top[0];
  dule2 = uedg.bottom[0] - usav.bottom[0];

  //---- set le and te ue sensitivities wrt all m values
  computeLeTeSensitivities(
    ipan.get(1)[0],
    ipan.get(2)[0],
    ipan.get(1)[iblte.top],
    ipan.get(2)[iblte.bottom],
    ule1_m, ule2_m, ute1_m, ute2_m
  );

  ule1_a = uinv_a.get(1)[0];
  ule2_a = uinv_a.get(2)[0];

  writeString(" \n");

  //*** go over each boundary layer/wake
  for (int is = 1; is <= 2; is++) {
    //---- there is no station "1" at similarity, so zero everything out
    clearDerivativeVectors(u1_m, d1_m);
    double u1_a = 0.0;
    double d1_a = 0.0;

    double due1 = 0.0;
    double dds1 = 0.0;

    //---- set forced transition arc length position
    xiforc = xifset(is);

    tran = false;
    turb = false;

    //**** sweep downstream setting up bl equation linearizations
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      
      int iv = isys.get(is)[ibl];

      simi = (ibl == 0);
      wake = (ibl > iblte.get(is));
      tran = (ibl == itran.get(is));
      turb = (ibl > itran.get(is));

      //---- set primary variables for current station
      xsi = xssi.get(is)[ibl];
      if (ibl < itran.get(is))
        ami = ctau.get(is)[ibl];
      else
        cti = ctau.get(is)[ibl];
      uei = uedg.get(is)[ibl];
      thi = thet.get(is)[ibl];
      mdi = mass.get(is)[ibl];

      dsi = mdi / uei;

      if (wake) {
        int iw = ibl - iblte.get(is);
        dswaki = wgap[iw - 1];
      } else
        dswaki = 0.0;

      //---- set derivatives of dsi (= d2)
      d2_m2 = 1.0 / uei;
      d2_u2 = -dsi / uei;

      for (int js = 1; js <= 2; js++) {
        for (int jbl = 0; jbl < nbl.get(js) - 1; ++jbl) {
          int jv = isys.get(js)[jbl];
          u2_m[jv] = -vti.get(is)[ibl] * vti.get(js)[jbl] *
                     dij(ipan.get(is)[ibl], ipan.get(js)[jbl]);
          d2_m[jv] = d2_u2 * u2_m[jv];
        }
      }
      d2_m[iv] = d2_m[iv] + d2_m2;

      u2_a = uinv_a.get(is)[ibl];
      d2_a = d2_u2 * u2_a;

      //---- "forced" changes due to mismatch between uedg and
      // usav=uinv+dij*mass
      due2 = uedg.get(is)[ibl] - usav.get(is)[ibl];
      dds2 = d2_u2 * due2;

      blprv(xsi, ami, cti, thi, dsi, dswaki, uei); // cti
      blkin();

      //---- check for transition and set tran, xt, etc. if found
      if (tran) {
        trchek();
        ami = blData2.param.amplz;
      }

      if (ibl == itran.get(is) && !tran) {
        // TRACE("setbl: xtr???  n1=%d n2=%d: \n", ampl1, ampl2);

        ss << "setbl: xtr???  n1=" << blData1.param.amplz
           << " n2=" << blData2.param.amplz << ":\n";
        writeString(ss.str());
        ss.str("");
      }

      //---- assemble 10x4 linearized system for dctau, dth, dds, due, dxi
      //	   at the previous "1" station and the current "2" station

      if (ibl == iblte.get(is) + 1) {
        //----- define quantities at start of wake, adding te base thickness to
        // dstar
        tte = thet.get(1)[iblte.top] +
              thet.get(2)[iblte.bottom];
        dte = dstr.get(1)[iblte.top] +
              dstr.get(2)[iblte.bottom] + ante;
        cte = (ctau.get(1)[iblte.top] *
                   thet.get(1)[iblte.top] +
               ctau.get(2)[iblte.bottom] *
                   thet.get(2)[iblte.bottom]) /
               tte;
        tesys(cte, tte, dte);

        tte_tte1 = 1.0;
        tte_tte2 = 1.0;
        dte_mte1 = 1.0 / uedg.top[iblte.top];
        dte_ute1 = -dstr.get(1)[iblte.top] /
                    uedg.top[iblte.top];
        dte_mte2 = 1.0 / uedg.bottom[iblte.bottom];
        dte_ute2 = -dstr.get(2)[iblte.bottom] /
                    uedg.bottom[iblte.bottom];
        cte_cte1 = thet.get(1)[iblte.top] / tte;
        cte_cte2 = thet.get(2)[iblte.bottom] / tte;
        cte_tte1 = (ctau.get(1)[iblte.top] - cte) / tte;
        cte_tte2 = (ctau.get(2)[iblte.bottom] - cte) / tte;

        //----- re-define d1 sensitivities wrt m since d1 depends on both te ds
        // values
      for (int js = 1; js <= 2; js++) {
        for (int jbl = 0; jbl < nbl.get(js) - 1; ++jbl) {
            int jv = isys.get(js)[jbl];
            d1_m[jv] = dte_ute1 * ute1_m[jv] + dte_ute2 * ute2_m[jv];
          }
        }
        d1_m[jvte1] = d1_m[jvte1] + dte_mte1;
        d1_m[jvte2] = d1_m[jvte2] + dte_mte2;

        //----- "forced" changes from  uedg --- usav=uinv+dij*mass	mismatch
        due1 = 0.0;
        dds1 =
            dte_ute1 *
                (uedg.top[iblte.top] -
                 usav.top[iblte.top]) +
            dte_ute2 *
                (uedg.bottom[iblte.bottom] -
                 usav.bottom[iblte.bottom]);
      } else {
        blsys();
      }

      //---- save wall shear and equil. max shear coefficient for plotting
      // output
      ctq.get(is)[ibl] = blData2.cqz.scalar;

      //---- set xi sensitivities wrt le ue changes
      if (is == 1) {
        xi_ule1 = sst_go;
        xi_ule2 = -sst_gp;
      } else {
        xi_ule1 = -sst_go;
        xi_ule2 = sst_gp;
      }

      //---- stuff bl system coefficients into main jacobian matrix

      for (int jv = 1; jv <= nsys; jv++) {
        vm[0][jv][iv] = blc.a1(0, 2) * d1_m[jv] + blc.a1(0, 3) * u1_m[jv] +
                        blc.a2(0, 2) * d2_m[jv] + blc.a2(0, 3) * u2_m[jv] +
                        (blc.a1(0, 4) + blc.a2(0, 4) + blc.d_xi[0]) *
                            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
      }

      vb[iv](0, 0) = blc.a1(0, 0);
      vb[iv](0, 1) = blc.a1(0, 1);

      va[iv](0, 0) = blc.a2(0, 0);
      va[iv](0, 1) = blc.a2(0, 1);

      if (lalfa)
        vdel[iv](0, 1) = blc.d_re[0] * re_clmr + blc.d_msq[0] * msq_clmr;
      else
        vdel[iv](0, 1) = (blc.a1(0, 3) * u1_a + blc.a1(0, 2) * d1_a) +
                         (blc.a2(0, 3) * u2_a + blc.a2(0, 2) * d2_a) +
                         (blc.a1(0, 4) + blc.a2(0, 4) + blc.d_xi[0]) *
                             (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

      vdel[iv](0, 0) = blc.rhs[0] + (blc.a1(0, 3) * due1 + blc.a1(0, 2) * dds1) +
                       (blc.a2(0, 3) * due2 + blc.a2(0, 2) * dds2) +
                       (blc.a1(0, 4) + blc.a2(0, 4) + blc.d_xi[0]) *
                           (xi_ule1 * dule1 + xi_ule2 * dule2);

      for (int jv = 1; jv <= nsys; jv++) {
        vm[1][jv][iv] = blc.a1(1, 2) * d1_m[jv] + blc.a1(1, 3) * u1_m[jv] +
                        blc.a2(1, 2) * d2_m[jv] + blc.a2(1, 3) * u2_m[jv] +
                        (blc.a1(1, 4) + blc.a2(1, 4) + blc.d_xi[1]) *
                            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
      }
      vb[iv](1, 0) = blc.a1(1, 0);
      vb[iv](1, 1) = blc.a1(1, 1);

      va[iv](1, 0) = blc.a2(1, 0);
      va[iv](1, 1) = blc.a2(1, 1);

      if (lalfa)
        vdel[iv](1, 1) = blc.d_re[1] * re_clmr + blc.d_msq[1] * msq_clmr;
      else
        vdel[iv](1, 1) = (blc.a1(1, 3) * u1_a + blc.a1(1, 2) * d1_a) +
                         (blc.a2(1, 3) * u2_a + blc.a2(1, 2) * d2_a) +
                         (blc.a1(1, 4) + blc.a2(1, 4) + blc.d_xi[1]) *
                             (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

      vdel[iv](1, 0) = blc.rhs[1] + (blc.a1(1, 3) * due1 + blc.a1(1, 2) * dds1) +
                       (blc.a2(1, 3) * due2 + blc.a2(1, 2) * dds2) +
                       (blc.a1(1, 4) + blc.a2(1, 4) + blc.d_xi[1]) *
                           (xi_ule1 * dule1 + xi_ule2 * dule2);

      // memory overlap problem
      for (int jv = 1; jv <= nsys; jv++) {
        vm[2][jv][iv] = blc.a1(2, 2) * d1_m[jv] + blc.a1(2, 3) * u1_m[jv] +
                        blc.a2(2, 2) * d2_m[jv] + blc.a2(2, 3) * u2_m[jv] +
                        (blc.a1(2, 4) + blc.a2(2, 4) + blc.d_xi[2]) *
                            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
      }

      vb[iv](2, 0) = blc.a1(2, 0);
      vb[iv](2, 1) = blc.a1(2, 1);

      va[iv](2, 0) = blc.a2(2, 0);
      va[iv](2, 1) = blc.a2(2, 1);

      if (lalfa)
        vdel[iv](2, 1) = blc.d_re[2] * re_clmr + blc.d_msq[2] * msq_clmr;
      else
        vdel[iv](2, 1) = (blc.a1(2, 3) * u1_a + blc.a1(2, 2) * d1_a) +
                         (blc.a2(2, 3) * u2_a + blc.a2(2, 2) * d2_a) +
                         (blc.a1(2, 4) + blc.a2(2, 4) + blc.d_xi[2]) *
                             (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

      vdel[iv](2, 0) = blc.rhs[2] + (blc.a1(2, 3) * due1 + blc.a1(2, 2) * dds1) +
                       (blc.a2(2, 3) * due2 + blc.a2(2, 2) * dds2) +
                       (blc.a1(2, 4) + blc.a2(2, 4) + blc.d_xi[2]) *
                           (xi_ule1 * dule1 + xi_ule2 * dule2);

      if (ibl == iblte.get(is) + 1) {
        //----- redefine coefficients for tte, dte, etc
        vz[0][0] = blc.a1(0, 0) * cte_cte1;
        vz[0][1] = blc.a1(0, 0) * cte_tte1 + blc.a1(0, 1) * tte_tte1;
        vb[iv](0, 0) = blc.a1(0, 0) * cte_cte2;
        vb[iv](0, 1) = blc.a1(0, 0) * cte_tte2 + blc.a1(0, 1) * tte_tte2;

        vz[1][0] = blc.a1(1, 0) * cte_cte1;
        vz[1][1] = blc.a1(1, 0) * cte_tte1 + blc.a1(1, 1) * tte_tte1;
        vb[iv](1, 0) = blc.a1(1, 0) * cte_cte2;
        vb[iv](1, 1) = blc.a1(1, 0) * cte_tte2 + blc.a1(1, 1) * tte_tte2;

        vz[2][0] = blc.a1(2, 0) * cte_cte1;
        vz[2][1] = blc.a1(2, 0) * cte_tte1 + blc.a1(2, 1) * tte_tte1;
        vb[iv](2, 0) = blc.a1(2, 0) * cte_cte2;
        vb[iv](2, 1) = blc.a1(2, 0) * cte_tte2 + blc.a1(2, 1) * tte_tte2;
      }

      //---- turbulent intervals will follow if currently at transition interval
      if (tran) {
        turb = true;

        //------ save transition location
        itran.get(is) = ibl;
      }

      tran = false;

      if (ibl == iblte.get(is)) {
        //----- set "2" variables at te to wake correlations for next station

        turb = true;
        wake = true;
        blvar(blData2, FlowRegimeEnum::Wake);
        blmid(FlowRegimeEnum::Wake);
      }
      u1_m = u2_m;
      d1_m = d2_m;

      u1_a = u2_a;
      d1_a = d2_a;

      due1 = due2;
      dds1 = dds2;

      //---- set bl variables for next station
      //			for (icom=1; icom<= ncom;icom++)
      // com1[icom] = com2[icom];
      stepbl();

      //---- next streamwise station
    }

    //---- next airfoil side
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

bool XFoil::stepbl() {
  blData1 = blData2;
  return true;
}

bool XFoil::stfind() {
  //-----------------------------------------
  //     locates stagnation point arc length
  //     location sst and panel index ist.
  //-----------------------------------------

  int i;
  bool bFound = false;
  for (i = 0; i < n - 1; i++) {
    if (surface_vortex(0, i) >= 0.0 && surface_vortex(0, i + 1) < 0.0) {
      bFound = true;
      break;
    }
  }

  if (!bFound) {
    writeString("stfind: Stagnation point not found. Continuing ...\n");
    i = n / 2;
  }

  i_stagnation = i;
  const double dgam = surface_vortex(0, i + 1) - surface_vortex(0, i);
  const double ds = spline_length[i + 1] - spline_length[i];

  //---- evaluate so as to minimize roundoff for very small gam[i] or gam[i+1]
  if (surface_vortex(0, i) < -surface_vortex(0, i + 1))
    sst = spline_length[i] - ds * (surface_vortex(0, i) / dgam);
  else
    sst = spline_length[i + 1] - ds * (surface_vortex(0, i + 1) / dgam);

  //---- tweak stagnation point if it falls right on a node (very unlikely)
  if (sst <= spline_length[i])
    sst = spline_length[i] + 0.0000001;
  if (sst >= spline_length[i + 1])
    sst = spline_length[i + 1] - 0.0000001;

  sst_go = (sst - spline_length[i + 1]) / dgam;
  sst_gp = (spline_length[i] - sst) / dgam;

  return true;
}

bool XFoil::stmove() {
  //--------------------------------------------------
  //    moves stagnation point location to new panel.
  //---------------------------------------------------
  int istold;
  //-- locate new stagnation point arc length sst from gam distribution
  istold = i_stagnation + INDEX_START_WITH;
  stfind();

  if (istold == i_stagnation + INDEX_START_WITH) {
    //--- recalculate new arc length array
    xicalc();
  } else {
    //--- set new bl position -> panel position  pointers
    iblpan();

    //--- set new inviscid bl edge velocity uinv from qinv
    uicalc();

    //--- recalculate new arc length array
    xicalc();

    //--- set  bl position -> system line  pointers
    iblsys();

    if (i_stagnation + INDEX_START_WITH > istold) {
      //---- increase in number of points on top side (is=1)
      int idif = i_stagnation + INDEX_START_WITH - istold;

      itran.top += idif;
      itran.bottom -= idif;

      //---- move top side bl variables downstream
      for (int ibl = nbl.top - 2; ibl >= idif; ibl--) {
        ctau.top[ibl] = ctau.top[ibl - idif];
        thet.top[ibl] = thet.top[ibl - idif];
        dstr.top[ibl] = dstr.top[ibl - idif];
        uedg.top[ibl] = uedg.top[ibl - idif];
      }

      //---- set bl variables between old and new stagnation point
      const double dudx =
          uedg.top[idif] / xssi.top[idif];
      for (int ibl = idif; ibl >= 1; ibl--) {
        ctau.top[ibl - 1] = ctau.top[idif];
        thet.top[ibl - 1] = thet.top[idif];
        dstr.top[ibl - 1] = dstr.top[idif];
        uedg.top[ibl - 1] = dudx * xssi.top[ibl - 1];
      }

      //---- move bottom side bl variables upstream
      for (int ibl = 0; ibl < nbl.bottom - 1; ibl++) {
        ctau.bottom[ibl] = ctau.bottom[ibl + idif];
        thet.bottom[ibl] = thet.bottom[ibl + idif];
        dstr.bottom[ibl] = dstr.bottom[ibl + idif];
        uedg.bottom[ibl] = uedg.bottom[ibl + idif];
      }
    } else {
      //---- increase in number of points on bottom side (is=2)
      int idif = istold - (i_stagnation + INDEX_START_WITH);

      itran.top = itran.top - idif;
      itran.bottom = itran.bottom + idif;

      //---- move bottom side bl variables downstream
      for (int ibl = nbl.bottom - 1; ibl >= idif + 1; ibl--) {
        ctau.bottom[ibl - 1] = ctau.bottom[(ibl - idif) - 1];
        thet.bottom[ibl - 1] = thet.bottom[(ibl - idif) - 1];
        dstr.bottom[ibl - 1] = dstr.bottom[(ibl - idif) - 1];
        uedg.bottom[ibl - 1] = uedg.bottom[(ibl - idif) - 1];
      }

      //---- set bl variables between old and new stagnation point
      const double dudx =
          uedg.bottom[idif] / xssi.bottom[idif];
      for (int ibl = idif; ibl >= 1; ibl--) {
        ctau.bottom[ibl - 1] = ctau.bottom[idif];
        thet.bottom[ibl - 1] = thet.bottom[idif];
        dstr.bottom[ibl - 1] = dstr.bottom[idif];
        uedg.bottom[ibl - 1] = dudx * xssi.bottom[ibl - 1];
      }

      //---- move top side bl variables upstream
      for (int ibl = 0; ibl < nbl.top - 1; ibl++) {
        ctau.top[ibl] = ctau.top[ibl + idif];
        thet.top[ibl] = thet.top[ibl + idif];
        dstr.top[ibl] = dstr.top[ibl + idif];
        uedg.top[ibl] = uedg.top[ibl + idif];
      }
    }
  }

  //-- set new mass array since ue has been tweaked
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      mass.get(is)[ibl] = dstr.get(is)[ibl] * uedg.get(is)[ibl];
    }
  }

  return true;
}

// moved to XFoil_geometry.cpp: tecalc()

bool XFoil::tesys(double cte, double tte, double dte) {
  //--------------------------------------------------------
  //	   sets up "dummy" bl system between airfoil te point
  //	   and first wake point infinitesimally behind te.
  //--------------------------------------------------------

  blc.clear();

  blvar(blData2, FlowRegimeEnum::Wake);

  blc.a1(0, 0) = -1.0;
  blc.a2(0, 0) = 1.0;
  blc.rhs[0] = cte - blData2.param.sz;

  blc.a1(1, 1) = -1.0;
  blc.a2(1, 1) = 1.0;
  blc.rhs[1] = tte - blData2.param.tz;

  blc.a1(2, 2) = -1.0;
  blc.a2(2, 2) = 1.0;
  blc.rhs[2] = dte - blData2.param.dz - blData2.param.dwz;

  return true;
}

bool XFoil::trchek() {
  //----------------------------------------------------------------
  //     new second-order version:  december 1994.
  //
  //     checks if transition occurs in the current interval x1..x2.
  //     if transition occurs, then set transition location xt, and
  //     its sensitivities to "1" and "2" variables.  if no transition,
  //     set amplification ampl2.
  //
  //     solves the implicit amplification equation for n2:
  //
  //       n2 - n1     n'(xt,nt) + n'(x1,n1)
  //       -------  =  ---------------------
  //       x2 - x1               2
  //
  //     in effect, a 2-point central difference is used between
  //     x1..x2 (no transition), or x1..xt (transition).  the switch
  //     is done by defining xt,nt in the equation above depending
  //     on whether n2 exceeds ncrit.
  //
  //  if n2<ncrit:  nt=n2    , xt=x2                  (no transition)
  //
  //  if n2>ncrit:  nt=ncrit , xt=(ncrit-n1)/(n2-n1)  (transition)
  //
  //----------------------------------------------------------------
  double amplt, sfa, sfa_a1, sfa_a2, sfx;
  double sfx_x1, sfx_x2, sfx_xf;
  double tt, dt, ut, amsave;
  double res = 0.0, res_a2 = 0.0;
  double da2 = 0.0, dxt = 0.0, tt_t1 = 0.0, dt_d1 = 0.0, ut_u1 = 0.0;
  double tt_t2 = 0.0, dt_d2 = 0.0, ut_u2 = 0.0, tt_a1 = 0.0, dt_a1 = 0.0;
  double ut_a1 = 0.0, tt_x1 = 0.0, dt_x1 = 0.0, ut_x1 = 0.0, tt_x2 = 0.0,
         dt_x2 = 0.0, ut_x2 = 0.0;
  double z_ax = 0.0, z_a1 = 0.0, z_t1 = 0.0, z_d1 = 0.0, z_u1 = 0.0, z_x1 = 0.0,
         z_a2 = 0.0, z_t2 = 0.0, z_d2 = 0.0, z_u2 = 0.0, z_x2 = 0.0, z_ms = 0.0,
         z_re = 0.0;
  double amplt_a2, wf, wf_a1, wf_a2, wf_xf, wf_x1, wf_x2;
  double xt_a2, dt_a2, tt_a2;
  double ut_a2;
  double daeps = 0.00005;

  amplt_a2 = 0.0;
  xt_a2 = dt_a2 = tt_a2 = 0.0;
  ut_a2 = 0.0;

  //---- save variables and sensitivities at ibl ("2") for future restoration
  saveblData(2);

  //---- calculate average amplification rate ax over x1..x2 interval
  AxResult ax_result =
      axset(blData1.hkz.scalar, blData1.param.tz, blData1.rtz.scalar,
            blData1.param.amplz, blData2.hkz.scalar, blData2.param.tz,
            blData2.rtz.scalar, blData2.param.amplz, amcrit);

  //---- set initial guess for iterate n2 (ampl2) at x2
  blData2.param.amplz = blData1.param.amplz +
                        ax_result.ax * (blData2.param.xz - blData1.param.xz);
  //---- solve implicit system for amplification ampl2
  auto iterateAmplification = [&]() -> bool {
    for (int itam = 0; itam < 30; itam++) {
      //---- define weighting factors wf1,wf2 for defining "t" quantities
      if (blData2.param.amplz <= amcrit) {
        //------ there is no transition yet,  "t" is the same as "2"
        amplt = blData2.param.amplz;
        amplt_a2 = 1.0;
        sfa = 1.0;
        sfa_a1 = 0.0;
        sfa_a2 = 0.0;
      } else {
        //------ there is transition in x1..x2, "t" is set from n1, n2
        amplt = amcrit;
        amplt_a2 = 0.0;
        sfa = (amplt - blData1.param.amplz) /
              (blData2.param.amplz - blData1.param.amplz);
        sfa_a1 = (sfa - 1.0) / (blData2.param.amplz - blData1.param.amplz);
        sfa_a2 = (-sfa) / (blData2.param.amplz - blData1.param.amplz);
      }

      if (xiforc < blData2.param.xz) {
        sfx = (xiforc - blData1.param.xz) /
              (blData2.param.xz - blData1.param.xz);
        sfx_x1 = (sfx - 1.0) / (blData2.param.xz - blData1.param.xz);
        sfx_x2 = (-sfx) / (blData2.param.xz - blData1.param.xz);
        sfx_xf = 1.0 / (blData2.param.xz - blData1.param.xz);
      } else {
        sfx = 1.0;
        sfx_x1 = 0.0;
        sfx_x2 = 0.0;
        sfx_xf = 0.0;
      }

      //---- set weighting factor from free or forced transition
      if (sfa < sfx) {
        wf = sfa;
        wf_a1 = sfa_a1;
        wf_a2 = sfa_a2;
        wf_x1 = 0.0;
        wf_x2 = 0.0;
        wf_xf = 0.0;
      } else {
        wf = sfx;
        wf_a1 = 0.0;
        wf_a2 = 0.0;
        wf_x1 = sfx_x1;
        wf_x2 = sfx_x2;
        wf_xf = sfx_xf;
      }

      //---- interpolate bl variables to xt
      xt = blData1.param.xz * (1 - wf) + blData2.param.xz * wf;
      tt = blData1.param.tz * (1 - wf) + blData2.param.tz * wf;
      dt = blData1.param.dz * (1 - wf) + blData2.param.dz * wf;
      ut = blData1.param.uz * (1 - wf) + blData2.param.uz * wf;

      xt_a2 = (blData2.param.xz - blData1.param.xz) * wf_a2;
      tt_a2 = (blData2.param.tz - blData1.param.tz) * wf_a2;
      dt_a2 = (blData2.param.dz - blData1.param.dz) * wf_a2;
      ut_a2 = (blData2.param.uz - blData1.param.uz) * wf_a2;

      //---- temporarily set "2" variables from "t" for blkin
      blData2.param.xz = xt;
      blData2.param.tz = tt;
      blData2.param.dz = dt;
      blData2.param.uz = ut;

      //---- calculate laminar secondary "t" variables hkt, rtt
      blkin();

      blData::blVector hkt = blData2.hkz;
      blData::blVector rtt = blData2.rtz;

      //---- restore clobbered "2" variables, except for ampl2
      amsave = blData2.param.amplz;

      restoreblData(2);

      blData2.param.amplz = amsave;

      //---- calculate amplification rate ax over current x1-xt interval
      ax_result = axset(blData1.hkz.scalar, blData1.param.tz,
                        blData1.rtz.scalar, blData1.param.amplz, hkt.scalar, tt, rtt.scalar,
                        amplt, amcrit);

      //---- punch out early if there is no amplification here
      if (ax_result.ax <= 0.0) {
        return true;
      }

      //---- set sensitivity of ax(a2)
      ax_result.ax_a2 =
          (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
           ax_result.ax_rt2 * rtt.t()) *
              tt_a2 +
          (ax_result.ax_hk2 * hkt.d()) * dt_a2 +
          (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_a2 +
          ax_result.ax_a2 * amplt_a2;

      //---- residual for implicit ampl2 definition (amplification equation)
      res = blData2.param.amplz - blData1.param.amplz -
            ax_result.ax * (blData2.param.xz - blData1.param.xz);
      res_a2 = 1.0 - ax_result.ax_a2 * (blData2.param.xz - blData1.param.xz);

      da2 = -res / res_a2;

      rlx = 1.0;
      dxt = xt_a2 * da2;

      if (rlx * fabs(dxt / (blData2.param.xz - blData1.param.xz)) > 0.05) {
        rlx = 0.05 * fabs((blData2.param.xz - blData1.param.xz) / dxt);
      }

      if (rlx * fabs(da2) > 1.0) {
        rlx = 1.0 * fabs(1.0 / da2);
      }

      //---- check if converged
      if (fabs(da2) < daeps) {
        return true;
      }

      if ((blData2.param.amplz > amcrit &&
           blData2.param.amplz + rlx * da2 < amcrit) ||
          (blData2.param.amplz < amcrit &&
           blData2.param.amplz + rlx * da2 > amcrit)) {
        //------ limited newton step so ampl2 doesn't step across amcrit either
        // way
        blData2.param.amplz = amcrit;
      } else {
        //------ regular newton step
        blData2.param.amplz = blData2.param.amplz + rlx * da2;
      }
    }
    return false;
  };

  if (!iterateAmplification()) {
    // TRACE("trchek2 - n2 convergence failed\n");
    writeString("trchek2 - n2 convergence failed\n");
    if (s_bCancel)
      return false;
  }

  //---- test for free or forced transition
  trfree = (blData2.param.amplz >= amcrit);
  trforc = (xiforc > blData1.param.xz) && (xiforc <= blData2.param.xz);

  //---- set transition interval flag
  tran = (trforc || trfree);

  if (!tran)
    return false;

  //---- resolve if both forced and free transition
  if (trfree && trforc) {
    trforc = xiforc < xt;
    trfree = xiforc >= xt;
  }

  if (trforc) {
    //----- if forced transition, then xt is prescribed,
    //-     no sense calculating the sensitivities, since we know them...
    xt = xiforc;
    xt_a1 = 0.0;
    xt_x1 = 0.0;
    xt_t1 = 0.0;
    xt_d1 = 0.0;
    xt_u1 = 0.0;
    xt_x2 = 0.0;
    xt_t2 = 0.0;
    xt_d2 = 0.0;
    xt_u2 = 0.0;
    xt_ms = 0.0;
    xt_re = 0.0;
    xt_xf = 1.0;
    return true;
  }

  //---- free transition ... set sensitivities of xt

  xt_x1 = (1 - wf);
  tt_t1 = (1 - wf);
  dt_d1 = (1 - wf);
  ut_u1 = (1 - wf);

  xt_x2 = wf;
  tt_t2 = wf;
  dt_d2 = wf;
  ut_u2 = wf;

  xt_a1 = (blData2.param.xz - blData1.param.xz) * wf_a1;
  tt_a1 = (blData2.param.tz - blData1.param.tz) * wf_a1;
  dt_a1 = (blData2.param.dz - blData1.param.dz) * wf_a1;
  ut_a1 = (blData2.param.uz - blData1.param.uz) * wf_a1;

  xt_x1 += (blData2.param.xz - blData1.param.xz) * wf_x1;
  tt_x1 = (blData2.param.tz - blData1.param.tz) * wf_x1;
  dt_x1 = (blData2.param.dz - blData1.param.dz) * wf_x1;
  ut_x1 = (blData2.param.uz - blData1.param.uz) * wf_x1;

  xt_x2 += (blData2.param.xz - blData1.param.xz) * wf_x2;
  tt_x2 = (blData2.param.tz - blData1.param.tz) * wf_x2;
  dt_x2 = (blData2.param.dz - blData1.param.dz) * wf_x2;
  ut_x2 = (blData2.param.uz - blData1.param.uz) * wf_x2;

  xt_xf = (blData2.param.xz - blData1.param.xz) * wf_xf;

  //---- at this point, ax = ax( hk1, t1, rt1, a1, hkt, tt, rtt, at )
  blData::blVector hkt = blData2.hkz;
  blData::blVector rtt = blData2.rtz;

  //---- set sensitivities of ax( t1 d1 u1 a1 t2 d2 u2 a2 ms re )
  double ax_t1 = ax_result.ax_hk1 * blData1.hkz.t() + ax_result.ax_t1 +
                 ax_result.ax_rt1 * blData1.rtz.t() +
                 (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
                  ax_result.ax_rt2 * rtt.t()) *
                     tt_t1;
  double ax_d1 =
      ax_result.ax_hk1 * blData1.hkz.d() + (ax_result.ax_hk2 * hkt.d()) * dt_d1;
  double ax_u1 =
      ax_result.ax_hk1 * blData1.hkz.u() + ax_result.ax_rt1 * blData1.rtz.u() +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_u1;
  double ax_a1 =
      ax_result.ax_a1 +
      (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
       ax_result.ax_rt2 * rtt.t()) *
          tt_a1 +
      (ax_result.ax_hk2 * hkt.d()) * dt_a1 +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_a1;
  double ax_x1 =
      (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
       ax_result.ax_rt2 * rtt.t()) *
          tt_x1 +
      (ax_result.ax_hk2 * hkt.d()) * dt_x1 +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_x1;

  double ax_t2 = (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
                  ax_result.ax_rt2 * rtt.t()) *
                 tt_t2;
  double ax_d2 = (ax_result.ax_hk2 * hkt.d()) * dt_d2;
  double ax_u2 =
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_u2;
  double ax_a2 =
      ax_result.ax_a2 * amplt_a2 +
      (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
       ax_result.ax_rt2 * rtt.t()) *
          tt_a2 +
      (ax_result.ax_hk2 * hkt.d()) * dt_a2 +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_a2;
  double ax_x2 =
      (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
       ax_result.ax_rt2 * rtt.t()) *
          tt_x2 +
      (ax_result.ax_hk2 * hkt.d()) * dt_x2 +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_x2;

  double ax_ms = ax_result.ax_hk2 * hkt.ms() + ax_result.ax_rt2 * rtt.ms() +
                 ax_result.ax_hk1 * blData1.hkz.ms() +
                 ax_result.ax_rt1 * blData1.rtz.ms();
  double ax_re =
      ax_result.ax_rt2 * rtt.re() + ax_result.ax_rt1 * blData1.rtz.re();

  //---- set sensitivities of residual res
  z_ax = -(blData2.param.xz - blData1.param.xz);

  z_a1 = z_ax * ax_a1 - 1.0;
  z_t1 = z_ax * ax_t1;
  z_d1 = z_ax * ax_d1;
  z_u1 = z_ax * ax_u1;
  z_x1 = z_ax * ax_x1 + ax_result.ax;

  z_a2 = z_ax * ax_a2 + 1.0;
  z_t2 = z_ax * ax_t2;
  z_d2 = z_ax * ax_d2;
  z_u2 = z_ax * ax_u2;
  z_x2 = z_ax * ax_x2 - ax_result.ax;

  z_ms = z_ax * ax_ms;
  z_re = z_ax * ax_re;

  //---- set sensitivities of xt, with res being stationary for a2 constraint
  xt_a1 = xt_a1 - (xt_a2 / z_a2) * z_a1;
  xt_t1 = -(xt_a2 / z_a2) * z_t1;
  xt_d1 = -(xt_a2 / z_a2) * z_d1;
  xt_u1 = -(xt_a2 / z_a2) * z_u1;
  xt_x1 = xt_x1 - (xt_a2 / z_a2) * z_x1;
  xt_t2 = -(xt_a2 / z_a2) * z_t2;
  xt_d2 = -(xt_a2 / z_a2) * z_d2;
  xt_u2 = -(xt_a2 / z_a2) * z_u2;
  xt_x2 = xt_x2 - (xt_a2 / z_a2) * z_x2;
  xt_ms = -(xt_a2 / z_a2) * z_ms;
  xt_re = -(xt_a2 / z_a2) * z_re;
  xt_xf = 0.0;

  return true;
}

bool XFoil::trdif() {
  //-----------------------------------------------
  //     sets up the newton system governing the
  //     transition interval.  equations governing
  //     the  laminar  part  x1 < xi < xt  and
  //     the turbulent part  xt < xi < x2
  //     are simply summed.
  //-----------------------------------------------
  Matrix<double, 4, 5> bl1, bl2, bt1, bt2;
  Vector<double, 4> blrez, blm, blr, blx, btrez, btm, btr, btx;

  double tt, tt_a1, tt_x1, tt_x2, tt_t1, tt_t2, tt_d1, tt_d2, tt_u1, tt_u2;
  double tt_ms, tt_re, tt_xf, dt, dt_a1, dt_x1, dt_x2, dt_t1, dt_t2;
  double dt_d1, dt_d2, dt_u1, dt_u2, dt_ms, dt_re, dt_xf;
  double ut, ut_a1, ut_x1, ut_x2, ut_t1, ut_t2, ut_d1, ut_d2, ut_u1, ut_u2;
  double ut_ms, ut_re, ut_xf;
  double st, st_tt, st_dt, st_ut, st_ms, st_re, st_a1, st_x1, st_x2, st_t1,
      st_t2;
  double st_d1, st_d2, st_u1, st_u2, st_xf;
  double ctr, ctr_hk2;

  saveblData(1);
  saveblData(2);

  //---- weighting factors for linear interpolation to transition point
  double wf2 = (xt - blData1.param.xz) / (blData2.param.xz - blData1.param.xz);
  double wf2_xt = 1.0 / (blData2.param.xz - blData1.param.xz);

  double wf2_a1 = wf2_xt * xt_a1;
  double wf2_x1 =
      wf2_xt * xt_x1 + (wf2 - 1.0) / (blData2.param.xz - blData1.param.xz);
  double wf2_x2 = wf2_xt * xt_x2 - wf2 / (blData2.param.xz - blData1.param.xz);
  double wf2_t1 = wf2_xt * xt_t1;
  double wf2_t2 = wf2_xt * xt_t2;
  double wf2_d1 = wf2_xt * xt_d1;
  double wf2_d2 = wf2_xt * xt_d2;
  double wf2_u1 = wf2_xt * xt_u1;
  double wf2_u2 = wf2_xt * xt_u2;
  double wf2_ms = wf2_xt * xt_ms;
  double wf2_re = wf2_xt * xt_re;
  double wf2_xf = wf2_xt * xt_xf;

  double wf1 = 1.0 - wf2;
  double wf1_a1 = -wf2_a1;
  double wf1_x1 = -wf2_x1;
  double wf1_x2 = -wf2_x2;
  double wf1_t1 = -wf2_t1;
  double wf1_t2 = -wf2_t2;
  double wf1_d1 = -wf2_d1;
  double wf1_d2 = -wf2_d2;
  double wf1_u1 = -wf2_u1;
  double wf1_u2 = -wf2_u2;
  double wf1_ms = -wf2_ms;
  double wf1_re = -wf2_re;
  double wf1_xf = -wf2_xf;

  //-----interpolate primary variables to transition point
  tt = blData1.param.tz * wf1 + blData2.param.tz * wf2;
  tt_a1 = blData1.param.tz * wf1_a1 + blData2.param.tz * wf2_a1;
  tt_x1 = blData1.param.tz * wf1_x1 + blData2.param.tz * wf2_x1;
  tt_x2 = blData1.param.tz * wf1_x2 + blData2.param.tz * wf2_x2;
  tt_t1 = blData1.param.tz * wf1_t1 + blData2.param.tz * wf2_t1 + wf1;
  tt_t2 = blData1.param.tz * wf1_t2 + blData2.param.tz * wf2_t2 + wf2;
  tt_d1 = blData1.param.tz * wf1_d1 + blData2.param.tz * wf2_d1;
  tt_d2 = blData1.param.tz * wf1_d2 + blData2.param.tz * wf2_d2;
  tt_u1 = blData1.param.tz * wf1_u1 + blData2.param.tz * wf2_u1;
  tt_u2 = blData1.param.tz * wf1_u2 + blData2.param.tz * wf2_u2;
  tt_ms = blData1.param.tz * wf1_ms + blData2.param.tz * wf2_ms;
  tt_re = blData1.param.tz * wf1_re + blData2.param.tz * wf2_re;
  tt_xf = blData1.param.tz * wf1_xf + blData2.param.tz * wf2_xf;

  dt = blData1.param.dz * wf1 + blData2.param.dz * wf2;
  dt_a1 = blData1.param.dz * wf1_a1 + blData2.param.dz * wf2_a1;
  dt_x1 = blData1.param.dz * wf1_x1 + blData2.param.dz * wf2_x1;
  dt_x2 = blData1.param.dz * wf1_x2 + blData2.param.dz * wf2_x2;
  dt_t1 = blData1.param.dz * wf1_t1 + blData2.param.dz * wf2_t1;
  dt_t2 = blData1.param.dz * wf1_t2 + blData2.param.dz * wf2_t2;
  dt_d1 = blData1.param.dz * wf1_d1 + blData2.param.dz * wf2_d1 + wf1;
  dt_d2 = blData1.param.dz * wf1_d2 + blData2.param.dz * wf2_d2 + wf2;
  dt_u1 = blData1.param.dz * wf1_u1 + blData2.param.dz * wf2_u1;
  dt_u2 = blData1.param.dz * wf1_u2 + blData2.param.dz * wf2_u2;
  dt_ms = blData1.param.dz * wf1_ms + blData2.param.dz * wf2_ms;
  dt_re = blData1.param.dz * wf1_re + blData2.param.dz * wf2_re;
  dt_xf = blData1.param.dz * wf1_xf + blData2.param.dz * wf2_xf;

  ut = blData1.param.uz * wf1 + blData2.param.uz * wf2;
  ut_a1 = blData1.param.uz * wf1_a1 + blData2.param.uz * wf2_a1;
  ut_x1 = blData1.param.uz * wf1_x1 + blData2.param.uz * wf2_x1;
  ut_x2 = blData1.param.uz * wf1_x2 + blData2.param.uz * wf2_x2;
  ut_t1 = blData1.param.uz * wf1_t1 + blData2.param.uz * wf2_t1;
  ut_t2 = blData1.param.uz * wf1_t2 + blData2.param.uz * wf2_t2;
  ut_d1 = blData1.param.uz * wf1_d1 + blData2.param.uz * wf2_d1;
  ut_d2 = blData1.param.uz * wf1_d2 + blData2.param.uz * wf2_d2;
  ut_u1 = blData1.param.uz * wf1_u1 + blData2.param.uz * wf2_u1 + wf1;
  ut_u2 = blData1.param.uz * wf1_u2 + blData2.param.uz * wf2_u2 + wf2;
  ut_ms = blData1.param.uz * wf1_ms + blData2.param.uz * wf2_ms;
  ut_re = blData1.param.uz * wf1_re + blData2.param.uz * wf2_re;
  ut_xf = blData1.param.uz * wf1_xf + blData2.param.uz * wf2_xf;

  //---- set primary "t" variables at xt  (really placed into "2" variables)
  blData2.param.xz = xt;
  blData2.param.tz = tt;
  blData2.param.dz = dt;
  blData2.param.uz = ut;

  blData2.param.amplz = amcrit;
  blData2.param.sz = 0.0;

  //---- calculate laminar secondary "t" variables
  blkin();
  blvar(blData2, FlowRegimeEnum::Laminar);

  //---- calculate x1-xt midpoint cfm value
  blmid(FlowRegimeEnum::Laminar);

  //=    at this point, all "2" variables are really "t" variables at xt

  //---- set up newton system for dam, dth, dds, due, dxi  at  x1 and xt
  bldif(1);

  //---- the current newton system is in terms of "1" and "t" variables,
  //-    so calculate its equivalent in terms of "1" and "2" variables.
  //-    in other words, convert residual sensitivities wrt "t" variables
  //-    into sensitivities wrt "1" and "2" variables.  the amplification
  //-    equation is unnecessary here, so the k=1 row is left empty.
  for (int k = 1; k < 3; k++) {
    blrez[k] = blc.rhs[k];
    blm[k] = blc.d_msq[k] + blc.a2(k, 1) * tt_ms + blc.a2(k, 2) * dt_ms +
             blc.a2(k, 3) * ut_ms + blc.a2(k, 4) * xt_ms;
    blr[k] = blc.d_re[k] + blc.a2(k, 1) * tt_re + blc.a2(k, 2) * dt_re +
             blc.a2(k, 3) * ut_re + blc.a2(k, 4) * xt_re;
    blx[k] = blc.d_xi[k] + blc.a2(k, 1) * tt_xf + blc.a2(k, 2) * dt_xf +
             blc.a2(k, 3) * ut_xf + blc.a2(k, 4) * xt_xf;

    bl1(k, 0) = blc.a1(k, 0) + blc.a2(k, 1) * tt_a1 + blc.a2(k, 2) * dt_a1 +
                blc.a2(k, 3) * ut_a1 + blc.a2(k, 4) * xt_a1;
    bl1(k, 1) = blc.a1(k, 1) + blc.a2(k, 1) * tt_t1 + blc.a2(k, 2) * dt_t1 +
                blc.a2(k, 3) * ut_t1 + blc.a2(k, 4) * xt_t1;
    bl1(k, 2) = blc.a1(k, 2) + blc.a2(k, 1) * tt_d1 + blc.a2(k, 2) * dt_d1 +
                blc.a2(k, 3) * ut_d1 + blc.a2(k, 4) * xt_d1;
    bl1(k, 3) = blc.a1(k, 3) + blc.a2(k, 1) * tt_u1 + blc.a2(k, 2) * dt_u1 +
                blc.a2(k, 3) * ut_u1 + blc.a2(k, 4) * xt_u1;
    bl1(k, 4) = blc.a1(k, 4) + blc.a2(k, 1) * tt_x1 + blc.a2(k, 2) * dt_x1 +
                blc.a2(k, 3) * ut_x1 + blc.a2(k, 4) * xt_x1;

    bl2(k, 0) = 0.0;
    bl2(k, 1) = blc.a2(k, 1) * tt_t2 + blc.a2(k, 2) * dt_t2 + blc.a2(k, 3) * ut_t2 +
                blc.a2(k, 4) * xt_t2;
    bl2(k, 2) = blc.a2(k, 1) * tt_d2 + blc.a2(k, 2) * dt_d2 + blc.a2(k, 3) * ut_d2 +
                blc.a2(k, 4) * xt_d2;
    bl2(k, 3) = blc.a2(k, 1) * tt_u2 + blc.a2(k, 2) * dt_u2 + blc.a2(k, 3) * ut_u2 +
                blc.a2(k, 4) * xt_u2;
    bl2(k, 4) = blc.a2(k, 1) * tt_x2 + blc.a2(k, 2) * dt_x2 + blc.a2(k, 3) * ut_x2 +
                blc.a2(k, 4) * xt_x2;
  }

  //**** second, set up turbulent part between xt and x2  ****

  //---- calculate equilibrium shear coefficient cqt at transition point
  blvar(blData2, FlowRegimeEnum::Turbulent);

  //---- set initial shear coefficient value st at transition point
  //-    ( note that cq2, cq2_t2, etc. are really "cqt", "cqt_tt", etc.)

  ctr = 1.8 * exp(-3.3 / (blData2.hkz.scalar - 1.0));
  ctr_hk2 = ctr * 3.3 / (blData2.hkz.scalar - 1.0) / (blData2.hkz.scalar - 1.0);

  st = ctr * blData2.cqz.scalar;
  st_tt =
      ctr * blData2.cqz.t() + blData2.cqz.scalar * ctr_hk2 * blData2.hkz.t();
  st_dt =
      ctr * blData2.cqz.d() + blData2.cqz.scalar * ctr_hk2 * blData2.hkz.d();
  st_ut =
      ctr * blData2.cqz.u() + blData2.cqz.scalar * ctr_hk2 * blData2.hkz.u();
  st_ms =
      ctr * blData2.cqz.ms() + blData2.cqz.scalar * ctr_hk2 * blData2.hkz.ms();
  st_re = ctr * blData2.cqz.re();

  //---- calculate st sensitivities wrt the actual "1" and "2" variables
  st_a1 = st_tt * tt_a1 + st_dt * dt_a1 + st_ut * ut_a1;
  st_x1 = st_tt * tt_x1 + st_dt * dt_x1 + st_ut * ut_x1;
  st_x2 = st_tt * tt_x2 + st_dt * dt_x2 + st_ut * ut_x2;
  st_t1 = st_tt * tt_t1 + st_dt * dt_t1 + st_ut * ut_t1;
  st_t2 = st_tt * tt_t2 + st_dt * dt_t2 + st_ut * ut_t2;
  st_d1 = st_tt * tt_d1 + st_dt * dt_d1 + st_ut * ut_d1;
  st_d2 = st_tt * tt_d2 + st_dt * dt_d2 + st_ut * ut_d2;
  st_u1 = st_tt * tt_u1 + st_dt * dt_u1 + st_ut * ut_u1;
  st_u2 = st_tt * tt_u2 + st_dt * dt_u2 + st_ut * ut_u2;
  st_ms = st_tt * tt_ms + st_dt * dt_ms + st_ut * ut_ms + st_ms;
  st_re = st_tt * tt_re + st_dt * dt_re + st_ut * ut_re + st_re;
  st_xf = st_tt * tt_xf + st_dt * dt_xf + st_ut * ut_xf;

  blData2.param.amplz = 0.0;
  blData2.param.sz = st;

  //---- recalculate turbulent secondary "t" variables using proper cti
  blvar(blData2, FlowRegimeEnum::Turbulent);

  stepbl();
  restoreblData(2);

  //---- calculate xt-x2 midpoint cfm value
  blmid(FlowRegimeEnum::Turbulent);

  //---- set up newton system for dct, dth, dds, due, dxi  at  xt and x2
  bldif(2);

  //---- convert sensitivities wrt "t" variables into sensitivities
  //-    wrt "1" and "2" variables as done before for the laminar part
  Matrix<double, 5, 5> bt1_right =
      Matrix<double, 5, 5>{{st_a1, st_t1, st_d1, st_u1, st_x1},
                           {tt_a1, tt_t1, tt_d1, tt_u1, tt_x1},
                           {dt_a1, dt_t1, dt_d1, dt_u1, dt_x1},
                           {ut_a1, ut_t1, ut_d1, ut_u1, ut_x1},
                           {xt_a1, xt_t1, xt_d1, xt_u1, xt_x1}};

  Matrix<double, 5, 5> bt2_right =
      Matrix<double, 5, 5>{{0, st_t2, st_d2, st_u2, st_x2},
                           {0, tt_t2, tt_d2, tt_u2, tt_x2},
                           {0, dt_t2, dt_d2, dt_u2, dt_x2},
                           {0, ut_t2, ut_d2, ut_u2, ut_x2},
                           {0, xt_t2, xt_d2, xt_u2, xt_x2}};
  bt1.block(0, 0, 3, 5) = blc.a1.block(0, 0, 3, 5) * bt1_right;
  bt2.block(0, 0, 3, 5) = blc.a1.block(0, 0, 3, 5) * bt2_right;
  bt2 += blc.a2;
  for (int k = 0; k < 3; k++) {
    btrez[k] = blc.rhs[k];
    btm[k] = blc.d_msq[k] + blc.a1(k, 0) * st_ms + blc.a1(k, 1) * tt_ms +
             blc.a1(k, 2) * dt_ms + blc.a1(k, 3) * ut_ms + blc.a1(k, 4) * xt_ms;
    btr[k] = blc.d_re[k] + blc.a1(k, 0) * st_re + blc.a1(k, 1) * tt_re +
             blc.a1(k, 2) * dt_re + blc.a1(k, 3) * ut_re + blc.a1(k, 4) * xt_re;
    btx[k] = blc.d_xi[k] + blc.a1(k, 0) * st_xf + blc.a1(k, 1) * tt_xf +
             blc.a1(k, 2) * dt_xf + blc.a1(k, 3) * ut_xf + blc.a1(k, 4) * xt_xf;
  }

  //---- add up laminar and turbulent parts to get final system
  //-    in terms of honest-to-god "1" and "2" variables.
  blc.rhs[0] = btrez[0];
  blc.rhs[1] = blrez[1] + btrez[1];
  blc.rhs[2] = blrez[2] + btrez[2];
  blc.d_msq[0] = btm[0];
  blc.d_msq[1] = blm[1] + btm[1];
  blc.d_msq[2] = blm[2] + btm[2];
  blc.d_re[0] = btr[0];
  blc.d_re[1] = blr[1] + btr[1];
  blc.d_re[2] = blr[2] + btr[2];
  blc.d_xi[0] = btx[0];
  blc.d_xi[1] = blx[1] + btx[1];
  blc.d_xi[2] = blx[2] + btx[2];
  blc.a1.row(0) = bt1.row(0);
  blc.a2.row(0) = bt2.row(0);
  blc.a1.middleRows(1, 2) = bl1.middleRows(1, 2) + bt1.middleRows(1, 2);
  blc.a2.middleRows(1, 2) = bl2.middleRows(1, 2) + bt2.middleRows(1, 2);

  //---- to be sanitary, restore "1" quantities which got clobbered
  //-    in all of the numerical gymnastics above.  the "2" variables
  //-    were already restored for the xt-x2 differencing part.
  //	for (icom=1; icom<=ncom;icom++){
  //		com1[icom] = c1sav[icom];
  //	}
  restoreblData(1);

  return true;
}

bool XFoil::ueset() {
  //---------------------------------------------------------
  //     sets ue from inviscid ue plus all source influence
  //---------------------------------------------------------
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      double dui = 0.0;
      for (int js = 1; js <= 2; js++) {
        for (int jbl = 0; jbl < nbl.get(js) - 1; ++jbl) {
          double ue_m = -vti.get(is)[ibl] * vti.get(js)[jbl] *
                        dij(ipan.get(is)[ibl],
                            ipan.get(js)[jbl]);
          dui += ue_m * mass.get(js)[jbl];
        }
      }
      uedg.get(is)[ibl] = uinv.get(is)[ibl] + dui;
    }
  }
  return true;
}

bool XFoil::uicalc() {
  //--------------------------------------------------------------
  //     sets inviscid ue from panel inviscid tangential velocity
  //--------------------------------------------------------------
  for (int is = 1; is <= 2; is++) {
    uinv.get(is)[0] = 0.0;
    uinv_a.get(is)[0] = 0.0;
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      int i = ipan.get(is)[ibl];
      uinv.get(is)[ibl] = vti.get(is)[ibl] * qinv[i];
      uinv_a.get(is)[ibl] = vti.get(is)[ibl] * qinv_a[i];
    }
  }

  return true;
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
                        const SidePair<VectorXd> &u_ac, double qnew[],
                        double q_ac[]) {
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
void XFoil::computeClFromQtan(const double qnew[], const double q_ac[],
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

bool XFoil::update() {
  //------------------------------------------------------------------
  //      adds on newton deltas to boundary layer variables.
  //      checks for excessive changes and underrelaxes if necessary.
  //      calculates max and rms changes.
  //      also calculates the change in the global variable "ac".
  //        if lalfa=true , "ac" is cl
  //        if lalfa=false, "ac" is alpha
  //------------------------------------------------------------------

  // int i = 0, is = 0, iv, iw, j, js, jv, ibl, jbl, kbl = 0;

  SidePair<VectorXd> unew, u_ac;
  unew.top = VectorXd::Zero(IVX);
  unew.bottom = VectorXd::Zero(IVX);
  u_ac.top = VectorXd::Zero(IVX);
  u_ac.bottom = VectorXd::Zero(IVX);

  double qnew[IQX], q_ac[IQX];
  memset(qnew, 0, IQX * sizeof(double));
  memset(q_ac, 0, IQX * sizeof(double));
  double dalmax = 0.0, dalmin = 0.0, dclmax = 0.0, dclmin = 0.0;
  double dac = 0.0, dhi = 0.0, dlo = 0.0;
  double dswaki, hklim, msq, dsw;
  double clnew = 0.0, cl_a = 0.0, cl_ms = 0.0, cl_ac = 0.0;

  //---- max allowable alpha changes per iteration
  dalmax = 0.5 * dtor;
  dalmin = -0.5 * dtor;
  //---- max allowable cl change per iteration
  dclmax = 0.5;
  dclmin = -0.5;
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

  if (lalfa) {
    //===== alpha is prescribed: ac is cl

    //---- set change in re to account for cl changing, since re = re(cl)
    dac = (clnew - cl) / (1.0 - cl_ac - cl_ms * 2.0 * minf * minf_cl);

    //---- set under-relaxation factor if re change is too large
    if (rlx * dac > dclmax)
      rlx = dclmax / dac;
    if (rlx * dac < dclmin)
      rlx = dclmin / dac;
  } else {
    //===== cl is prescribed: ac is alpha

    //---- set change in alpha to drive cl to prescribed value
    dac = (clnew - clspec) / (0.0 - cl_ac - cl_a);

    //---- set under-relaxation factor if alpha change is too large
    if (rlx * dac > dalmax)
      rlx = dalmax / dac;
    if (rlx * dac < dalmin)
      rlx = dalmin / dac;
  }
  rmsbl = 0.0;
  rmxbl = 0.0;
  dhi = 1.5;
  dlo = -.5;

  SidePair<VectorXd> dctau_seg, dthet_seg, ddstr_seg, duedg_seg;
  for (int is = 1; is <= 2; ++is) {
    // Use 0-based BL indexing for segments (skip ibl=0)
    int start = 0;
    int len = nbl.get(is) - 1;
    if (len <= 0)
      continue;

    dctau_seg.get(is) = VectorXd(len);
    dthet_seg.get(is) = VectorXd(len);
    ddstr_seg.get(is) = VectorXd(len);
    duedg_seg.get(is) = VectorXd(len);

    // Build 0-based local slices for calculations
    VectorXd ctau_slice(len), thet_slice(len), dstr_slice(len), uedg_slice(len);
    VectorXd unew_slice = unew.get(is).segment(0, len);
    VectorXd uac_slice = u_ac.get(is).segment(0, len);
    for (int j = 0; j < len; ++j) {
      ctau_slice[j] = ctau.get(is)[j];
      thet_slice[j] = thet.get(is)[j];
      dstr_slice[j] = dstr.get(is)[j];
      uedg_slice[j] = uedg.get(is)[j];
    }

    VectorXi iv = isys.get(is).segment(start, len);
    VectorXd dmass(len);
    for (int j = 0; j < len; ++j) {
      int idx = iv[j];
      dctau_seg.get(is)[j] = vdel[idx](0, 0) - dac * vdel[idx](0, 1);
      dthet_seg.get(is)[j] = vdel[idx](1, 0) - dac * vdel[idx](1, 1);
      dmass[j] = vdel[idx](2, 0) - dac * vdel[idx](2, 1);
    }
    duedg_seg.get(is) = unew_slice + dac * uac_slice - uedg_slice;
    ddstr_seg.get(is) = (dmass - dstr_slice.cwiseProduct(duedg_seg.get(is)))
                            .cwiseQuotient(uedg_slice);

    VectorXi iblSeq = VectorXi::LinSpaced(len, 0, len - 1);
    // Compare using 0-based logical indices
    VectorXd dn1 =
        ((iblSeq.array()) < itran.get(is))
            .select(dctau_seg.get(is).array() / 10.0,
                    dctau_seg.get(is).cwiseQuotient(ctau_slice).array())
            .matrix();
    VectorXd dn2 = dthet_seg.get(is).cwiseQuotient(thet_slice);
    VectorXd dn3 = ddstr_seg.get(is).cwiseQuotient(dstr_slice);
    VectorXd dn4 = duedg_seg.get(is).array().abs() / 0.25;

    rmsbl += (dn1.array().square() + dn2.array().square() +
              dn3.array().square() + dn4.array().square())
                 .sum();

    auto relax = [&](const VectorXd &dn) {
      double max_pos = dn.cwiseMax(0.0).maxCoeff();
      if (max_pos > 0.0)
        rlx = std::min(rlx, dhi / max_pos);
      double min_neg = dn.cwiseMin(0.0).minCoeff();
      if (min_neg < 0.0)
        rlx = std::min(rlx, dlo / min_neg);
    };
    relax(dn1);
    relax(dn2);
    relax(dn3);
    relax(dn4);

    double local_max = dn1.cwiseAbs().maxCoeff();
    local_max = std::max(local_max, dn2.cwiseAbs().maxCoeff());
    local_max = std::max(local_max, dn3.cwiseAbs().maxCoeff());
    local_max = std::max(local_max, dn4.cwiseAbs().maxCoeff());
    rmxbl = std::max(rmxbl, local_max);
  }

  //--- set true rms change
  rmsbl = sqrt(rmsbl / (4.0 * double(nbl.top + nbl.bottom)));

  if (lalfa) {
    //---- set underrelaxed change in reynolds number from change in lift
    cl = cl + rlx * dac;
  } else {
    //---- set underrelaxed change in alpha
    alfa = alfa + rlx * dac;
  }

  //--- update bl variables with underrelaxed changes
  for (int is = 1; is <= 2; ++is) {
    // 0-based segments starting at ibl=1 (skip 0)
    int start = 1;
    int len = nbl.get(is) - 1;
    if (len <= 0)
      continue;

    // Local 0-based slices, then write back to arrays
    VectorXd ctau_slice(len), thet_slice(len), dstr_slice(len), uedg_slice(len);
    for (int j = 0; j < len; ++j) {
      ctau_slice[j] = ctau.get(is)[j];
      thet_slice[j] = thet.get(is)[j];
      dstr_slice[j] = dstr.get(is)[j];
      uedg_slice[j] = uedg.get(is)[j];
    }

    ctau_slice += rlx * dctau_seg.get(is);
    thet_slice += rlx * dthet_seg.get(is);
    dstr_slice += rlx * ddstr_seg.get(is);
    uedg_slice += rlx * duedg_seg.get(is);

    VectorXi iblSeq = VectorXi::LinSpaced(len, 0, len - 1);
    ctau_slice =
        ((iblSeq.array()) >= itran.get(is))
            .select(ctau_slice.array().cwiseMin(0.25), ctau_slice.array())
            .matrix();

    // Write back updated slices to arrays via 0-based accessors
    for (int j = 0; j < len; ++j) {
      ctau.get(is)[j] = ctau_slice[j];
      thet.get(is)[j] = thet_slice[j];
      dstr.get(is)[j] = dstr_slice[j];
      uedg.get(is)[j] = uedg_slice[j];
    }

    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      if (ibl > iblte.get(is)) {
        dswaki = wgap[ibl - (iblte.get(is) + 1)];
      } else
        dswaki = 0.0;

      if (ibl <= iblte.get(is))
        hklim = 1.02;
      else
        hklim = 1.00005;

      msq = uedg.get(is)[ibl] *
            uedg.get(is)[ibl] * hstinv /
            (gamm1 *
             (1.0 - 0.5 * uedg.get(is)[ibl] *
                        uedg.get(is)[ibl] * hstinv));
      dsw = dstr.get(is)[ibl] - dswaki;
      dslim(dsw, thet.get(is)[ibl], msq, hklim);
      dstr.get(is)[ibl] = dsw + dswaki;

      //------- set new mass defect (nonlinear update)
      mass.get(is)[ibl] =
          dstr.get(is)[ibl] * uedg.get(is)[ibl];
    }
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

bool XFoil::xicalc() {
  //-------------------------------------------------------------
  //     sets bl arc length array on each airfoil side and wake
  //-------------------------------------------------------------

  
    for (int ibl = 0; ibl <= iblte.top; ++ibl) {
      xssi.top[ibl] = sst - spline_length[ipan.get(1)[ibl]];
    }
  
    for (int ibl = 0; ibl <= iblte.bottom; ++ibl) {
      xssi.bottom[ibl] = spline_length[ipan.get(2)[ibl]] - sst;
    }

    // Wake: start from TE, duplicate TE value at first wake station
    xssi.bottom[iblte.bottom + 1] = xssi.bottom[iblte.bottom];
    for (int ibl = iblte.bottom + 2; ibl < nbl.bottom; ++ibl) {
      xssi.bottom[ibl] = xssi.bottom[ibl - 1] +
                          (points.col(ipan.get(2)[ibl]) -
                           points.col(ipan.get(2)[ibl - 1]))
                              .norm();
    }
  

  //---- trailing edge flap length to te gap ratio
  const double telrat = 2.50;

  //---- set up parameters for te flap cubics

  const double crosp = cross2(dpoints_ds.col(n - 1).normalized(),
                              dpoints_ds.col(0).normalized());
  double dwdxte = crosp / sqrt(1.0 - crosp * crosp);

  //---- limit cubic to avoid absurd te gap widths
  dwdxte = std::max(dwdxte, -3.0 / telrat);
  dwdxte = std::min(dwdxte, 3.0 / telrat);

  const double aa = 3.0 + telrat * dwdxte;
  const double bb = -2.0 - telrat * dwdxte;

  if (sharp) {
    for (int iw0 = 0; iw0 < nw; iw0++)
      wgap[iw0] = 0.0;
  }

  else {
    //----- set te flap (wake gap) array (0-based: iw0=0..nw-1)
    for (int iw0 = 0; iw0 < nw; iw0++) {
      const int te_bot_0b = iblte.bottom; // 0-based TE for array indexing
      const double zn = 1.0 - (xssi.bottom[te_bot_0b + (iw0 + 1)] -
                               xssi.bottom[te_bot_0b]) /
                                  (telrat * ante);
      wgap[iw0] = 0.0;
      if (zn >= 0.0)
        wgap[iw0] = ante * (aa + bb * zn) * zn * zn;
    }
  }
  return true;
}

/** -----------------------------------------------------
 * 	   sets forced-transition bl coordinate locations.
 * ----------------------------------------------------- */
double XFoil::xifset(int is) {
  std::stringstream ss;
  VectorXd w1 = VectorXd::Zero(6 * IQX);
  VectorXd w2 = VectorXd::Zero(6 * IQX);
  VectorXd w3 = VectorXd::Zero(6 * IQX);
  VectorXd w4 = VectorXd::Zero(6 * IQX);
  double str;

  if (xstrip.get(is) >= 1.0) {
    return xssi.get(is)[iblte.get(is)];
  }

  Vector2d point_chord = point_te - point_le;

  //---- calculate chord-based x/c, y/c
  for (int i = 0; i < n; i++) {
    w1[i] = (points.col(i) - point_le).dot(point_chord.normalized());
    w2[i] = cross2(points.col(i) - point_le, point_chord.normalized());
  }

  w3 = spline::splind(w1, spline_length.head(n));
  w4 = spline::splind(w2, spline_length.head(n));

  if (is == 1) {
    str = sle + (spline_length[0] - sle) * xstrip.top;
  } else {
    str = sle + (spline_length[n - 1] - sle) * xstrip.bottom;
  }
  str = spline::sinvrt(str, xstrip.get(is), w1, w3, spline_length.head(n), n);
  xiforc = std::min((str - sst), xssi.get(is)[iblte.get(is)]);
  if (xiforc < 0.0) {
    ss << " ***  stagnation point is past trip on side " << is << "\n";
    writeString(ss.str());

    xiforc = xssi.get(is)[iblte.get(is)];
  }

  return xiforc;
}

// moved to XFoil_geometry.cpp: xyWake()

bool XFoil::isValidFoilAngles(Matrix2Xd points) {

  double max_angle = cang(points);
  return max_angle <= angtol;
}

bool XFoil::isValidFoilPointSize(Matrix2Xd points) {
  return points.cols() >= 3;
}
