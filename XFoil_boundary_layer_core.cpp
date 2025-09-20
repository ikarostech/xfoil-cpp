#include "XFoil.h"
#include "model/coefficient/skin_friction.hpp"
#include "model/coefficient/dissipation.hpp"
#include <algorithm>
#include <cmath>
using namespace Eigen;


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
  ref = computeShapeParameters(ref, flowRegimeType);
  ref = computeCoefficients(ref, flowRegimeType);
  ref = computeDissipationAndThickness(ref, flowRegimeType);
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
