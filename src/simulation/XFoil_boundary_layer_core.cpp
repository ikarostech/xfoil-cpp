#include "XFoil.h"
#include "domain/coefficient/skin_friction.hpp"
#include "domain/coefficient/dissipation.hpp"
#include <algorithm>
#include <cmath>
using namespace Eigen;

/** ----------------------------------------------------------
 *      returns average amplification ax over interval 1..2
 * ----------------------------------------------------------- */
XFoil::AxResult XFoil::axset(double hk1, double t1, double rt1, double a1,
                             double hk2, double t2, double rt2, double a2,
                             double acrit) const {
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

bool XFoil::blkin(BoundaryLayerState& state) {
  //----------------------------------------------------------
  //     calculates turbulence-independent secondary "2"
  //     variables from the primary "2" variables.
  //----------------------------------------------------------
  double tr2, herat, he_u2, he_ms, v2_he;
  blData& current = state.current();
  //---- set edge mach number ** 2
  current.param.mz =
      current.param.uz * current.param.uz * hstinv /
      (gm1bl * (1.0 - 0.5 * current.param.uz * current.param.uz * hstinv));
  tr2 = 1.0 + 0.5 * gm1bl * current.param.mz;
  current.param.mz_uz = 2.0 * current.param.mz * tr2 / current.param.uz;
  current.param.mz_ms =
      current.param.uz * current.param.uz * tr2 /
      (gm1bl * (1.0 - 0.5 * current.param.uz * current.param.uz * hstinv)) *
      hstinv_ms;

  //---- set edge density (isentropic relation)
  current.param.rz = rstbl * pow(tr2, (-1.0 / gm1bl));
  current.param.rz_uz = -current.param.rz / tr2 * 0.5 * current.param.mz_uz;
  current.param.rz_ms = -current.param.rz / tr2 * 0.5 * current.param.mz_ms +
                        rstbl_ms * pow(tr2, (-1.0 / gm1bl));

  //---- set shape parameter
  current.param.hz = current.param.dz / current.param.tz;
  current.param.hz_dz = 1.0 / current.param.tz;
  current.param.hz_tz = -current.param.hz / current.param.tz;

  //---- set edge static/stagnation enthalpy
  herat = 1.0 - 0.5 * current.param.uz * current.param.uz * hstinv;
  he_u2 = -current.param.uz * hstinv;
  he_ms = -0.5 * current.param.uz * current.param.uz * hstinv_ms;
  //---- set molecular viscosity
  v2_he = (1.5 / herat - 1.0 / (herat + hvrat));

  //---- set kinematic shape parameter
  boundary_layer::KineticShapeParameterResult hkin_result =
      boundary_layer::hkin(current.param.hz, current.param.mz);
  current.hkz.scalar = hkin_result.hk;

  current.hkz.u() = hkin_result.hk_msq * current.param.mz_uz;
  current.hkz.t() = hkin_result.hk_h * current.param.hz_tz;
  current.hkz.d() = hkin_result.hk_h * current.param.hz_dz;
  current.hkz.ms() = hkin_result.hk_msq * current.param.mz_ms;

  //---- set momentum thickness reynolds number
  current.rtz.scalar =
      current.param.rz * current.param.uz * current.param.tz /
      (sqrt(herat * herat * herat) * (1.0 + hvrat) / (herat + hvrat) / reybl);
  current.rtz.u() = current.rtz.scalar *
                    (1.0 / current.param.uz +
                     current.param.rz_uz / current.param.rz - v2_he * he_u2);
  current.rtz.t() = current.rtz.scalar / current.param.tz;
  current.rtz.ms() =
      current.rtz.scalar * (current.param.rz_ms / current.param.rz +
                            (1 / reybl * reybl_ms - v2_he * he_ms));
  current.rtz.re() = current.rtz.scalar * (reybl_re / reybl);

  return true;
}


bool XFoil::blmid(BoundaryLayerState& state, FlowRegimeEnum flowRegimeType) {

  blData& previous = state.previous();
  blData& current = state.current();

  //---- set similarity variables if not defined
  if (simi) {
    previous.hkz = current.hkz;
    previous.rtz = current.rtz;
    previous.param.mz = current.param.mz;
    previous.param.mz_uz = current.param.mz_uz;
    previous.param.mz_ms = current.param.mz_ms;
  }

  //---- define stuff for midpoint cf
  double hka = 0.5 * (previous.hkz.scalar + current.hkz.scalar);
  double rta = 0.5 * (previous.rtz.scalar + current.rtz.scalar);
  double ma = 0.5 * (previous.param.mz + current.param.mz);

  //---- compute midpoint skin friction coefficient
  skin_friction::C_f cf_res = skin_friction::getSkinFriction(hka, rta, ma, flowRegimeType);

  cfm = cf_res.cf;
  double cfm_hka = cf_res.hk;
  double cfm_rta = cf_res.rt;
  double cfm_ma = cf_res.msq;

  cfm_u1 = 0.5 * (cfm_hka * previous.hkz.u() + cfm_ma * previous.param.mz_uz +
                  cfm_rta * previous.rtz.u());
  cfm_t1 = 0.5 * (cfm_hka * previous.hkz.t() + cfm_rta * previous.rtz.t());
  cfm_d1 = 0.5 * (cfm_hka * previous.hkz.d());

  cfm_u2 = 0.5 * (cfm_hka * current.hkz.u() + cfm_ma * current.param.mz_uz +
                  cfm_rta * current.rtz.u());
  cfm_t2 = 0.5 * (cfm_hka * current.hkz.t() + cfm_rta * current.rtz.t());
  cfm_d2 = 0.5 * (cfm_hka * current.hkz.d());

  cfm_ms = 0.5 * (cfm_hka * previous.hkz.ms() + cfm_ma * previous.param.mz_ms +
                  cfm_rta * previous.rtz.ms() + cfm_hka * current.hkz.ms() +
                  cfm_ma * current.param.mz_ms + cfm_rta * current.rtz.ms());
  cfm_re = 0.5 * (cfm_rta * previous.rtz.re() + cfm_rta * current.rtz.re());

  return true;
}

bool XFoil::blmid(FlowRegimeEnum flowRegimeType) {
  return blmid(boundaryLayerState, flowRegimeType);
}


/** ----------------------------------------------------------
 *     set bl primary "2" variables from parameter list
 *  ---------------------------------------------------------- */
bool XFoil::blprv(BoundaryLayerState& state, double xsi, double ami, double cti,
                  double thi, double dsi, double dswaki, double uei) {
  blData& current = state.current();
  current.param.xz = xsi;
  current.param.amplz = ami;
  current.param.sz = cti;
  current.param.tz = thi;
  current.param.dz = dsi - dswaki;
  current.param.dwz = dswaki;

  current.param.uz =
      uei * (1.0 - tkbl) / (1.0 - tkbl * (uei / qinfbl) * (uei / qinfbl));
  current.param.uz_uei =
      (1.0 + tkbl * (2.0 * current.param.uz * uei / qinfbl / qinfbl - 1.0)) /
      (1.0 - tkbl * (uei / qinfbl) * (uei / qinfbl));
  current.param.uz_ms =
      (current.param.uz * (uei / qinfbl) * (uei / qinfbl) - uei) * tkbl_ms /
      (1.0 - tkbl * (uei / qinfbl) * (uei / qinfbl));
  return true;
}

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
bool XFoil::blsys(BoundaryLayerState& state, [[maybe_unused]] BoundaryLayerLattice& lattice) {
  blData& previous = state.previous();
  blData& current = state.current();

  //---- calculate secondary bl variables and their sensitivities
  if (wake) {
    blvar(current, FlowRegimeEnum::Wake);
    blmid(state, FlowRegimeEnum::Wake);
  } else {
    if (turb || tran) {
      blvar(current, FlowRegimeEnum::Turbulent);
      blmid(state, FlowRegimeEnum::Turbulent);
    } else {
      blvar(current, FlowRegimeEnum::Laminar);
      blmid(state, FlowRegimeEnum::Laminar);
    }
  }

  //---- for the similarity station, "1" and "2" variables are the same
  if (simi) {
    //		for(int icom=1;icom<= ncom;icom++) com1[icom] = com2[icom];
    stepbl(state);
  }

  //---- set up appropriate finite difference system for current interval
  if (tran)
    trdif();
  else if (simi)
    blc = bldif(0, boundaryLayerState);
  else if (!turb)
    blc = bldif(1, boundaryLayerState);
  else if (wake)
    blc = bldif(3, boundaryLayerState);
  else
    blc = bldif(2, boundaryLayerState);

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
    blc.a1(k, 3) *= previous.param.uz_uei;
    blc.a2(k, 3) *= current.param.uz_uei;
    blc.d_msq[k] =
        res_u1 * previous.param.uz_ms + res_u2 * current.param.uz_ms + res_ms;
  }
  return true;
}

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
XFoil::EnvEnResult XFoil::dampl(double hk, double th, double rt) const {
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
