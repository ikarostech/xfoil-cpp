#include "boundary_layer_util.hpp"
#include <cmath>
BoundaryLayerUtil::AxResult BoundaryLayerUtil::axset(double hk1, double t1, double rt1, double a1,
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
BoundaryLayerUtil::EnvEnResult BoundaryLayerUtil::dampl(double hk, double th, double rt) {
  BoundaryLayerUtil::EnvEnResult result;
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