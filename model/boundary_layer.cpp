#include "boundary_layer.hpp"

#include "math_util.hpp"
#include <cmath>

bool boundary_layer::hct(double hk, double msq, double &hc, double &hc_hk, double &hc_msq) {
    //---- density shape parameter    (from whitfield)
    hc = msq * (0.064 / (hk - 0.8) + 0.251);
    hc_hk = msq * (-.064 / (hk - 0.8) / (hk - 0.8));
    hc_msq = 0.064 / (hk - 0.8) + 0.251;

    return true;
}

bool boundary_layer::hkin(double h, double msq, double &hk, double &hk_h, double &hk_msq) {
  //---- calculate kinematic shape parameter (assuming air)
  //     (from Whitfield )
  hk = (h - 0.29 * msq) / (1.0 + 0.113 * msq);
  hk_h = 1.0 / (1.0 + 0.113 * msq);
  hk_msq = (-.29 - 0.113 * (hk)) / (1.0 + 0.113 * msq);

  return true;
}

bool boundary_layer::hsl(double hk, double &hs, double &hs_hk, double &hs_rt,
                double &hs_msq) {
  //---- laminar hs correlation
  if (hk < 4.35) {
    double tmp = hk - 4.35;
    hs = 0.0111 * tmp * tmp / (hk + 1.0) -
         0.0278 * tmp * tmp * tmp / (hk + 1.0) + 1.528 -
         0.0002 * (tmp * hk) * (tmp * hk);
    hs_hk =
        0.0111 * (2.0 * tmp - tmp * tmp / (hk + 1.0)) / (hk + 1.0) -
        0.0278 * (3.0 * tmp * tmp - tmp * tmp * tmp / (hk + 1.0)) / (hk + 1.0) -
        0.0002 * 2.0 * tmp * hk * (tmp + hk);
  } else {
    hs = 0.015 * (hk - 4.35) * (hk - 4.35) / hk + 1.528;
    hs_hk = 0.015 * 2.0 * (hk - 4.35) / hk -
            0.015 * (hk - 4.35) * (hk - 4.35) / hk / hk;
  }

  hs_rt = 0.0;
  hs_msq = 0.0;

  return true;
}

bool boundary_layer::hst(double hk, double rt, double msq, double &hs, double &hs_hk,
                double &hs_rt, double &hs_msq) {

    //---- turbulent hs correlation
    const double hsmin = 1.5;
    const double dhsinf = 0.015;

    //---- ###  12/4/94
    //---- limited rtheta dependence for rtheta < 200

    double ho, ho_rt;
    if (rt > 400.0) {
        ho = 3.0 + 400.0 / rt;
        ho_rt = -400.0 / MathUtil::pow(rt, 2);
    } else {
        ho = 4.0;
        ho_rt = 0.0;
    }

    double rtz, rtz_rt;
    if (rt > 200.0) {
        rtz = rt;
        rtz_rt = 1.0;
    } else {
        rtz = 200.0;
        rtz_rt = 0.0;
    }

    if (hk < ho) {
        //----- attached branch
        //----- new correlation  29 nov 91
        //-     (from  arctan(y+) + schlichting  profiles)
        const double hr = (ho - hk) / (ho - 1.0);
        const double hr_hk = -1.0 / (ho - 1.0);
        const double hr_rt = (1.0 - hr) / (ho - 1.0) * ho_rt;
        hs = (2.0 - hsmin - 4.0 / rtz) * hr * hr * 1.5 / (hk + 0.5) + hsmin +
            4.0 / rtz;
        hs_hk =
            -(2.0 - hsmin - 4.0 / rtz) * MathUtil::pow(hr, 2) * 1.5 / (hk + 0.5) / (hk + 0.5) +
            (2.0 - hsmin - 4.0 / rtz) * hr * 2.0 * 1.5 / (hk + 0.5) * hr_hk;
        hs_rt = (2.0 - hsmin - 4.0 / rtz) * hr * 2.0 * 1.5 / (hk + 0.5) * hr_rt +
                (hr * hr * 1.5 / (hk + 0.5) - 1.0) * 4.0 / rtz / rtz * rtz_rt;
    } else {
        //----- separated branch
        const double grt = log(rtz);
        const double hdif = hk - ho;
        const double rtmp = hk - ho + 4.0 / grt;
        const double htmp = 0.007 * grt / rtmp / rtmp + dhsinf / hk;
        const double htmp_hk = -.014 * grt / rtmp / rtmp / rtmp - dhsinf / hk / hk;
        const double htmp_rt = -.014 * grt / rtmp / rtmp / rtmp *
                    (-ho_rt - 4.0 / grt / grt / rtz * rtz_rt) +
                0.007 / rtmp / rtmp / rtz * rtz_rt;
        hs = hdif * hdif * htmp + hsmin + 4.0 / rtz;
        hs_hk = hdif * 2.0 * htmp + hdif * hdif * htmp_hk;
        hs_rt = hdif * hdif * htmp_rt - 4.0 / rtz / rtz * rtz_rt +
                hdif * 2.0 * htmp * (-ho_rt);
    }

    //---- whitfield's minor additional compressibility correction
    double fm = 1.0 + 0.014 * msq;
    hs = (hs + 0.028 * msq) / fm;
    hs_hk = (hs_hk) / fm;
    hs_rt = (hs_rt) / fm;
    hs_msq = 0.028 / fm - 0.014 * (hs) / fm;

    return true;
}