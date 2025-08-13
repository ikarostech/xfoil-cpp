#include "skin_friction.hpp"
#include <numeric>
#include <algorithm>
#include <cmath>
using namespace std;
skin_friction::C_f skin_friction::cfl(double hk, double rt) {
  skin_friction::C_f c_f = skin_friction::C_f();
  if (hk < 5.5) {
    double tmp = pow(5.5 - hk, 3) / (hk + 1.0);
    c_f.cf = (0.0727 * tmp - 0.07) / rt;
    c_f.hk =
        (-0.0727 * tmp * 3.0 / (5.5 - hk) - 0.0727 * tmp / (hk + 1.0)) / rt;
  } else {
    double tmp = 1.0 - 1.0 / (hk - 4.5);
    c_f.cf = (0.015 * tmp * tmp - 0.07) / rt;
    c_f.hk = 0.015 * tmp * 2.0 / pow(hk - 4.5, 2.0) / rt;
  }
  c_f.rt = -c_f.cf / rt;
  c_f.msq = 0.0;
  return c_f;
}

skin_friction::C_f skin_friction::cft(double hk, double rt, double msq) {
  C_f c_f = C_f();

  //---- turbulent skin friction function  ( cf )    (coles)
  double gm1 = 1.4 - 1.0;
  double fc = sqrt(1.0 + 0.5 * gm1 * msq);
  double grt = std::max(log(rt / fc), 3.0);

  double gex = -1.74 - 0.31 * hk;

  double f_arg = std::max(-1.33 * hk, -20.0);

  double tanh_hk = tanh(4.0 - hk / 0.875);

  double cfo = 0.3 * exp(f_arg) * pow((grt / 2.3026), gex);
  c_f.cf = (cfo + 0.00011 * (tanh_hk - 1.0)) / fc;
  c_f.hk = (-1.33 * cfo - 0.31 * log(grt / 2.3026) * cfo -
            0.00011 * (1.0 - pow(tanh_hk, 2)) / 0.875) /
           fc;
  c_f.rt = gex * cfo / (fc * grt) / rt;
  c_f.msq = gex * cfo / (fc * grt) * (-0.25 * gm1 / fc / fc) -
            0.25 * gm1 * (c_f.cf) / pow(fc, 2);

  return c_f;
}

skin_friction::C_f skin_friction::getSkinFriction(double hk, double rt, double msq, FlowRegimeEnum flowRegimeType) {
    if (flowRegimeType == FlowRegimeEnum::Wake) {
        return C_f(); // Wake flow has no skin friction
    }
    else if (flowRegimeType == FlowRegimeEnum::Laminar) {
        return skin_friction::cfl(hk, rt);
    } else {
        C_f cfl = skin_friction::cfl(hk, rt);
        C_f cft = skin_friction::cft(hk, rt, msq);
        return cfl.cf > cft.cf ? cfl : cft;
    }
}