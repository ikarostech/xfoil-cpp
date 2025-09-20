#include "dissipation.hpp"
#include <cmath>

//---- laminar dissipation function

dissipation::DissipationResult dissipation::dil(double hk, double rt) {
  DissipationResult result;
  if (hk < 4.0) {
    result.di = (0.00205 * pow((4.0 - hk), 5.5) + 0.207) / rt;
    result.di_hk = (-.00205 * 5.5 * pow((4.0 - hk), 4.5)) / rt;
  } else {
    double hkb = hk - 4.0;
    double den = 1.0 + 0.02 * hkb * hkb;
    result.di = (-.0016 * hkb * hkb / den + 0.207) / rt;
    result.di_hk =
        (-.0016 * 2.0 * hkb * (1.0 / den - 0.02 * hkb * hkb / den / den)) / rt;
  }
  result.di_rt = -(result.di) / rt;

  return result;
}

//---- laminar wake dissipation function

dissipation::DissipationResult dissipation::dilw(double hk, double rt) {
  DissipationResult result;
  boundary_layer::ThicknessShapeParameterResult hsl_result = boundary_layer::hsl(hk);
  double rcd = 1.10 * (1.0 - 1.0 / hk) * (1.0 - 1.0 / hk) / hk;
  double rcd_hk = -1.10 * (1.0 - 1.0 / hk) * 2.0 / hk / hk / hk - rcd / hk;

  result.di = 2.0 * rcd / (hsl_result.hs * rt);
  result.di_hk = 2.0 * rcd_hk / (hsl_result.hs * rt) -
                 ((result.di) / hsl_result.hs) * hsl_result.hs_hk;
  result.di_rt = -(result.di) / rt -
                 ((result.di) / hsl_result.hs) * hsl_result.hs_rt;

  return result;
}

dissipation::DissipationResult dissipation::getDissipation(double hk, double rt, FlowRegimeEnum flowRegimeType) {
  if (flowRegimeType == FlowRegimeEnum::Wake) {
    return dilw(hk, rt);
  }
  return dil(hk, rt);
}

