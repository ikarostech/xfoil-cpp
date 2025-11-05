#include "XFoil.h"

#include <algorithm>
bool XFoil::dslim(double &dstr, double thet, double msq, double hklim) {
  const double h = (dstr) / thet;

  boundary_layer::KineticShapeParameterResult hkin_result =
      boundary_layer::hkin(h, msq);

  const double dh = std::max(0.0, hklim - hkin_result.hk) / hkin_result.hk_h;
  dstr = (dstr) + dh * thet;

  return true;
}
