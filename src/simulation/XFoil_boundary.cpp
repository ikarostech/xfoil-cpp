#include "XFoil.h"

#include <algorithm>
bool XFoil::dslim(double &displacementThickness, double momentumThickness, double msq, double hklim) {
  const double h = (displacementThickness) / momentumThickness;

  boundary_layer::KineticShapeParameterResult hkin_result =
      boundary_layer::hkin(h, msq);

  const double dh = std::max(0.0, hklim - hkin_result.hk) / hkin_result.hk_h;
  displacementThickness = (displacementThickness) + dh * momentumThickness;

  return true;
}
