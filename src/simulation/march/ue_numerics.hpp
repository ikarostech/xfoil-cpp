#pragma once

#include <algorithm>
#include <cmath>

#include "core/math_util.hpp"

namespace march::ue_numerics
{
inline double adjustDisplacementForHkLimit(double displacementThickness,
                                           double momentumThickness,
                                           double msq, double hklim)
{
  const double h = displacementThickness / momentumThickness;
  const auto hkinResult = boundary_layer::hkin(h, msq);
  const double dh = std::max(0.0, hklim - hkinResult.hk) / hkinResult.hk_h;
  return displacementThickness + dh * momentumThickness;
}

inline double computeMachSquared(double edgeVelocity, double hstinv,
                                 double gm1bl)
{
  const double ueiSq = edgeVelocity * edgeVelocity;
  return ueiSq * hstinv / (gm1bl * (1.0 - 0.5 * ueiSq * hstinv));
}

inline double computeProposedHk(double displacementThickness,
                                double momentumThickness, double relaxation,
                                double displacementDelta,
                                double momentumDelta, double msq)
{
  const double htest = (displacementThickness + relaxation * displacementDelta) /
                       (momentumThickness + relaxation * momentumDelta);
  return boundary_layer::hkin(htest, msq).hk;
}

inline double computeDmax(double rhs0, double rhs1, double rhs2, double thi,
                          double dsi, double cti, bool includeLaminarAmpTerm,
                          bool postTransition)
{
  double dmaxLocal = std::max(std::fabs(rhs1 / thi), std::fabs(rhs2 / dsi));
  if (includeLaminarAmpTerm && !postTransition)
  {
    dmaxLocal = std::max(dmaxLocal, std::fabs(rhs0 / 10.0));
  }
  if (postTransition)
  {
    dmaxLocal = std::max(dmaxLocal, std::fabs(rhs0 / cti));
  }
  return dmaxLocal;
}

inline double computeRelaxation(double dmaxLocal)
{
  if (dmaxLocal > 0.3)
  {
    return 0.3 / dmaxLocal;
  }
  return 1.0;
}

inline double computeInverseTarget(double currentHk, double targetHk,
                                   double minimumHk)
{
  (void)currentHk;
  (void)targetHk;
  return minimumHk;
}
} // namespace march::ue_numerics
