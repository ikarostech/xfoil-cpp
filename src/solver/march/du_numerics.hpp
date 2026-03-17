#pragma once

#include <algorithm>

namespace march::du_numerics
{
struct SensitivityUpdateMode
{
  bool reset = false;
  bool average = false;
};

inline double computeInitialDisplacement(double displacementThickness,
                                         double wakeGap,
                                         double momentumThickness,
                                         int stationIndex,
                                         int trailingEdgeIndex)
{
  const double thicknessLimit =
      (stationIndex <= trailingEdgeIndex) ? 1.02 : 1.00005;
  return std::max(displacementThickness - wakeGap,
                  thicknessLimit * momentumThickness) +
         wakeGap;
}

inline SensitivityUpdateMode computeSensitivityUpdateMode(int iteration)
{
  SensitivityUpdateMode mode;
  mode.reset = (iteration <= 5);
  mode.average = (iteration > 5 && iteration <= 15);
  return mode;
}
} // namespace march::du_numerics
