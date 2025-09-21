#pragma once

#include <Eigen/Core>

/**
 * @brief Holds the geometric definition of an airfoil.
 */
struct FoilShape {
  Eigen::Matrix2Xd points;  ///< Airfoil and wake points.
  int n = 0;                ///< Number of airfoil points.

  void reset() {
    points.resize(2, 0);
    n = 0;
  }
};

