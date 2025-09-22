#pragma once

#include <Eigen/Core>

/**
 * @brief Holds the geometric definition of an airfoil.
 */
class FoilShape {
  private:
  Eigen::VectorXd calcSpline() {
    Eigen::VectorXd s = Eigen::VectorXd::Zero(n);
    
    for (int i = 1; i < n; i++) {
      s[i] = s[i - 1] + (points.col(i) - points.col(i - 1)).norm();
    }

    return s;
  }
  public:
  Eigen::Matrix2Xd points;  ///< Airfoil and wake points.
  int n = 0;                ///< Number of airfoil points.
  Eigen::VectorXd spline_length;  ///< Spline lengths for airfoil points.
  void setFoilShape(Eigen::Matrix2Xd points, int n) {
    this->points = points;
    this->n = n;
    this->spline_length = calcSpline();
  }
};

class Foil {
  public:
    FoilShape foil_shape;  // geometric domain (airfoil points and count)
    FoilShape wake_shape;  // geometric domain (wake points and count)
    

};