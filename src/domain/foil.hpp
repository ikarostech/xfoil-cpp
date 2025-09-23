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
  Eigen::Matrix2Xd calcNormalVector() {
  Eigen::Matrix2Xd normal_vector = Eigen::Matrix2Xd::Zero(2, n);
  if (n == 0) {
    return normal_vector;
  }

  const Eigen::VectorXd dx_ds = spline::splind(points.row(0), spline_length);
  const Eigen::VectorXd dy_ds = spline::splind(points.row(1), spline_length);

  for (int i = 0; i < n; i++) {
    const Eigen::Vector2d tangent{dx_ds[i], dy_ds[i]};
    const double tangent_norm = tangent.norm();
    if (tangent_norm == 0.0) {
      continue;
    }
    const Eigen::Vector2d normal{tangent.y(), -tangent.x()};
    normal_vector.col(i) = normal / tangent_norm;
  }
  return normal_vector;
}
  public:
  Eigen::Matrix2Xd points;  ///< Airfoil and wake points.
  int n = 0;                ///< Number of airfoil points.
  Eigen::VectorXd spline_length;  ///< Spline lengths for airfoil points.
  Eigen::Matrix2Xd normal_vector; ///< Normal vectors for airfoil points.
  void setFoilShape(Eigen::Matrix2Xd points, int n) {
    this->points = points;
    this->n = n;
    this->spline_length = calcSpline();
    this->normal_vector = calcNormalVector();
  }
};

class Foil {
  public:
    FoilShape foil_shape;  // geometric domain (airfoil points and count)
    FoilShape wake_shape;  // geometric domain (wake points and count)
    

};
