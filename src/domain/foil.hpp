#pragma once

#include <Eigen/Core>
#include <algorithm>

#include "core/spline.hpp"

/**
 * @brief Holds the geometric definition of an airfoil.
 */
class FoilShape {
 private:
  Eigen::VectorXd calcSpline() const {
    Eigen::VectorXd s = Eigen::VectorXd::Zero(n);

    for (int i = 1; i < n; i++) {
      s[i] = s[i - 1] + (points.col(i) - points.col(i - 1)).norm();
    }

    return s;
  }

  Eigen::Matrix2Xd calcDPointsDs() const {
    Eigen::Matrix2Xd result = Eigen::Matrix2Xd::Zero(2, n);
    if (n < 2) {
      return result;
    }

    const Eigen::VectorXd s_local = spline_length.head(n);
    result.row(0) = spline::splind(points.row(0).head(n).transpose(), s_local);
    result.row(1) = spline::splind(points.row(1).head(n).transpose(), s_local);
    return result;
  }

  Eigen::Matrix2Xd calcNormalVector(const Eigen::Matrix2Xd &dpoints_ds) const {
    Eigen::Matrix2Xd normal_vector = Eigen::Matrix2Xd::Zero(2, n);
    if (n == 0) {
      return normal_vector;
    }

    for (int i = 0; i < n; i++) {
      const Eigen::Vector2d tangent = dpoints_ds.col(i);
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
  Eigen::Matrix2Xd dpoints_ds;    ///< First derivative of airfoil coordinates with respect to arc length.
  Eigen::Matrix2Xd normal_vector; ///< Normal vectors for airfoil points.
  Eigen::Vector2d point_le = Eigen::Vector2d::Zero(); ///< Leading edge point.
  Eigen::Vector2d point_te = Eigen::Vector2d::Zero(); ///< Trailing edge point.
  double chord = 0.0;              ///< Chord length.

  void setFoilShape(const Eigen::Matrix2Xd &points, int n) {
    this->points = points;
    this->n = n;
    this->spline_length = calcSpline();
    this->dpoints_ds = calcDPointsDs();
    this->normal_vector = calcNormalVector(this->dpoints_ds);
  }

  void updateEdges(double sle) {
    if (n < 2) {
      point_le = Eigen::Vector2d::Zero();
      point_te = Eigen::Vector2d::Zero();
      chord = 0.0;
      return;
    }

    const int cols = n;
    const Eigen::VectorXd s_local = spline_length.head(cols);
    const double sle_clamped = std::min(std::max(sle, s_local[0]), s_local[cols - 1]);

    const Eigen::VectorXd x = points.row(0).head(cols).transpose();
    const Eigen::VectorXd y = points.row(1).head(cols).transpose();
    const Eigen::VectorXd xs = dpoints_ds.row(0).head(cols).transpose();
    const Eigen::VectorXd ys = dpoints_ds.row(1).head(cols).transpose();

    point_le.x() = spline::seval(sle_clamped, x, xs, s_local, cols);
    point_le.y() = spline::seval(sle_clamped, y, ys, s_local, cols);
    point_te = 0.5 * (points.col(0) + points.col(cols - 1));
    chord = (point_le - point_te).norm();
  }
};



class Foil {
  public:
    FoilShape foil_shape;  // geometric domain (airfoil points and count)
    FoilShape wake_shape;  // geometric domain (wake points and count)
    

};
