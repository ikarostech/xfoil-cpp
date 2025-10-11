#pragma once

#include <Eigen/Core>

#include "core/spline.hpp"
#include "domain/foil/edge.hpp"

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
    if (n == 0) {
      return result;
    }

    result.row(0) = spline::splind(points.row(0), spline_length);
    result.row(1) = spline::splind(points.row(1), spline_length);
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

  FoilShape() = default;
  FoilShape(const Eigen::Matrix2Xd& points, int n) { setFoilShape(points, n); }

  void setFoilShape(const Eigen::Matrix2Xd& points, int n) {
    this->points = points;
    this->n = n;
    this->spline_length = calcSpline();
    this->dpoints_ds = calcDPointsDs();
    this->normal_vector = calcNormalVector(this->dpoints_ds);
  }
};
class Foil {
  private:
  public:
    FoilShape foil_shape;  // geometric domain (airfoil points and count)
    FoilShape wake_shape;  // geometric domain (wake points and count)

    Edge edge;

    Foil() = default;
    Foil(const Eigen::Matrix2Xd& points, int n) {
      foil_shape = FoilShape(points, n);
      edge = Edge(foil_shape);
      // xyWake();
    }
};
