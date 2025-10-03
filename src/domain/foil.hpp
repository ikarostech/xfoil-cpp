#pragma once

#include <Eigen/Core>

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
  void setFoilShape(Eigen::Matrix2Xd points, int n) {
    this->points = points;
    this->n = n;
    this->spline_length = calcSpline();
    this->dpoints_ds = calcDPointsDs();
    this->normal_vector = calcNormalVector(this->dpoints_ds);
  }
};

class Foil {
  public:
    FoilShape foil_shape;  // geometric domain (airfoil points and count)
    FoilShape wake_shape;  // geometric domain (wake points and count)
    class EdgeData {
      public:
      Eigen::Vector2d point_le;  // leading edge point
      Eigen::Vector2d point_te;  // trailing edge point
      double chord;       // chord length
    };
    EdgeData edge_data;

    double lefind(const Eigen::Matrix2Xd &points, const Eigen::Matrix2Xd &dpoints_ds, const Eigen::VectorXd &s, int n) const {
      const double dseps = (s[n - 1] - s[0]) * 0.00001;
      const Eigen::Vector2d point_te_local = 0.5 * (points.col(0) + points.col(n - 1));

      int i = 2;
      for (; i < n - 2; i++) {
        const Eigen::Vector2d dpoint_te = points.col(i) - point_te_local;
        const Eigen::Vector2d dpoint = points.col(i + 1) - points.col(i);
        if (dpoint_te.dot(dpoint) < 0.0)
          break;
      }

      const double sle_initial = s[i];
      if (s[i] == s[i - 1])
        return sle_initial;

      double sle_candidate = sle_initial;
      Eigen::Vector2d point_le_local;
      for (int iter = 1; iter <= 50; iter++) {
        point_le_local.x() = spline::seval(sle_candidate, points.row(0), dpoints_ds.row(0), s, n);
        point_le_local.y() = spline::seval(sle_candidate, points.row(1), dpoints_ds.row(1), s, n);

        Eigen::Vector2d dpoint_ds_vec;
        dpoint_ds_vec.x() = spline::deval(sle_candidate, points.row(0), dpoints_ds.row(0), s, n);
        dpoint_ds_vec.y() = spline::deval(sle_candidate, points.row(1), dpoints_ds.row(1), s, n);

        Eigen::Vector2d dpoint_dd_vec;
        dpoint_dd_vec.x() = spline::d2val(sle_candidate, points.row(0), dpoints_ds.row(0), s, n);
        dpoint_dd_vec.y() = spline::d2val(sle_candidate, points.row(1), dpoints_ds.row(1), s, n);

        const Eigen::Vector2d chord_v = point_le_local - point_te_local;
        const double res = chord_v.dot(dpoint_ds_vec);
        const double ress = dpoint_ds_vec.dot(dpoint_ds_vec) + chord_v.dot(dpoint_dd_vec);
        double dsle = -res / ress;
        const double chord_sum = fabs(chord_v.x() + chord_v.y());
        dsle = std::max(dsle, -0.02 * chord_sum);
        dsle = std::min(dsle, 0.02 * chord_sum);
        sle_candidate += dsle;
        if (fabs(dsle) < dseps)
          return sle_candidate;
      }

      return sle_initial;
    }

  
    EdgeData getEdgeData(FoilShape foilShape) {
      EdgeData edgeData;
      double sle = lefind(foilShape.points, foilShape.dpoints_ds, foilShape.spline_length, foilShape.n);
      edgeData.point_le.x() = spline::seval(sle, foil_shape.points.row(0), foil_shape.dpoints_ds.row(0), foil_shape.spline_length.head(foil_shape.n), foil_shape.n);
      edgeData.point_le.y() = spline::seval(sle, foil_shape.points.row(1), foil_shape.dpoints_ds.row(1), foil_shape.spline_length.head(foil_shape.n), foil_shape.n);
      edgeData.point_te = 0.5 * (foil_shape.points.col(0) + foil_shape.points.col(foil_shape.n - 1));
      edgeData.chord = (edgeData.point_le - edgeData.point_te).norm();
      return edgeData;
    }
  
    Foil() = default;
    Foil(Eigen::Matrix2Xd points, int n) {
      foil_shape.setFoilShape(points, n);
      edge_data = getEdgeData(foil_shape);
    }
};
