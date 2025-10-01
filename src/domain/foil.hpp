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
  private:
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
    Eigen::VectorXd setexp(double ds, double chord, int nw) {
    //........................................................
    //     sets geometriy stretched array s:
    //
    //       s(i+1) - s(i)  =  r * [s(i) - s(i-1)]
    //
    //       ds1   (input)   first s increment:  spline_length[2] -
    //       spline_length[1] smax  (input)   final s value:      s(nn) nn (input)
    //       number of points
    //........................................................
    const int nex = nw - 1;

    const double sigma = chord / ds;
    const double rnex = static_cast<double>(nex);
    const double rni = 1.0 / rnex;

    //-- solve quadratic for initial geometric ratio guess
    const double aaa = rnex * (rnex - 1.0) * (rnex - 2.0) / 6.0;
    const double bbb = rnex * (rnex - 1.0) / 2.0;
    const double ccc = rnex - sigma;

    double disc = std::max(0.0, bbb * bbb - 4.0 * aaa * ccc);
    double ratio = 1.0;
    if (nex == 2) {
      ratio = -ccc / bbb + 1.0;
    } else {
      ratio = (-bbb + sqrt(disc)) / (2.0 * aaa) + 1.0;
    }

    //-- newton iteration for actual geometric ratio
    for (int iter = 0; iter < 100; iter++) {
      const double sigman = (pow(ratio, static_cast<double>(nex)) - 1.0) / (ratio - 1.0);
      const double sigman_rni = pow(sigman, rni);
      const double sigma_rni = pow(sigma, rni);
      const double res = sigman_rni - sigma_rni;
      const double numerator = rnex * pow(ratio, static_cast<double>(nex - 1)) - sigman;
      const double denominator = pow(ratio, static_cast<double>(nex)) - 1.0;
      const double dresdr = rni * sigman_rni * numerator / denominator;

      const double dratio = -res / dresdr;
      ratio += dratio;

      if (fabs(dratio) < 1.0e-5) {
        break;
      }
    }

    Eigen::VectorXd spline_length(nw);
    spline_length[0] = 0.0;
    double ds_i = ds;
    for (int i = 1; i < nw; i++) {
      spline_length[i] = spline_length[i - 1] + ds_i;
      ds_i *= ratio;
    }
    return spline_length;
  }
    bool xyWake() {
      //-----------------------------------------------------
      //     sets wake coordinate array for current surface
      //     vorticity and/or mass source distributions.
      //-----------------------------------------------------

      int nw = foil_shape.n / 8 + 2;
      wake_shape.points = Eigen::Matrix2Xd::Zero(2, foil_shape.n + nw);
      wake_shape.n = foil_shape.n + nw;
      wake_shape.normal_vector = Eigen::Matrix2Xd::Zero(2, foil_shape.n + nw);
      wake_shape.spline_length = Eigen::VectorXd::Zero(foil_shape.n + nw);
      wake_shape.points.block(0, 0, 2, foil_shape.n) = foil_shape.points;
      wake_shape.normal_vector.block(0, 0, 2, foil_shape.n) = foil_shape.normal_vector;
      wake_shape.spline_length.head(foil_shape.n) = foil_shape.spline_length;

      double ds1 = 0.5 * (foil_shape.spline_length[1] - foil_shape.spline_length[0] + foil_shape.spline_length[foil_shape.n - 1] - foil_shape.spline_length[foil_shape.n - 2]);
      const auto wake_spacing = setexp(ds1, edge_data.chord, nw);

      //snew.segment(foil_shape.n, nw) = *wake_spacing;
      //-- set first wake point a tiny distance behind te
      Eigen::Vector2d spacing = {
        (foil_shape.dpoints_ds.col(foil_shape.n - 1).y() - foil_shape.dpoints_ds.col(0).y()) / 2,
        (foil_shape.dpoints_ds.col(0).x() - foil_shape.dpoints_ds.col(foil_shape.n - 1).x()) / 2
      };
      wake_shape.normal_vector.col(foil_shape.n) = spacing.normalized();

      wake_shape.points.col(foil_shape.n).x() = edge_data.point_te.x() - 0.0001 * wake_shape.normal_vector.col(foil_shape.n).y();
      wake_shape.points.col(foil_shape.n).y() = edge_data.point_te.y() + 0.0001 * wake_shape.normal_vector.col(foil_shape.n).x();
      wake_shape.spline_length[foil_shape.n] = wake_shape.spline_length[foil_shape.n - 1];
      //---- calculate streamfunction gradient components at first point
      //Vector2d psi = {psilin(foil_shape.points, foil_shape.n, wake_shape.points.col(foil_shape.n), {1.0, 0.0}, false).psi_ni,
      //                psilin(foil_shape.points, foil_shape.n, wake_shape.points.col(foil_shape.n), {0.0, 1.0}, false).psi_ni};
      //---- set unit vector normal to wake at first point
      //wake_shape.normal_vector.col(foil_shape.n + 1) = -psi.normalized();
      //---- set angle of wake panel normal
      //apanel[n] = atan2(psi.y(), psi.x());
      //---- set rest of wake points
      for (int i = foil_shape.n + 1; i < wake_shape.n; i++) {
        const double ds = wake_spacing[i - foil_shape.n] - wake_spacing[i - 1 - foil_shape.n];
        //------ set new point ds downstream of last point
        wake_shape.points.col(i).x() = wake_shape.points.col(i - 1).x() - ds * wake_shape.normal_vector.col(i - 1).y();
        wake_shape.points.col(i).y() = wake_shape.points.col(i - 1).y() + ds * wake_shape.normal_vector.col(i - 1).x();
        wake_shape.spline_length[i] = wake_shape.spline_length[i - 1] + ds;
        if (i != wake_shape.n - 1) {
          //Vector2d psi2 = {psilin(wake_shape.points, i, wake_shape.points.col(i), {1.0, 0.0}, false).psi_ni,
          //                psilin(wake_shape.points, i, wake_shape.points.col(i), {0.0, 1.0}, false).psi_ni};
          //wake_shape.normal_vector.col(i + 1) = -psi2.normalized();
          //wake_shape.apanel[i] = atan2(psi2.y(), psi2.x());
        }
      }
      //---- set wake presence flag and corresponding alpha
      //lwake = true;
      //awake = alfa;
      //---- old source influence matrix is invalid for the new wake geometry
      //lwdij = false;
      return true;
    }
    Foil() = default;
    Foil(Eigen::Matrix2Xd points, int n) {
      foil_shape.setFoilShape(points, n);
      edge_data = getEdgeData(foil_shape);
      xyWake();
    }
};
