#include "domain/foil/edge.hpp"

#include "core/spline.hpp"

#include <algorithm>
#include <cmath>

namespace {
// 2D cross product (z-component)
inline double cross2(const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
  return a[0] * b[1] - a[1] * b[0];
}
} // namespace

Edge::Edge()
    : point_le(Eigen::Vector2d::Zero()),
      point_te(Eigen::Vector2d::Zero()),
      chord(0.0),
      sle(0.0),
      ante(0.0),
      aste(0.0),
      dste(0.0),
      sharp(false) {}

Edge::Edge(const FoilShape& foilShape) : Edge() {

  sle = lefind(foilShape.points,
               foilShape.dpoints_ds,
               foilShape.spline_length.head(foilShape.n),
               foilShape.n);

  point_le.x() = spline::seval(sle,
                               foilShape.points.row(0),
                               foilShape.dpoints_ds.row(0),
                               foilShape.spline_length.head(foilShape.n),
                               foilShape.n);
  point_le.y() = spline::seval(sle,
                               foilShape.points.row(1),
                               foilShape.dpoints_ds.row(1),
                               foilShape.spline_length.head(foilShape.n),
                               foilShape.n);
  point_te = 0.5 * (foilShape.points.col(0) + foilShape.points.col(foilShape.n - 1));
  chord = (point_le - point_te).norm();

  const Eigen::Vector2d tevec =
      foilShape.points.col(0) - foilShape.points.col(foilShape.n - 1);
  const Eigen::Vector2d dpoint_ds_te =
      0.5 * (-foilShape.dpoints_ds.col(0) +
             foilShape.dpoints_ds.col(foilShape.n - 1));

  ante = cross2(dpoint_ds_te, tevec);
  aste = tevec.dot(dpoint_ds_te);
  dste = tevec.norm();

  sharp = dste < 0.0001 * chord;
}

double Edge::lefind(const Eigen::Matrix2Xd& points,
                    const Eigen::Matrix2Xd& dpoints_ds,
                    const Eigen::VectorXd& s,
                    int n) {
    const double dseps = (s[n - 1] - s[0]) * 0.00001;
    const Eigen::Vector2d point_te_local = 0.5 * (points.col(0) + points.col(n - 1));

    int i = 2;
    for (; i < n - 2; i++) {
    const Eigen::Vector2d dpoint_te = points.col(i) - point_te_local;
    const Eigen::Vector2d dpoint = points.col(i + 1) - points.col(i);
    if (dpoint_te.dot(dpoint) < 0.0) {
      break;
    }
    }

    const double sle_initial = s[i];
    if (s[i] == s[i - 1]) {
      return sle_initial;
    }

    double sle_candidate = sle_initial;
    Eigen::Vector2d point_le_local;
    for (int iter = 1; iter <= 50; iter++) {
      point_le_local.x() =
          spline::seval(sle_candidate, points.row(0), dpoints_ds.row(0), s, n);
      point_le_local.y() =
          spline::seval(sle_candidate, points.row(1), dpoints_ds.row(1), s, n);

      Eigen::Vector2d dpoint_ds_vec;
      dpoint_ds_vec.x() =
          spline::deval(sle_candidate, points.row(0), dpoints_ds.row(0), s, n);
      dpoint_ds_vec.y() =
          spline::deval(sle_candidate, points.row(1), dpoints_ds.row(1), s, n);

      Eigen::Vector2d dpoint_dd_vec;
      dpoint_dd_vec.x() =
          spline::d2val(sle_candidate, points.row(0), dpoints_ds.row(0), s, n);
      dpoint_dd_vec.y() =
          spline::d2val(sle_candidate, points.row(1), dpoints_ds.row(1), s, n);

      const Eigen::Vector2d chord_v = point_le_local - point_te_local;
      const double res = chord_v.dot(dpoint_ds_vec);
      const double ress =
          dpoint_ds_vec.dot(dpoint_ds_vec) + chord_v.dot(dpoint_dd_vec);
      double dsle = -res / ress;
      const double chord_sum = std::fabs(chord_v.x() + chord_v.y());
      dsle = std::max(dsle, -0.02 * chord_sum);
      dsle = std::min(dsle, 0.02 * chord_sum);
      sle_candidate += dsle;
      if (std::fabs(dsle) < dseps) {
        return sle_candidate;
      }
    }

    return sle_initial;
}
