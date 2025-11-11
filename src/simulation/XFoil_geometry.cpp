#include "XFoil.h"
#include "Eigen/Core"
#include <algorithm>
#include <numbers>

using Eigen::Matrix2Xd;
using Eigen::Vector2d;
using Eigen::VectorXd;

bool XFoil::abcopy(Matrix2Xd copyFrom) {
  const int original_point_count = foil.foil_shape.n;
  int point_count = static_cast<int>(copyFrom.cols());
  if (original_point_count != point_count)
    lblini = false;

  //---- strip out doubled points
  int r = 1;
  while (r < point_count) {
    // FIXME double型の==比較
    if (copyFrom.col(r - 1) == copyFrom.col(r)) {
      for (int j = r; j < point_count - 1; j++) {
        copyFrom.col(j) = copyFrom.col(j + 1);
      }
      point_count -= 1;
    } else {
      r++;
    }
  }
  //--- number of wake points
  nw = point_count / 8 + 2;
  if (nw > IWX) {
    writeString(" XYWake: array size (IWX) too small.\n  Last wake point index reduced.");
    nw = IWX;
  }
  Matrix2Xd foil_points = Matrix2Xd::Zero(2, IZX);
  foil_points.leftCols(point_count) = copyFrom.leftCols(point_count);

  foil.foil_shape.n = point_count;
  initialize();  
  
  foil = Foil(foil_points, point_count);
  updateTrailingEdgeState();
  apanel.head(point_count) = apcalc(foil.foil_shape.points);

  lgamu = false;
  lwake = false;
  lqaij = false;
  ladij = false;
  lwdij = false;
  lipan = false;
  lvconv = false;

  return true;
}

VectorXd XFoil::apcalc(Matrix2Xd points) {
  const int point_count = foil.foil_shape.n;
  VectorXd result = VectorXd::Zero(point_count);
  if (point_count == 0) {
    return result;
  }
  //---- set angles of airfoil panels
  for (int i = 0; i < point_count; i++) {
    Vector2d diff = points.col((i + 1) % point_count) - points.col(i);
    result[i] = atan2(diff.x(), -diff.y());
  }
  //---- TE panel
  if (foil.edge.sharp) {
    result[point_count - 1] = std::numbers::pi;
  }
  return result;
}

double XFoil::atanc(double y, double x, double thold) {
  double tpi, thnew, dthet, dtcorr;
  tpi = 6.2831853071795864769;
  thnew = atan2(y, x);
  dthet = thnew - thold;
  dtcorr = dthet - tpi * int((dthet + sign(std::numbers::pi, dthet)) / tpi);
  return thold + dtcorr;
}

double XFoil::cang(Matrix2Xd points) {
  double max_angle = 0;
  //---- go over each point, calculating corner angle
  for (int i = 1; i < points.cols() - 1; i++) {
    Vector2d delta_former = points.col(i) - points.col(i - 1);
    Vector2d delta_later = points.col(i) - points.col(i + 1);
    double sin = MathUtil::cross2(delta_later, delta_former) / delta_former.norm() / delta_later.norm();
    double delta_angle = asin(sin) * 180.0 / std::numbers::pi;
    max_angle = std::max(fabs(delta_angle), max_angle);
  }
  return max_angle;
}

double XFoil::lefind(const Matrix2Xd &points, const Matrix2Xd &dpoints_ds, const VectorXd &s, int n) const {
  const double dseps = (s[n - 1] - s[0]) * 0.00001;
  const Vector2d point_te_local = 0.5 * (points.col(0) + points.col(n - 1));

  int i = 2;
  for (; i < n - 2; i++) {
    const Vector2d dpoint_te = points.col(i) - point_te_local;
    const Vector2d dpoint = points.col(i + 1) - points.col(i);
    if (dpoint_te.dot(dpoint) < 0.0)
      break;
  }

  const double sle_initial = s[i];
  if (s[i] == s[i - 1])
    return sle_initial;

  double sle_candidate = sle_initial;
  Vector2d point_le_local;
  for (int iter = 1; iter <= 50; iter++) {
    point_le_local.x() = spline::seval(sle_candidate, points.row(0), dpoints_ds.row(0), s, n);
    point_le_local.y() = spline::seval(sle_candidate, points.row(1), dpoints_ds.row(1), s, n);

    Vector2d dpoint_ds_vec;
    dpoint_ds_vec.x() = spline::deval(sle_candidate, points.row(0), dpoints_ds.row(0), s, n);
    dpoint_ds_vec.y() = spline::deval(sle_candidate, points.row(1), dpoints_ds.row(1), s, n);

    Vector2d dpoint_dd_vec;
    dpoint_dd_vec.x() = spline::d2val(sle_candidate, points.row(0), dpoints_ds.row(0), s, n);
    dpoint_dd_vec.y() = spline::d2val(sle_candidate, points.row(1), dpoints_ds.row(1), s, n);

    const Vector2d chord_v = point_le_local - point_te_local;
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

void XFoil::updateTrailingEdgeState() {
  const int node_count = foil.foil_shape.n;

  if (node_count < 2) {
    sigte = 0.0;
    gamte = 0.0;
    return;
  }

  double scs = 0.0;
  double sds = 0.0;
  if (foil.edge.sharp) {
    scs = 1.0;
    sds = 0.0;
  } else if (foil.edge.dste != 0.0) {
    const double inv_dste = 1.0 / foil.edge.dste;
    scs = foil.edge.ante * inv_dste;
    sds = foil.edge.aste * inv_dste;
  }

  if (surface_vortex.rows() > 0 && surface_vortex.cols() >= node_count) {
    const double surface_delta =
        surface_vortex(0, 0) - surface_vortex(0, node_count - 1);
    sigte = 0.5 * surface_delta * scs;
    gamte = -0.5 * surface_delta * sds;
  } else {
    sigte = 0.0;
    gamte = 0.0;
  }
}

bool XFoil::xyWake() {
  if (!foil.xyWake(nw, apanel, gamu, surface_vortex, alfa, qinf)) {
    return false;
  }
  lwake = true;
  lwdij = false;
  return true;
}

double XFoil::sign(double a, double b) {
  if (b >= 0.0)
    return fabs(a);
  else
    return -fabs(a);
}
