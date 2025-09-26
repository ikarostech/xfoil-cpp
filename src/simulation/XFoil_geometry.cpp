#include "XFoil.h"
#include "Eigen/Core"
#include <algorithm>
#include <numbers>

using namespace Eigen;

namespace {
// 2D cross product (z-component)
inline double cross2(const Eigen::Vector2d &a, const Eigen::Vector2d &b) {
  return a[0] * b[1] - a[1] * b[0];
}
}

bool XFoil::abcopy(Matrix2Xd copyFrom) {
  if (n != copyFrom.cols())
    lblini = false;

  n = copyFrom.cols();

  //---- strip out doubled points
  int r = 1;
  while (r < n) {
    // FIXME double型の==比較
    if (copyFrom.col(r - 1) == copyFrom.col(r)) {
      for (int j = r; j < n - 1; j++) {
        copyFrom.col(j) = copyFrom.col(j + 1);
      }
      n = n - 1;
    } else {
      r++;
    }
  }
  //--- number of wake points
  nw = n / 8 + 2;
  if (nw > IWX) {
    writeString(" XYWake: array size (IWX) too small.\n  Last wake point index reduced.");
    nw = IWX;
  }
  points = Matrix2Xd::Zero(2, IZX);
  for (int i = 0; i < n; i++) {
    points.col(i) = copyFrom.col(i);
  }

  initialize();
  //TODO foil側に寄せて最終的にfoilを返すようにする
  
  
  foil.foil_shape.setFoilShape(points, n);
  spline_length.head(n) = foil.foil_shape.spline_length;
  normal_vectors.block(0, 0, 2, n) = foil.foil_shape.normal_vector;
  sle = lefind(points, dpoints_ds, spline_length, n);
  point_le.x() = spline::seval(sle, points.row(0), dpoints_ds.row(0), spline_length.head(n), n);
  point_le.y() = spline::seval(sle, points.row(1), dpoints_ds.row(1), spline_length.head(n), n);
  point_te = 0.5 * (points.col(0) + points.col(n - 1));
  chord = (point_le - point_te).norm();
  tecalc();
  apanel.head(n) = apcalc(points);

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
  VectorXd result = VectorXd::Zero(n);
  //---- set angles of airfoil panels
  for (int i = 0; i < n; i++) {
    Vector2d diff = points.col((i + 1) % n) - points.col(i);
    result[i] = atan2(diff.x(), -diff.y());
  }
  //---- TE panel
  if (sharp) {
    result[n - 1] = std::numbers::pi;
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
    double sin = cross2(delta_later, delta_former) / delta_former.norm() / delta_later.norm();
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

bool XFoil::tecalc() {
  //-------------------------------------------
  //     calculates total and projected TE
  //     areas and TE panel strengths.
  //-------------------------------------------
  double scs, sds;
  //---- set te base vector and te bisector components
  Vector2d tevec = points.col(0) - points.col(n - 1);
  Vector2d dpoint_ds_te = 0.5 * (-dpoints_ds.col(0) + dpoints_ds.col(n - 1));
  //---- normal and streamwise projected TE gap areas
  ante = cross2(dpoint_ds_te, tevec);
  aste = tevec.dot(dpoint_ds_te);
  //---- total TE gap area
  dste = tevec.norm();
  sharp = dste < 0.0001 * chord;
  if (sharp) {
    scs = 1.0;
    sds = 0.0;
  } else {
    scs = ante / dste;
    sds = aste / dste;
  }
  //---- TE panel source and vorticity strengths
  sigte = 0.5 * (surface_vortex(0, 0) - surface_vortex(0, n - 1)) * scs;
  gamte = -.5 * (surface_vortex(0, 0) - surface_vortex(0, n - 1)) * sds;
  return true;
}

bool XFoil::xyWake() {
  //-----------------------------------------------------
  //     sets wake coordinate array for current surface
  //     vorticity and/or mass source distributions.
  //-----------------------------------------------------
  double ds1, sx, sy, smod;
  writeString("   Calculating wake trajectory ...\n");
  ds1 = 0.5 * (spline_length[1] - spline_length[0] + spline_length[n - 1] - spline_length[n - 2]);
  setexp(snew.data() + n, ds1, waklen * chord, nw);
  point_te = 0.5 * (points.col(0) + points.col(n - 1));
  //-- set first wake point a tiny distance behind te
  sx = 0.5 * (dpoints_ds.col(n - 1).y() - dpoints_ds.col(0).y());
  sy = 0.5 * (dpoints_ds.col(0).x() - dpoints_ds.col(n - 1).x());
  smod = sqrt(sx * sx + sy * sy);
  normal_vectors.col(n).x() = sx / smod;
  normal_vectors.col(n).y() = sy / smod;
  points.col(n).x() = point_te.x() - 0.0001 * normal_vectors.col(n).y();
  points.col(n).y() = point_te.y() + 0.0001 * normal_vectors.col(n).x();
  spline_length[n] = spline_length[n - 1];
  //---- calculate streamfunction gradient components at first point
  Vector2d psi = {psilin(points, n, points.col(n), {1.0, 0.0}, false).psi_ni,
                  psilin(points, n, points.col(n), {0.0, 1.0}, false).psi_ni};
  //---- set unit vector normal to wake at first point
  normal_vectors.col(n + 1) = -psi.normalized();
  //---- set angle of wake panel normal
  apanel[n] = atan2(psi.y(), psi.x());
  //---- set rest of wake points
  for (int i = n + 1; i < n + nw; i++) {
    const double ds = snew[i] - snew[i - 1];
    //------ set new point ds downstream of last point
    points.col(i).x() = points.col(i - 1).x() - ds * normal_vectors.col(i - 1).y();
    points.col(i).y() = points.col(i - 1).y() + ds * normal_vectors.col(i - 1).x();
    spline_length[i] = spline_length[i - 1] + ds;
    if (i != n + nw - 1) {
      Vector2d psi2 = {psilin(points, i, points.col(i), {1.0, 0.0}, false).psi_ni,
                       psilin(points, i, points.col(i), {0.0, 1.0}, false).psi_ni};
      normal_vectors.col(i + 1) = -psi2.normalized();
      apanel[i] = atan2(psi2.y(), psi2.x());
    }
  }
  //---- set wake presence flag and corresponding alpha
  lwake = true;
  awake = alfa;
  //---- old source influence matrix is invalid for the new wake geometry
  lwdij = false;
  return true;
}

double XFoil::sign(double a, double b) {
  if (b >= 0.0)
    return fabs(a);
  else
    return -fabs(a);
}
