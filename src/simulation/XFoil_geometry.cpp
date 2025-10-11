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
  if (sharp) {
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

XFoil::TrailingEdgeData XFoil::tecalc(const Matrix2Xd& points,
                                      const Matrix2Xd& dpoints_ds,
                                      const Matrix2Xd& surface_vortex,
                                      int n,
                                      double chord) {
  TrailingEdgeData data{};
  if (n < 2 || points.cols() < n || dpoints_ds.cols() < n) {
    return data;
  }

  const Vector2d tevec = points.col(0) - points.col(n - 1);
  const Vector2d dpoint_ds_te =
      0.5 * (-dpoints_ds.col(0) + dpoints_ds.col(n - 1));

  data.ante = cross2(dpoint_ds_te, tevec);
  data.aste = tevec.dot(dpoint_ds_te);
  data.dste = tevec.norm();

  const bool is_sharp = data.dste < 0.0001 * chord;
  double scs = 0.0;
  double sds = 0.0;
  if (is_sharp) {
    scs = 1.0;
    sds = 0.0;
  } else if (data.dste != 0.0) {
    const double inv_dste = 1.0 / data.dste;
    scs = data.ante * inv_dste;
    sds = data.aste * inv_dste;
  }

  double surface_delta = 0.0;
  if (surface_vortex.rows() > 0 && surface_vortex.cols() >= n) {
    surface_delta = surface_vortex(0, 0) - surface_vortex(0, n - 1);
  }

  data.sharp = is_sharp;
  data.sigte = 0.5 * surface_delta * scs;
  data.gamte = -0.5 * surface_delta * sds;

  return data;
}

void XFoil::updateTrailingEdgeState() {
  const auto data = tecalc(foil.foil_shape.points, foil.foil_shape.dpoints_ds, surface_vortex,
                           foil.foil_shape.n, foil.edge_data.chord);
  ante = data.ante;
  aste = data.aste;
  dste = data.dste;
  sharp = data.sharp;
  sigte = data.sigte;
  gamte = data.gamte;
}

VectorXd XFoil::setexp(double ds1, double smax, int nn) const {
  //........................................................
  //     sets geometriy stretched array s:
  //
  //       s(i+1) - s(i)  =  r * [s(i) - s(i-1)]
  //
  //       ds1   (input)   first s increment:  spline_length[2] -
  //       spline_length[1] smax  (input)   final s value:      s(nn) nn (input)
  //       number of points
  //........................................................
  const int nex = nn - 1;

  const double sigma = smax / ds1;
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

  VectorXd spline_length(nn);
  spline_length[0] = 0.0;
  double ds = ds1;
  for (int i = 1; i < nn; i++) {
    spline_length[i] = spline_length[i - 1] + ds;
    ds *= ratio;
  }
  return spline_length;
}

bool XFoil::xyWake() {
  //-----------------------------------------------------
  //     sets wake coordinate array for current surface
  //     vorticity and/or mass source distributions.
  //-----------------------------------------------------
  const int point_count = foil.foil_shape.n;
  foil.wake_shape.points = Eigen::Matrix2Xd::Zero(2, foil.foil_shape.n + nw);
  foil.wake_shape.n = foil.foil_shape.n + nw;
  foil.wake_shape.normal_vector = Eigen::Matrix2Xd::Zero(2, foil.foil_shape.n + nw);
  foil.wake_shape.spline_length = Eigen::VectorXd::Zero(foil.foil_shape.n + nw);
  foil.wake_shape.points.block(0, 0, 2, foil.foil_shape.n) = foil.foil_shape.points;
  foil.wake_shape.normal_vector.block(0, 0, 2, foil.foil_shape.n) = foil.foil_shape.normal_vector;
  foil.wake_shape.spline_length.head(foil.foil_shape.n) = foil.foil_shape.spline_length;

  double ds1 = 0.5 * (foil.foil_shape.spline_length[1] - foil.foil_shape.spline_length[0] +
                      foil.foil_shape.spline_length[point_count - 1] -
                      foil.foil_shape.spline_length[point_count - 2]);
  const auto wake_spacing = setexp(ds1, foil.edge_data.chord, nw);

  Vector2d tangent_vector =
      (foil.foil_shape.dpoints_ds.col(point_count - 1) - foil.foil_shape.dpoints_ds.col(0)).normalized();
  foil.wake_shape.normal_vector.col(point_count) = Vector2d{tangent_vector.y(), -tangent_vector.x()};
  
  foil.wake_shape.points.col(point_count) = foil.edge_data.point_te + 0.0001 * tangent_vector.col(point_count);
  foil.foil_shape.points.col(point_count) = foil.wake_shape.points.col(point_count);
  foil.wake_shape.spline_length[point_count] = foil.wake_shape.spline_length[point_count - 1];
  //---- calculate streamfunction gradient components at first point
  Vector2d psi = {
      psilin(foil.wake_shape.points, point_count, foil.wake_shape.points.col(point_count), {1.0, 0.0}, false,
             foil.foil_shape.spline_length, point_count,
             gamu, surface_vortex, alfa, qinf, apanel, sharp, ante, dste, aste)
          .psi_ni,
      psilin(foil.wake_shape.points, point_count, foil.wake_shape.points.col(point_count), {0.0, 1.0}, false,
             foil.foil_shape.spline_length, point_count,
             gamu, surface_vortex, alfa, qinf, apanel, sharp, ante, dste, aste)
          .psi_ni};
  //---- set unit vector normal to wake at first point
  foil.wake_shape.normal_vector.col(point_count + 1) = -psi.normalized();
  //---- set angle of wake panel normal
  apanel[point_count] = atan2(psi.y(), psi.x());
  //---- set rest of wake points
  for (int i = point_count + 1; i < point_count + nw; i++) {
    const double ds = wake_spacing[i - point_count] - wake_spacing[i - point_count - 1];
    //------ set new point ds downstream of last point
    foil.wake_shape.points.col(i).x() = foil.wake_shape.points.col(i - 1).x() - ds * foil.wake_shape.normal_vector.col(i - 1).y();
    foil.wake_shape.points.col(i).y() = foil.wake_shape.points.col(i - 1).y() + ds * foil.wake_shape.normal_vector.col(i - 1).x();
    foil.foil_shape.points.col(i) = foil.wake_shape.points.col(i);
    foil.wake_shape.spline_length[i] = foil.wake_shape.spline_length[i - 1] + ds;
    if (i == point_count + nw - 1) {
      break;
    }

    Vector2d psi2 = {
        psilin(foil.wake_shape.points, i, foil.wake_shape.points.col(i), {1.0, 0.0}, false,
                foil.wake_shape.spline_length, point_count,
                gamu, surface_vortex, alfa, qinf, apanel, sharp, ante, dste,
                aste)
            .psi_ni,
        psilin(foil.wake_shape.points, i, foil.wake_shape.points.col(i), {0.0, 1.0}, false,
                foil.wake_shape.spline_length, point_count,
                gamu, surface_vortex, alfa, qinf, apanel, sharp, ante, dste,
                aste)
            .psi_ni};
    foil.wake_shape.normal_vector.col(i + 1) = -psi2.normalized();
    apanel[i] = atan2(psi2.y(), psi2.x());
    
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
