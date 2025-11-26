#include "foil_shape.hpp"

Eigen::VectorXd FoilShape::calcSpline() const {
  Eigen::VectorXd s = Eigen::VectorXd::Zero(n);

  for (int i = 1; i < n; i++) {
    s[i] = s[i - 1] + (points.col(i) - points.col(i - 1)).norm();
  }

  return s;
}

Eigen::Matrix2Xd FoilShape::calcDPointsDs() const {
  Eigen::Matrix2Xd result = Eigen::Matrix2Xd::Zero(2, n);
  if (n == 0) {
    return result;
  }

  result.row(0) = spline::splind(points.row(0), spline_length);
  result.row(1) = spline::splind(points.row(1), spline_length);
  return result;
}

Eigen::Matrix2Xd FoilShape::calcNormalVector(const Eigen::Matrix2Xd &dpoints_ds) const {
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

Eigen::VectorXd FoilShape::calcAnglePanel() const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(n);
  if (n == 0) {
    return result;
  }
  //---- set angles of airfoil panels
  for (int i = 0; i < n; i++) {
    Eigen::Vector2d diff = points.col((i + 1) % n) - points.col(i);
    result[i] = std::atan2(diff.x(), -diff.y());
  }
  return result;
}
FoilShape::FoilShape(const Eigen::Matrix2Xd& points, int n) {
  setFoilShape(points, n);
}

void FoilShape::setFoilShape(const Eigen::Matrix2Xd& points, int n) {
  this->points = points;
  this->n = n;
  this->spline_length = calcSpline();
  this->dpoints_ds = calcDPointsDs();
  this->normal_vector = calcNormalVector(this->dpoints_ds);
  this->angle_panel = calcAnglePanel();
}
