#pragma once

#include <Eigen/Core>

struct SkinFrictionCoefficients {
  double cfm = 0.0;
  Eigen::Vector2d global = Eigen::Vector2d::Zero();
  Eigen::Matrix<double, 2, 3> station = Eigen::Matrix<double, 2, 3>::Zero();
};
