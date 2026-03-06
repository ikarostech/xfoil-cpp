#pragma once

#include <vector>

#include "Eigen/Core"
#include "simulation/march/context.hpp"

class XFoil;

class Marcher
{
public:
  using Matrix3x2d = Eigen::Matrix<double, 3, 2>;
  using Matrix3x2dVector = std::vector<Matrix3x2d>;
};
