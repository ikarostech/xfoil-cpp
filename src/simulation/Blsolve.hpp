#pragma once

#include <array>
#include <vector>

#include "Eigen/Core"
#include "core/side_pair.hpp"
#include "simulation/XFoil.h"

class Blsolve {
 public:
  using Matrix3x2d = Eigen::Matrix<double, 3, 2>;
  using Matrix3x2dVector = std::vector<Matrix3x2d>;
  struct Output {
    XFoil::VmMatrix vm;
    Matrix3x2dVector vdel;
  };

  Output solve(int nsys,
               const SidePair<int>& ivte,
               double vaccel,
               const XFoil::BlNewtonSystem& bl_newton_system) const;
};
