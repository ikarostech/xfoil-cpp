#pragma once

#include <array>
#include <vector>

#include "Eigen/Core"
#include "core/side_pair.hpp"

class Blsolve {
 public:
  using Matrix3x2d = Eigen::Matrix<double, 3, 2>;
  using Matrix3x2dVector = std::vector<Matrix3x2d>;
  struct Output {
    std::vector<double> vm;
    Matrix3x2dVector vdel;
  };

  Output solve(int nsys,
               const SidePair<int>& ivte,
               double vaccel,
               const Matrix3x2dVector& va,
               const Matrix3x2dVector& vb,
               const std::vector<double>& vm,
               int vm_size,
               Matrix3x2dVector vdel,
               const std::array<std::array<double, 2>, 3>& vz) const;
};
