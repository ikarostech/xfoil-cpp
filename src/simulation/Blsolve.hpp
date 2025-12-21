#pragma once

#include <array>
#include <vector>

#include "Eigen/Core"
#include "infrastructure/xfoil_params.h"

class Blsolve {
 public:
  using Matrix3x2d = Eigen::Matrix<double, 3, 2>;
  using Matrix3x2dVector = std::vector<Matrix3x2d>;
  struct Output {
    std::vector<double> vm;
    Matrix3x2dVector vdel;
  };

  Output solve(int nsys,
               int ivte1,
               int ivz,
               double vaccel,
               const Matrix3x2dVector& va,
               const Matrix3x2dVector& vb,
               const double (&vm)[3][IZX][IZX],
               Matrix3x2dVector vdel,
               const double (&vz)[3][2]) const;
};
