#pragma once

#include <array>
#include <vector>

#include "Eigen/Core"
#include "core/side_pair.hpp"

class Blsolve {
 public:
  using Matrix3x2d = Eigen::Matrix<double, 3, 2>;
  using Matrix3x2dVector = std::vector<Matrix3x2d>;
  using VzMatrix = std::array<std::array<double, 2>, 3>;
  struct VmMatrix {
    int size = 0;
    std::vector<double> data;

    void resize(int new_size) {
      size = new_size;
      data.assign(3 * size * size, 0.0);
    }

    double& at(int k, int i, int j) {
      return data[(k * size + i) * size + j];
    }

    const double& at(int k, int i, int j) const {
      return data[(k * size + i) * size + j];
    }
  };

  struct BlNewtonSystem {
    Matrix3x2dVector va;
    Matrix3x2dVector vb;
    Matrix3x2dVector vdel;
    VmMatrix vm;
    VzMatrix vz{};
  };

  std::vector<Eigen::Matrix<double, 3, 2>> solve(int nsys,
               const SidePair<int>& ivte,
               double vaccel,
               const BlNewtonSystem& bl_newton_system) const;
};
