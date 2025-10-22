#pragma once

#include <Eigen/Core>

struct BlSystemCoeffs {
  Eigen::Matrix<double, 4, 5> a1;     // coefficients at station 1
  Eigen::Matrix<double, 4, 5> a2;     // coefficients at station 2
  Eigen::Vector<double, 4> rhs;       // residual vector
  Eigen::Vector<double, 4> d_msq;     // d(res)/d(Ma^2)
  Eigen::Vector<double, 4> d_re;      // d(res)/d(Re)
  Eigen::Vector<double, 4> d_xi;      // d(res)/d(xi)
  void clear() {
    a1.setZero();
    a2.setZero();
    rhs.setZero();
    d_msq.setZero();
    d_re.setZero();
    d_xi.setZero();
  }
};
