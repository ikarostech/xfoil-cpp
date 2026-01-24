#pragma once

#include <Eigen/Core>
#include <Eigen/LU>


/**
 * Holds aerodynamic influence data derived purely from a foil geometry.
 * These matrices/vectors can be reused as long as the geometry does not change.
 */
class FoilAerodynamicCache {
 public:
  Eigen::Matrix2Xd gamu;
  Eigen::Matrix2Xd qinvu;
  Eigen::FullPivLU<Eigen::MatrixXd> psi_gamma_lu;
  Eigen::MatrixXd bij;
  Eigen::MatrixXd dij;
};
