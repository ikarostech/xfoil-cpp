#pragma once

#include "Eigen/Core"

struct XFoilClComputation {
  double cl = 0.0;
  double cm = 0.0;
  double cl_alf = 0.0;
  double cl_msq = 0.0;
  double xcp = 0.0;
};

struct XFoilViscalEndResult {
  Eigen::VectorXd inviscidCp;
  Eigen::VectorXd viscousCp;
};

struct XFoilCompressibilityParams {
  double beta = 0.0;
  double beta_msq = 0.0;
  double karmanTsienFactor = 0.0;
  double karmanTsienFactor_msq = 0.0;
  double prandtlGlauertFactor = 0.0;
  double prandtlGlauertFactor_msq = 0.0;
};

struct XFoilPressureCoefficientResult {
  double cp = 0.0;
  double cp_msq = 0.0;
  double cp_velocity_derivative = 0.0;
};
