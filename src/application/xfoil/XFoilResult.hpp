#pragma once

#include "Eigen/Core"

#include "model/coefficient/aero_coefficients.hpp"
#include "numerics/side_pair.hpp"

struct XFoilResult {
  AeroCoefficients aeroCoefficients;
  Eigen::VectorXd inviscidCp;
  Eigen::VectorXd viscousCp;
  SidePair<double> transitionLocations{1.0, 1.0};
  bool converged = false;
};
