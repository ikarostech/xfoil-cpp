#pragma once

#include "Eigen/Core"

class XFoil;

class InviscidSolver {
 public:
  static bool specal(XFoil& xfoil);
  static bool speccl(XFoil& xfoil);
  static Eigen::Matrix2Xd qiset(const XFoil& xfoil);
  static Eigen::VectorXd cpcalc(XFoil& xfoil, int n, Eigen::VectorXd q,
                                double qinf, double minf);
};
