#pragma once

#include "Eigen/Core"

class XFoilAnalysis;
enum class SpecTarget { AngleOfAttack, LiftCoefficient };

class InviscidSolver {
public:
  static bool specal(XFoilAnalysis &xfoil);
  static bool speccl(XFoilAnalysis &xfoil);
  static Eigen::Matrix2Xd qiset(double alpha, const Eigen::Matrix2Xd &qinvu);
  static Eigen::VectorXd cpcalc(int n, Eigen::VectorXd q, double qinf,
                                double minf);

private:
  static bool specConverge(XFoilAnalysis &xfoil, SpecTarget target);
};
