#pragma once

#include "Eigen/Core"

class XFoil;
enum class SpecTarget { AngleOfAttack, LiftCoefficient };

class InviscidSolver {
 public:
  static bool specal(XFoil& xfoil);
  static bool speccl(XFoil& xfoil);
  static Eigen::Matrix2Xd qiset(double alpha,
                                const Eigen::Matrix2Xd& qinvu);
  static Eigen::VectorXd cpcalc(int n, Eigen::VectorXd q,
                                double qinf, double minf);
  private:
    static bool specConverge(XFoil &xfoil, SpecTarget target);
};
