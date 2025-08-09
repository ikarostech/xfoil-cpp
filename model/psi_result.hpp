#ifndef MODEL_PSI_RESULT_HPP
#define MODEL_PSI_RESULT_HPP

#include <Eigen/Core>
#include "../xfoil_params.h"

using namespace Eigen;

class PsiResult {
 public:
  double psi;
  double psi_ni;
  Vector2d qtan;
  Vector<double, IQX> dzdg;
  Vector<double, IQX> dqdg;
  Vector<double, IZX> dzdm;
  Vector<double, IZX> dqdm;

  PsiResult() {
    psi = 0;
    psi_ni = 0;
    qtan = Vector2d::Zero();
    dzdg = Vector<double, IQX>::Zero();
    dqdg = Vector<double, IQX>::Zero();
    dzdm = Vector<double, IZX>::Zero();
    dqdm = Vector<double, IZX>::Zero();
  }

  static PsiResult sum(PsiResult a, PsiResult b) {
    PsiResult result;
    result.psi = a.psi + b.psi;
    result.psi_ni = a.psi_ni + b.psi_ni;
    result.qtan = a.qtan + b.qtan;
    result.dzdg = a.dzdg + b.dzdg;
    result.dqdg = a.dqdg + b.dqdg;
    result.dzdm = a.dzdm + b.dzdm;
    result.dqdm = a.dqdm + b.dqdm;
    return result;
  }
};

#endif  // MODEL_PSI_RESULT_HPP
