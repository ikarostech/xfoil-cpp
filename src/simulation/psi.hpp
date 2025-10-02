#pragma once
#include <Eigen/Core>
#include "infrastructure/xfoil_params.h"

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

PsiResult psilin(Matrix2Xd points, int i, Vector2d point, Vector2d normal_vector,
                 bool siglin, const VectorXd& spline_length, int total_nodes,
                 const Matrix2Xd& gamu, const Matrix2Xd& surface_vortex,
                 double alfa, double qinf, const VectorXd& apanel, bool sharp,
                 double ante, double dste, double aste);
PsiResult pswlin(Matrix2Xd points, int i, Vector2d point, Vector2d normal_vector,
                 int total_nodes, int wake_nodes, const VectorXd& apanel);
PsiResult psisig(Matrix2Xd points, int iNode, int jNode, Vector2d point,
                 Vector2d normal_vector, int total_nodes,
                 const VectorXd& apanel);
PsiResult psi_te(Matrix2Xd points, int i, Vector2d normal_vector,
                 int total_nodes, const VectorXd& apanel, bool sharp,
                 double ante, double dste, double aste, const Matrix2Xd& gamu,
                 const Matrix2Xd& surface_vortex);
