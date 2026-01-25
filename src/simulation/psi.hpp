#pragma once

#include <Eigen/Core>

#include "domain/foil/foil.hpp"

using Eigen::Matrix2Xd;
using Eigen::Vector;
using Eigen::Vector2d;
using Eigen::VectorXd;

class PsiResult {
 public:
  double psi = 0.0;
  double psi_ni = 0.0;
  Vector2d qtan = Vector2d::Zero();
  VectorXd dzdg;
  VectorXd dqdg;
  VectorXd dzdm;
  VectorXd dqdm;

  PsiResult() = default;
  PsiResult(int dzdg_size, int dzdm_size) {
    dzdg = VectorXd::Zero(dzdg_size);
    dqdg = VectorXd::Zero(dzdg_size);
    dzdm = VectorXd::Zero(dzdm_size);
    dqdm = VectorXd::Zero(dzdm_size);
  }

  static PsiResult sum(const PsiResult& a, const PsiResult& b) {
    PsiResult result(static_cast<int>(a.dzdg.size()),
                     static_cast<int>(a.dzdm.size()));
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

PsiResult psilin(const Foil& foil, int i, Vector2d point,
                 Vector2d normal_vector, bool siglin,
                 const Matrix2Xd& gamu, const Matrix2Xd& surface_vortex,
                 double alfa, double qinf, const VectorXd& apanel);
PsiResult pswlin(const Foil& foil, int i, Vector2d point,
                 Vector2d normal_vector, const VectorXd& apanel);
PsiResult psisig(const Foil& foil, int iNode, int jNode, Vector2d point,
                 Vector2d normal_vector, const VectorXd& apanel);
PsiResult psi_te(const Foil& foil, int i, Vector2d normal_vector,
                 const VectorXd& apanel, const Matrix2Xd& gamu,
                 const Matrix2Xd& surface_vortex);
