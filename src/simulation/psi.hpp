#pragma once
#include "psi_result.hpp"

PsiResult psilin(Matrix2Xd points, int i, Vector2d point, Vector2d normal_vector,
                 bool siglin, const VectorXd& spline_length, int total_nodes,
                 const Matrix2Xd& gamu, const Matrix2Xd& surface_vortex,
                 double alfa, double qinf);
PsiResult pswlin(Matrix2Xd points, int i, Vector2d point, Vector2d normal_vector);
PsiResult psisig(Matrix2Xd points, int iNode, int jNode, Vector2d point,
                 Vector2d normal_vector);
PsiResult psi_te(Matrix2Xd points, int i, Vector2d normal_vector);
