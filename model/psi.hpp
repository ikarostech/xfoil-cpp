#pragma once
#include "psi_result.hpp"

PsiResult psilin(Matrix2Xd points, int i, Vector2d point, Vector2d normal_vector, bool siglin);
PsiResult pswlin(Matrix2Xd points, int i, Vector2d point, Vector2d normal_vector);
PsiResult psisig(Matrix2Xd points, int iNode, int jNode, Vector2d point, Vector2d normal_vector);
PsiResult psi_te(Matrix2Xd points, int i, Vector2d normal_vector);

