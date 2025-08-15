#pragma once

PsiResult psilin(int i, Vector2d point, Vector2d normal_vector, bool siglin);
PsiResult pswlin(int i, Vector2d point, Vector2d normal_vector);
PsiResult psisig(int iNode, int jNode, Vector2d point, Vector2d normal_vector);
PsiResult psi_te(int i, Vector2d point, Vector2d normal_vector);

