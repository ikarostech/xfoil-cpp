/**
 * BL Newton system custom solver split from XFoil.cpp
 */

#include "XFoil.h"

namespace {
inline void plu3x3(double m[3][3], int piv[3]) {
  piv[0] = 0;
  piv[1] = 1;
  piv[2] = 2;
  int p = 0;
  if (std::abs(m[1][0]) > std::abs(m[0][0]))
    p = 1;
  if (std::abs(m[2][0]) > std::abs(m[p][0]))
    p = 2;

  if (p != 0) {
    std::swap(m[0], m[p]);
    std::swap(piv[0], piv[p]);
  }
  double inv = 1.0 / m[0][0];
  m[1][0] *= inv;
  m[2][0] *= inv;
  m[1][1] -= m[1][0] * m[0][1];
  m[1][2] -= m[1][0] * m[0][2];
  m[2][1] -= m[2][0] * m[0][1];
  m[2][2] -= m[2][0] * m[0][2];
  if (std::abs(m[2][1]) > std::abs(m[1][1])) {
    std::swap(m[1][0], m[2][0]);
    std::swap(m[1][1], m[2][1]);
    std::swap(m[1][2], m[2][2]);
    std::swap(piv[1], piv[2]);
  }
  inv = 1.0 / m[1][1];
  m[2][1] *= inv;
  m[2][2] -= m[2][1] * m[1][2];
}

inline void luSolve3x3(const double m[3][3], const int piv[3], double b[3]) {
  double x0 = b[piv[0]];
  double x1 = b[piv[1]] - m[1][0] * x0;
  double x2 = b[piv[2]] - m[2][0] * x0 - m[2][1] * x1;
  x2 /= m[2][2];
  x1 = (x1 - m[1][2] * x2) / m[1][1];
  x0 = (x0 - m[0][1] * x1 - m[0][2] * x2) / m[0][0];
  b[0] = x0;
  b[1] = x1;
  b[2] = x2;
}

inline void luSolve3x3x6(const double m[3][3], const int piv[3],
                         double b[3][6]) {
  for (int j = 0; j < 6; ++j) {
    double x0 = b[piv[0]][j];
    double x1 = b[piv[1]][j] - m[1][0] * x0;
    double x2 = b[piv[2]][j] - m[2][0] * x0 - m[2][1] * x1;
    x2 /= m[2][2];
    x1 = (x1 - m[1][2] * x2) / m[1][1];
    x0 = (x0 - m[0][1] * x1 - m[0][2] * x2) / m[0][0];
    b[0][j] = x0;
    b[1][j] = x1;
    b[2][j] = x2;
  }
}
} // namespace

bool XFoil::blsolve() {

  auto eliminateVaBlock = [&](int iv, int ivp) {
    double D[3][3] = {{va[iv](0, 0), va[iv](0, 1), vm[0][iv][iv]},
                      {va[iv](1, 0), va[iv](1, 1), vm[1][iv][iv]},
                      {va[iv](2, 0), va[iv](2, 1), vm[2][iv][iv]}};
    int piv[3];
    plu3x3(D, piv);

    double rhs[3][6] = {};
    int cols = 0;
    for (int offset = 0; offset < 3 && iv + offset <= nsys; ++offset, ++cols) {
      rhs[0][cols] = vm[0][iv + offset][iv];
      rhs[1][cols] = vm[1][iv + offset][iv];
      rhs[2][cols] = vm[2][iv + offset][iv];
    }
    rhs[0][cols] = vdel[iv](0, 0);
    rhs[1][cols] = vdel[iv](1, 0);
    rhs[2][cols] = vdel[iv](2, 0);
    ++cols;
    rhs[0][cols] = vdel[iv](0, 1);
    rhs[1][cols] = vdel[iv](1, 1);
    rhs[2][cols] = vdel[iv](2, 1);

    luSolve3x3x6(D, piv, rhs);

    int idx = 0;
    for (int offset = 0; offset < 3 && iv + offset <= nsys; ++offset, ++idx) {
      vm[0][iv + offset][iv] = rhs[0][idx];
      vm[1][iv + offset][iv] = rhs[1][idx];
      vm[2][iv + offset][iv] = rhs[2][idx];
    }
    vdel[iv](0, 0) = rhs[0][idx];
    vdel[iv](1, 0) = rhs[1][idx];
    vdel[iv](2, 0) = rhs[2][idx];
    ++idx;
    vdel[iv](0, 1) = rhs[0][idx];
    vdel[iv](1, 1) = rhs[1][idx];
    vdel[iv](2, 1) = rhs[2][idx];

    for (int l = iv + 3; l <= nsys; ++l) {
      double col[3] = {vm[0][l][iv], vm[1][l][iv], vm[2][l][iv]};
      luSolve3x3(D, piv, col);
      vm[0][l][iv] = col[0];
      vm[1][l][iv] = col[1];
      vm[2][l][iv] = col[2];
    }
  };

  auto eliminateVbBlock = [&](int iv, int ivp, int ivte1) {
    double D[3][3] = {{vb[ivp](0, 0), vb[ivp](0, 1), vm[0][iv][ivp]},
                      {vb[ivp](1, 0), vb[ivp](1, 1), vm[1][iv][ivp]},
                      {vb[ivp](2, 0), vb[ivp](2, 1), vm[2][iv][ivp]}};

    double col[3];
    for (int l = ivp; l <= nsys; ++l) {
      col[0] = vm[0][l][iv];
      col[1] = vm[1][l][iv];
      col[2] = vm[2][l][iv];
      for (int k = 0; k < 3; ++k)
        vm[k][l][ivp] -= D[k][0] * col[0] + D[k][1] * col[1] + D[k][2] * col[2];
    }

    col[0] = vdel[iv](0, 0);
    col[1] = vdel[iv](1, 0);
    col[2] = vdel[iv](2, 0);
    for (int k = 0; k < 3; ++k)
      vdel[ivp](k, 0) -= D[k][0] * col[0] + D[k][1] * col[1] + D[k][2] * col[2];

    col[0] = vdel[iv](0, 1);
    col[1] = vdel[iv](1, 1);
    col[2] = vdel[iv](2, 1);
    for (int k = 0; k < 3; ++k)
      vdel[ivp](k, 1) -= D[k][0] * col[0] + D[k][1] * col[1] + D[k][2] * col[2];

    if (iv == ivte1) {
      int ivz = boundaryLayerWorkflow.lattice.bottom.stationToSystem[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex];
      double Dz[3][2] = {{vz[0][0], vz[0][1]},
                         {vz[1][0], vz[1][1]},
                         {vz[2][0], vz[2][1]}};

      for (int l = ivp; l <= nsys; ++l) {
        col[0] = vm[0][l][iv];
        col[1] = vm[1][l][iv];
        for (int k = 0; k < 3; ++k)
          vm[k][l][ivz] -= Dz[k][0] * col[0] + Dz[k][1] * col[1];
      }

      col[0] = vdel[iv](0, 0);
      col[1] = vdel[iv](1, 0);
      for (int k = 0; k < 3; ++k)
        vdel[ivz](k, 0) -= Dz[k][0] * col[0] + Dz[k][1] * col[1];

      col[0] = vdel[iv](0, 1);
      col[1] = vdel[iv](1, 1);
      for (int k = 0; k < 3; ++k)
        vdel[ivz](k, 1) -= Dz[k][0] * col[0] + Dz[k][1] * col[1];
    }
  };

  auto eliminateLowerVmColumn = [&](int iv, int ivp) {
    for (int kv = iv + 2; kv <= nsys; kv++) {
      double vtmp1 = vm[0][iv][kv];
      double vtmp2 = vm[1][iv][kv];
      double vtmp3 = vm[2][iv][kv];
      if (fabs(vtmp1) > VAccel()) {
        for (int l = ivp; l <= nsys; l++)
          vm[0][l][kv] -= vtmp1 * vm[2][l][iv];
        vdel[kv](0, 0) -= vtmp1 * vdel[iv](2, 0);
        vdel[kv](0, 1) -= vtmp1 * vdel[iv](2, 1);
      }
      if (fabs(vtmp2) > VAccel()) {
        for (int l = ivp; l <= nsys; l++)
          vm[1][l][kv] -= vtmp2 * vm[2][l][iv];
        vdel[kv](1, 0) -= vtmp2 * vdel[iv](2, 0);
        vdel[kv](1, 1) -= vtmp2 * vdel[iv](2, 1);
      }
      if (fabs(vtmp3) > VAccel()) {
        for (int l = ivp; l <= nsys; l++)
          vm[2][l][kv] -= vtmp3 * vm[2][l][iv];
        vdel[kv](2, 0) -= vtmp3 * vdel[iv](2, 0);
        vdel[kv](2, 1) -= vtmp3 * vdel[iv](2, 1);
      }
    }
  };

  auto backSubstitute = [&]() {
    for (int iv = nsys; iv >= 2; iv--) {
      double vtmp = vdel[iv](2, 0);
      for (int kv = iv - 1; kv >= 1; kv--) {
        vdel[kv](0, 0) -= vm[0][iv][kv] * vtmp;
        vdel[kv](1, 0) -= vm[1][iv][kv] * vtmp;
        vdel[kv](2, 0) -= vm[2][iv][kv] * vtmp;
      }
      vtmp = vdel[iv](2, 1);
      for (int kv = iv - 1; kv >= 1; kv--) {
        vdel[kv](0, 1) -= vm[0][iv][kv] * vtmp;
        vdel[kv](1, 1) -= vm[1][iv][kv] * vtmp;
        vdel[kv](2, 1) -= vm[2][iv][kv] * vtmp;
      }
    }
  };

  int ivte1 = boundaryLayerWorkflow.lattice.top.stationToSystem[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex];
  for (int iv = 1; iv <= nsys; iv++) {
    int ivp = iv + 1;
    eliminateVaBlock(iv, ivp);
    if (iv != nsys) {
      eliminateVbBlock(iv, ivp, ivte1);
      if (ivp != nsys)
        eliminateLowerVmColumn(iv, ivp);
    }
  }

  backSubstitute();
  return true;
}
