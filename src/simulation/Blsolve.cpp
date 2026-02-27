/**
 * BL Newton system custom solver extracted from XFoil.cpp
 */

#include "simulation/Blsolve.hpp"

#include <cmath>
#include <utility>

namespace {
constexpr double kMinPivot = 1.0e-30;

inline double safePivot(double value) {
  if (std::abs(value) < kMinPivot) {
    return (value < 0.0) ? -kMinPivot : kMinPivot;
  }
  return value;
}

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
  double inv = 1.0 / safePivot(m[0][0]);
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
  inv = 1.0 / safePivot(m[1][1]);
  m[2][1] *= inv;
  m[2][2] -= m[2][1] * m[1][2];
}

inline void luSolve3x3(const double m[3][3], const int piv[3], double b[3]) {
  double x0 = b[piv[0]];
  double x1 = b[piv[1]] - m[1][0] * x0;
  double x2 = b[piv[2]] - m[2][0] * x0 - m[2][1] * x1;
  x2 /= safePivot(m[2][2]);
  x1 = (x1 - m[1][2] * x2) / safePivot(m[1][1]);
  x0 = (x0 - m[0][1] * x1 - m[0][2] * x2) / safePivot(m[0][0]);
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
    x2 /= safePivot(m[2][2]);
    x1 = (x1 - m[1][2] * x2) / safePivot(m[1][1]);
    x0 = (x0 - m[0][1] * x1 - m[0][2] * x2) / safePivot(m[0][0]);
    b[0][j] = x0;
    b[1][j] = x1;
    b[2][j] = x2;
  }
}

inline int vmIndex(int size, int k, int i, int j) {
  return (k * size + i) * size + j;
}
} // namespace

std::vector<Eigen::Matrix<double, 3, 2>> Blsolve::solve(int nsys,
                               const SidePair<int>& ivte,
                               double vaccel,
                               const BlNewtonSystem& bl_newton_system) const {
  const Matrix3x2dVector& va = bl_newton_system.va;
  const Matrix3x2dVector& vb = bl_newton_system.vb;
  VmMatrix vm = bl_newton_system.vm;
  Matrix3x2dVector vdel = bl_newton_system.vdel;
  VzMatrix vz = bl_newton_system.vz;
  auto eliminateVaBlock = [&](int iv) {
    double D[3][3] = {{va[iv](0, 0), va[iv](0, 1), vm.data[vmIndex(vm.size, 0, iv, iv)]},
                      {va[iv](1, 0), va[iv](1, 1), vm.data[vmIndex(vm.size, 1, iv, iv)]},
                      {va[iv](2, 0), va[iv](2, 1), vm.data[vmIndex(vm.size, 2, iv, iv)]}};
    int piv[3];
    plu3x3(D, piv);

    double rhs[3][6] = {};
    int cols = 0;
    for (int offset = 0; offset < 3 && iv + offset < nsys; ++offset, ++cols) {
      rhs[0][cols] = vm.data[vmIndex(vm.size, 0, iv + offset, iv)];
      rhs[1][cols] = vm.data[vmIndex(vm.size, 1, iv + offset, iv)];
      rhs[2][cols] = vm.data[vmIndex(vm.size, 2, iv + offset, iv)];
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
    for (int offset = 0; offset < 3 && iv + offset < nsys; ++offset, ++idx) {
      vm.data[vmIndex(vm.size, 0, iv + offset, iv)] = rhs[0][idx];
      vm.data[vmIndex(vm.size, 1, iv + offset, iv)] = rhs[1][idx];
      vm.data[vmIndex(vm.size, 2, iv + offset, iv)] = rhs[2][idx];
    }
    vdel[iv](0, 0) = rhs[0][idx];
    vdel[iv](1, 0) = rhs[1][idx];
    vdel[iv](2, 0) = rhs[2][idx];
    ++idx;
    vdel[iv](0, 1) = rhs[0][idx];
    vdel[iv](1, 1) = rhs[1][idx];
    vdel[iv](2, 1) = rhs[2][idx];

    for (int l = iv + 3; l < nsys; ++l) {
      double col[3] = {vm.data[vmIndex(vm.size, 0, l, iv)],
                       vm.data[vmIndex(vm.size, 1, l, iv)],
                       vm.data[vmIndex(vm.size, 2, l, iv)]};
      luSolve3x3(D, piv, col);
      vm.data[vmIndex(vm.size, 0, l, iv)] = col[0];
      vm.data[vmIndex(vm.size, 1, l, iv)] = col[1];
      vm.data[vmIndex(vm.size, 2, l, iv)] = col[2];
    }
  };

  auto eliminateVbBlock = [&](int iv, int ivp) {
    double D[3][3] = {{vb[ivp](0, 0), vb[ivp](0, 1), vm.data[vmIndex(vm.size, 0, iv, ivp)]},
                      {vb[ivp](1, 0), vb[ivp](1, 1), vm.data[vmIndex(vm.size, 1, iv, ivp)]},
                      {vb[ivp](2, 0), vb[ivp](2, 1), vm.data[vmIndex(vm.size, 2, iv, ivp)]}};

    double col[3];
    for (int l = ivp; l < nsys; ++l) {
      col[0] = vm.data[vmIndex(vm.size, 0, l, iv)];
      col[1] = vm.data[vmIndex(vm.size, 1, l, iv)];
      col[2] = vm.data[vmIndex(vm.size, 2, l, iv)];
      for (int k = 0; k < 3; ++k)
        vm.data[vmIndex(vm.size, k, l, ivp)] -=
            D[k][0] * col[0] + D[k][1] * col[1] + D[k][2] * col[2];
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

    if (iv == ivte.top) {
      double Dz[3][2] = {{vz[0][0], vz[0][1]},
                         {vz[1][0], vz[1][1]},
                         {vz[2][0], vz[2][1]}};

      for (int l = ivp; l < nsys; ++l) {
        col[0] = vm.data[vmIndex(vm.size, 0, l, iv)];
        col[1] = vm.data[vmIndex(vm.size, 1, l, iv)];
        for (int k = 0; k < 3; ++k)
          vm.data[vmIndex(vm.size, k, l, ivte.bottom)] -=
              Dz[k][0] * col[0] + Dz[k][1] * col[1];
      }

      col[0] = vdel[iv](0, 0);
      col[1] = vdel[iv](1, 0);
      for (int k = 0; k < 3; ++k)
        vdel[ivte.bottom](k, 0) -= Dz[k][0] * col[0] + Dz[k][1] * col[1];

      col[0] = vdel[iv](0, 1);
      col[1] = vdel[iv](1, 1);
      for (int k = 0; k < 3; ++k)
        vdel[ivte.bottom](k, 1) -= Dz[k][0] * col[0] + Dz[k][1] * col[1];
    }
  };

  auto eliminateLowerVmColumn = [&](int iv, int ivp) {
    for (int kv = iv + 2; kv < nsys; kv++) {
      double vtmp1 = vm.data[vmIndex(vm.size, 0, iv, kv)];
      double vtmp2 = vm.data[vmIndex(vm.size, 1, iv, kv)];
      double vtmp3 = vm.data[vmIndex(vm.size, 2, iv, kv)];
      if (std::fabs(vtmp1) > vaccel) {
        for (int l = ivp; l < nsys; l++)
          vm.data[vmIndex(vm.size, 0, l, kv)] -=
              vtmp1 * vm.data[vmIndex(vm.size, 2, l, iv)];
        vdel[kv](0, 0) -= vtmp1 * vdel[iv](2, 0);
        vdel[kv](0, 1) -= vtmp1 * vdel[iv](2, 1);
      }
      if (std::fabs(vtmp2) > vaccel) {
        for (int l = ivp; l < nsys; l++)
          vm.data[vmIndex(vm.size, 1, l, kv)] -=
              vtmp2 * vm.data[vmIndex(vm.size, 2, l, iv)];
        vdel[kv](1, 0) -= vtmp2 * vdel[iv](2, 0);
        vdel[kv](1, 1) -= vtmp2 * vdel[iv](2, 1);
      }
      if (std::fabs(vtmp3) > vaccel) {
        for (int l = ivp; l < nsys; l++)
          vm.data[vmIndex(vm.size, 2, l, kv)] -=
              vtmp3 * vm.data[vmIndex(vm.size, 2, l, iv)];
        vdel[kv](2, 0) -= vtmp3 * vdel[iv](2, 0);
        vdel[kv](2, 1) -= vtmp3 * vdel[iv](2, 1);
      }
    }
  };

  auto backSubstitute = [&]() {
    for (int iv = nsys - 1; iv >= 2; iv--) {
      double vtmp = vdel[iv](2, 0);
      for (int kv = iv - 1; kv >= 1; kv--) {
        vdel[kv](0, 0) -= vm.data[vmIndex(vm.size, 0, iv, kv)] * vtmp;
        vdel[kv](1, 0) -= vm.data[vmIndex(vm.size, 1, iv, kv)] * vtmp;
        vdel[kv](2, 0) -= vm.data[vmIndex(vm.size, 2, iv, kv)] * vtmp;
      }
      vtmp = vdel[iv](2, 1);
      for (int kv = iv - 1; kv >= 1; kv--) {
        vdel[kv](0, 1) -= vm.data[vmIndex(vm.size, 0, iv, kv)] * vtmp;
        vdel[kv](1, 1) -= vm.data[vmIndex(vm.size, 1, iv, kv)] * vtmp;
        vdel[kv](2, 1) -= vm.data[vmIndex(vm.size, 2, iv, kv)] * vtmp;
      }
    }
  };

  for (int iv = 1; iv < nsys; iv++) {
    int ivp = iv + 1;
    eliminateVaBlock(iv);
    if (iv != nsys - 1) {
      eliminateVbBlock(iv, ivp);
      if (ivp != nsys - 1)
        eliminateLowerVmColumn(iv, ivp);
    }
  }

  backSubstitute();
  return vdel;
}
