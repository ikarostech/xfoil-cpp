#pragma once

#include "simulation/boundary_layer_state.hpp"

struct BlSystemCoeffs;

class MrchueLinearSystemOps {
 public:
  static double readNewtonRhs(const BlSystemCoeffs &blc, int row) {
    return blc.rhs[row];
  }

  static void solveDirect(BlSystemCoeffs &blc) {
    blc.a2(3, 0) = 0.0;
    blc.a2(3, 1) = 0.0;
    blc.a2(3, 2) = 0.0;
    blc.a2(3, 3) = 1.0;
    blc.rhs[3] = 0.0;
    blc.rhs = blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);
  }

  static void solveInverse(BlSystemCoeffs &blc, const BoundaryLayerState &state,
                           double htarg) {
    blc.a2(3, 0) = 0.0;
    blc.a2(3, 1) = state.station2.hkz.vector[0];
    blc.a2(3, 2) = state.station2.hkz.vector[1];
    blc.a2(3, 3) = state.station2.hkz.vector[2];
    blc.rhs[3] = htarg - state.station2.hkz.scalar;
    blc.rhs = blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);
  }
};
