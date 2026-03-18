#pragma once

#include "model/boundary_layer/bl_transition_params.hpp"
#include "model/boundary_layer/runtime_types.hpp"
#include "model/boundary_layer/state.hpp"
#include "model/foil/foil.hpp"
#include "solver/boundary_layer/boundary_layer_geometry.hpp"

class BoundaryLayerRuntimeStateOps {
 public:
  static int resetSideState(SidePair<BoundaryLayerLattice> &lattice,
                            BlTransitionParams &bl_transition,
                            FlowRegimeEnum &flow_regime, int side,
                            const Foil &foil,
                            const StagnationResult &stagnation);
  static int readSideStationCount(const SidePair<BoundaryLayerLattice> &lattice,
                                  int side);
  static BoundaryLayerStationReadModel readStationModel(
      const SidePair<BoundaryLayerLattice> &lattice, const Eigen::VectorXd &wgap,
      int side, int stationIndex);
  static BoundaryLayerTrailingEdgeReadModel readTrailingEdgeModel(
      const SidePair<BoundaryLayerLattice> &lattice);
  static double xifset(const SidePair<BoundaryLayerLattice> &lattice,
                       const Foil &foil, const StagnationResult &stagnation,
                       int side);
};
