#pragma once

#include "simulation/BoundaryLayer.hpp"

class BoundaryLayerRuntimeStateOps {
 public:
  static int resetSideState(SidePair<BoundaryLayerLattice> &lattice,
                            BlTransitionParams &bl_transition,
                            FlowRegimeEnum &flow_regime, int side,
                            const Foil &foil,
                            const StagnationResult &stagnation);
  static int readSideStationCount(const SidePair<BoundaryLayerLattice> &lattice,
                                  int side);
  static BoundaryLayerWorkflow::StationReadModel readStationModel(
      const SidePair<BoundaryLayerLattice> &lattice, const Eigen::VectorXd &wgap,
      int side, int stationIndex);
  static BoundaryLayerWorkflow::TrailingEdgeReadModel readTrailingEdgeModel(
      const SidePair<BoundaryLayerLattice> &lattice);
  static double xifset(const SidePair<BoundaryLayerLattice> &lattice,
                       const Foil &foil, const StagnationResult &stagnation,
                       int side);
};
