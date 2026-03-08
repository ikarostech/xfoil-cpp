#pragma once

#include "simulation/BoundaryLayer.hpp"
#include "simulation/boundary_layer_runtime_state.hpp"
#include "simulation/boundary_layer_setbl.hpp"

class BoundaryLayerInitializer {
 public:
  static SetblOutputView run(
      BoundaryLayerWorkflow &workflow,
      SidePairRef<const BoundaryLayerSideProfiles> profiles,
      const FlowState &analysis_state, const AeroCoefficients &aero_coeffs,
      double acrit, const Foil &foil, const StagnationResult &stagnation,
      const Eigen::MatrixXd &dij, bool bl_initialized) {
    return runBoundaryLayerSetbl(workflow, profiles, analysis_state,
                                 aero_coeffs, acrit, foil, stagnation, dij,
                                 bl_initialized);
  }

  static void applyOutput(BoundaryLayerWorkflow &workflow,
                          SetblOutputView &output) {
    workflow.blCompressibility = output.blCompressibility;
    workflow.blReynolds = output.blReynolds;
    workflow.lattice.top.profiles = std::move(output.profiles.top);
    workflow.lattice.bottom.profiles = std::move(output.profiles.bottom);
    workflow.flowRegime = output.flowRegime;
    workflow.blTransition = output.blTransition;
  }

  static int resetSideState(BoundaryLayerWorkflow &workflow, int side,
                            const Foil &foil,
                            const StagnationResult &stagnation) {
    return BoundaryLayerRuntimeStateOps::resetSideState(
        workflow.lattice, workflow.blTransition, workflow.flowRegime, side,
        foil, stagnation);
  }

  static double xifset(const BoundaryLayerWorkflow &workflow, const Foil &foil,
                       const StagnationResult &stagnation, int side) {
    return BoundaryLayerRuntimeStateOps::xifset(workflow.lattice, foil,
                                                stagnation, side);
  }
};
