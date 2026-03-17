#pragma once

#include "solver/boundary_layer/boundary_layer_workflow.hpp"
#include "solver/boundary_layer/boundary_layer_runtime_state.hpp"
#include "solver/boundary_layer/boundary_layer_setbl.hpp"

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
    auto &state_store = workflow.stateStore();
    state_store.blCompressibility = output.blCompressibility;
    state_store.blReynolds = output.blReynolds;
    state_store.lattice.top.profiles = std::move(output.profiles.top);
    state_store.lattice.bottom.profiles = std::move(output.profiles.bottom);
    state_store.flowRegime = output.flowRegime;
    state_store.blTransition = output.blTransition;
  }

  static int resetSideState(BoundaryLayerWorkflow &workflow, int side,
                            const Foil &foil,
                            const StagnationResult &stagnation) {
    auto &state_store = workflow.stateStore();
    return BoundaryLayerRuntimeStateOps::resetSideState(
        state_store.lattice, state_store.blTransition, state_store.flowRegime, side,
        foil, stagnation);
  }

  static double xifset(const BoundaryLayerWorkflow &workflow, const Foil &foil,
                       const StagnationResult &stagnation, int side) {
    return BoundaryLayerRuntimeStateOps::xifset(workflow.stateStore().lattice, foil,
                                                stagnation, side);
  }
};
