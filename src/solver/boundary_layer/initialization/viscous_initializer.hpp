#pragma once

#include "solver/boundary_layer/workflow/workflow.hpp"
#include "solver/boundary_layer/runtime/state.hpp"
#include "solver/boundary_layer/initialization/setbl.hpp"

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
    workflow.compressibility() = output.blCompressibility;
    workflow.reynolds() = output.blReynolds;
    workflow.applyProfiles(std::move(output.profiles));
    workflow.flowRegime() = output.flowRegime;
    workflow.transition() = output.blTransition;
  }

  static int resetSideState(BoundaryLayerWorkflow &workflow, int side,
                            const Foil &foil,
                            const StagnationResult &stagnation) {
    return workflow.resetSideState(side, foil, stagnation);
  }

  static double xifset(const BoundaryLayerWorkflow &workflow, const Foil &foil,
                       const StagnationResult &stagnation, int side) {
    return workflow.computeForcedTransitionArcLength(foil, stagnation, side);
  }
};
