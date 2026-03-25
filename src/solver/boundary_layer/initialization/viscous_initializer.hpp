#pragma once

#include "solver/boundary_layer/boundary_layer.hpp"
#include "solver/boundary_layer/runtime/state.hpp"
#include "solver/boundary_layer/initialization/setbl.hpp"

class BoundaryLayerInitializer {
 public:
  static SetblOutputView run(
      BoundaryLayer &boundaryLayer, const FlowState &analysis_state,
      const AeroCoefficients &aero_coeffs,
      double acrit, const Foil &foil, const StagnationFeature &stagnation,
      const Eigen::MatrixXd &dij, bool bl_initialized) {
    return BoundaryLayerSetblUseCase{}.run(
                                           boundaryLayer,
                                           boundaryLayer.currentProfiles(),
                                           analysis_state,
                                           aero_coeffs, acrit, foil,
                                           stagnation, dij, bl_initialized);
  }

  static void applyOutput(BoundaryLayer &boundaryLayer,
                          SetblOutputView &output) {
    boundaryLayer.applyInitializationState(output.blCompressibility,
                                           output.blReynolds,
                                           output.blTransition,
                                           output.flowRegime,
                                           std::move(output.profiles));
  }

  static int resetSideState(BoundaryLayer &boundaryLayer, int side,
                            const Foil &foil,
                            const StagnationFeature &stagnation) {
    return boundaryLayer.resetSideState(side, foil, stagnation);
  }

  static double xifset(const BoundaryLayer &boundaryLayer, const Foil &foil,
                       const StagnationFeature &stagnation, int side) {
    return boundaryLayer.computeForcedTransitionArcLength(foil, stagnation,
                                                          side);
  }
};
