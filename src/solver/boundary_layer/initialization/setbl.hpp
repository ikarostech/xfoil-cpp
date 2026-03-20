#pragma once

#include "model/boundary_layer/state.hpp"
#include "solver/boundary_layer/boundary_layer_builder.hpp"
#include "solver/boundary_layer/boundary_layer.hpp"
#include "model/coefficient/aero_coefficients.hpp"
#include "model/flow_state.hpp"

class BoundaryLayerSetblUseCase {
 public:
  SetblOutputView run(
      BoundaryLayer &boundaryLayer,
      SidePairRef<const BoundaryLayerSideProfiles> profiles,
      const FlowState &analysis_state, const AeroCoefficients &aero_coeffs,
      double acrit, const Foil &foil, const StagnationResult &stagnation,
      const Eigen::MatrixXd &dij, bool bl_initialized) const;
};
