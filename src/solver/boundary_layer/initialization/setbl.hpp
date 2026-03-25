#pragma once

#include "model/boundary_layer/state.hpp"
#include "model/coefficient/aero_coefficients.hpp"
#include "model/flow_state.hpp"
#include "solver/boundary_layer/boundary_layer.hpp"
#include "solver/boundary_layer/boundary_layer_builder.hpp"

class BoundaryLayerSetblUseCase {
  public:
    SetblOutputView run(BoundaryLayer &boundaryLayer, SidePairRef<const BoundaryLayerSideState> profiles,
                        const FlowState &analysis_state, const AeroCoefficients &aero_coeffs, double acrit,
                        const Foil &foil, const StagnationFeature &stagnation, const Eigen::MatrixXd &dij,
                        bool bl_initialized) const;
};
