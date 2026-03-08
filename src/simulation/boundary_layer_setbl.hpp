#pragma once

#include "simulation/boundary_layer_state.hpp"
#include "domain/boundary_layer/boundary_layer_builder.hpp"
#include "domain/coefficient/aero_coefficients.hpp"
#include "domain/flow_state.hpp"
#include "simulation/BoundaryLayer.hpp"

SetblOutputView runBoundaryLayerSetbl(
    BoundaryLayerWorkflow &workflow,
    SidePairRef<const BoundaryLayerSideProfiles> profiles,
    const FlowState &analysis_state, const AeroCoefficients &aero_coeffs,
    double acrit, const Foil &foil, const StagnationResult &stagnation,
    const Eigen::MatrixXd &dij, bool bl_initialized);
