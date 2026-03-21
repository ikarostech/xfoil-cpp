#pragma once

#include "model/coefficient/aero_coefficients.hpp"
#include "model/flow_state.hpp"
#include "model/foil/foil.hpp"
#include "solver/boundary_layer/boundary_layer.hpp"
#include "solver/boundary_layer/boundary_layer_geometry.hpp"
#include "solver/xfoil/viscous_update.hpp"

struct SetblOutputView;

namespace boundary_layer_iteration_ops {

SetblOutputView initializeNewtonSystem(BoundaryLayer &boundary_layer,
                                       const FlowState &analysis_state,
                                       const AeroCoefficients &aero_coeffs,
                                       double acrit, const Foil &foil,
                                       const StagnationResult &stagnation,
                                       const Eigen::MatrixXd &dij,
                                       bool bl_initialized);

BoundaryLayerMatrix3x2dVector solveNewtonStep(
    const BoundaryLayer &boundary_layer, double vaccel,
    const SetblOutputView &setbl_output);

void applyIterationUpdate(BoundaryLayer &boundary_layer,
                          const ViscousUpdateResult &update_result);

} // namespace boundary_layer_iteration_ops
