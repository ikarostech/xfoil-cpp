#pragma once

#include "Eigen/Core"

#include "application/xfoil/XFoilInternalState.hpp"
#include "application/xfoil/XFoilSharedTypes.hpp"
#include "model/coefficient/aero_coefficients.hpp"
#include "model/flow_state.hpp"
#include "model/foil/foil.hpp"
#include "solver/boundary_layer/boundary_layer.hpp"

namespace xfoil_postprocess {

XFoilCompressibilityParams buildCompressibilityParams(
    const FlowState &analysis_state);

XFoilPressureCoefficientResult computePressureCoefficient(
    double tangential_velocity, double velocity_derivative, double qinf,
    const XFoilCompressibilityParams &params);

XFoilClComputation computeCl(const Foil &foil, const FlowState &analysis_state,
                             const XFoilInternalState::InviscidFieldState &state,
                             const Eigen::Vector2d &ref);

double computeCd(const FlowState &analysis_state,
                 const BoundaryLayer &boundary_layer, bool bl_initialized);

void applyClComputation(AeroCoefficients &aero_coeffs,
                        const XFoilClComputation &result);

} // namespace xfoil_postprocess
