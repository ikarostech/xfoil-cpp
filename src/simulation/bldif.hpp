#pragma once

#include "domain/coefficient/bl_newton.hpp"
#include "domain/flow_regime.hpp"
#include "simulation/boundary_layer_state.hpp"
#include "simulation/skin_friction_coefficients.hpp"

BlSystemCoeffs AssembleBoundaryLayerSystem(FlowRegimeEnum flowRegimeType,
                                           BoundaryLayerState boundaryLayerState,
                                           const SkinFrictionCoefficients& skinFriction,
                                           double amcrit);
