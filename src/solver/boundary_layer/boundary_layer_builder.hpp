#pragma once
#include "Eigen/Core"
#include "numerics/side_pair.hpp"
#include "model/boundary_layer.hpp"
#include "model/boundary_layer/state.hpp"
#include "model/boundary_layer/bl_compressibility_params.hpp"
#include "model/boundary_layer/bl_reynolds_params.hpp"
#include "model/boundary_layer/bl_transition_params.hpp"
#include "model/flow_regime.hpp"
#include "solver/boundary_layer/blsolve.hpp"

struct SetblOutputView {
  BlCompressibilityParams blCompressibility{};
  BlReynoldsParams blReynolds{};
  BlTransitionParams blTransition{};
  SidePair<BoundaryLayerSideProfiles> profiles{};
  Blsolve::BlNewtonSystem bl_newton_system{};
  FlowRegimeEnum flowRegime = FlowRegimeEnum::Laminar;
};
