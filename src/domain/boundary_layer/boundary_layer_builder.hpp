#pragma once
#include "core/side_pair.hpp"
#include "Eigen/Core"
#include "domain/flow_regime.hpp"
#include "domain/boundary_layer/bl_compressibility_params.hpp"
#include "domain/boundary_layer/bl_reynolds_params.hpp"
#include "domain/boundary_layer/bl_transition_params.hpp"
#include "XFoil.h"

struct SetblOutputView {
  bool lblini = false;
  BlCompressibilityParams blCompressibility{};
  BlReynoldsParams blReynolds{};
  BlTransitionParams blTransition{};
  SidePair<BoundaryLayerSideProfiles> profiles{};
  Blsolve::BlNewtonSystem bl_newton_system{};
  FlowRegimeEnum flowRegime = FlowRegimeEnum::Laminar;

  void applyToXFoil(XFoil& xfoil);
};
