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
  std::vector<Eigen::Matrix<double, 3, 2>> va;
  std::vector<Eigen::Matrix<double, 3, 2>> vb;
  std::vector<Eigen::Matrix<double, 3, 2>> vdel;
  XFoil::VmMatrix vm;
  XFoil::VzMatrix vz{};
  FlowRegimeEnum flowRegime = FlowRegimeEnum::Laminar;

  void applyToXFoil(XFoil& xfoil);
};
