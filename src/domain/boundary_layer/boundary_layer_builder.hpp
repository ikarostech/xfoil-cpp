#pragma once
#include "core/side_pair.hpp"
#include "Eigen/Core"
#include "domain/flow_regime.hpp"
#include "XFoil.h"

struct SetblOutputView {
  bool lblini = false;
  XFoil::BlCompressibilityParams blCompressibility{};
  XFoil::BlReynoldsParams blReynolds{};
  XFoil::BlTransitionParams blTransition{};
  SidePair<BoundaryLayerSideProfiles> profiles{};
  std::vector<Eigen::Matrix<double, 3, 2>> va;
  std::vector<Eigen::Matrix<double, 3, 2>> vb;
  std::vector<Eigen::Matrix<double, 3, 2>> vdel;
  XFoil::VmMatrix vm;
  XFoil::VzMatrix vz{};
  FlowRegimeEnum flowRegime = FlowRegimeEnum::Laminar;

  static SetblOutputView fromXFoil(const XFoil& xfoil);
  void applyToXFoil(XFoil& xfoil);
};
