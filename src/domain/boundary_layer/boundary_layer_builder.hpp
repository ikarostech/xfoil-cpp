#pragma once
#include "core/side_pair.hpp"
#include "Eigen/Core"
#include "domain/flow_regime.hpp"
#include "XFoil.h"

struct SetblOutputView {
  bool& lblini;
  XFoil::BlCompressibilityParams& blCompressibility;
  XFoil::BlReynoldsParams& blReynolds;
  XFoil::BlTransitionParams& blTransition;
  SidePairRef<BoundaryLayerSideProfiles> profiles;
  std::vector<Eigen::Matrix<double, 3, 2>>& va;
  std::vector<Eigen::Matrix<double, 3, 2>>& vb;
  std::vector<Eigen::Matrix<double, 3, 2>>& vdel;
  XFoil::VmMatrix& vm;
  XFoil::VzMatrix& vz;
  FlowRegimeEnum& flowRegime;

  static SetblOutputView fromXFoil(XFoil& xfoil);
};
