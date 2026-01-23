#pragma once
#include "core/side_pair.hpp"
#include "Eigen/Core"
#include "infrastructure/xfoil_params.h"
#include "domain/flow_regime.hpp"
#include "XFoil.h"

struct SetblInputView {
  const bool& lblini;
  SidePairRef<const BoundaryLayerSideProfiles> profiles;

  static SetblInputView fromXFoil(const XFoil& xfoil);
};

struct SetblOutputView {
  bool& lblini;
  XFoil::BlCompressibilityParams& blCompressibility;
  XFoil::BlReynoldsParams& blReynolds;
  XFoil::BlTransitionParams& blTransition;
  SidePairRef<BoundaryLayerSideProfiles> profiles;
  std::vector<Eigen::Matrix<double, 3, 2>>& va;
  std::vector<Eigen::Matrix<double, 3, 2>>& vb;
  std::vector<Eigen::Matrix<double, 3, 2>>& vdel;
  double (&vm)[3][IZX][IZX];
  double (&vz)[3][2];
  FlowRegimeEnum& flowRegime;

  static SetblOutputView fromXFoil(XFoil& xfoil);
};
