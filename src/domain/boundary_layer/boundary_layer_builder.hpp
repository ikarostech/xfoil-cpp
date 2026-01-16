#pragma once
#include "core/side_pair.hpp"
#include "Eigen/Core"
#include "infrastructure/xfoil_params.h"
#include "domain/flow_regime.hpp"
#include "XFoil.h"

struct SetblInputView {
  const bool& lblini;
  SidePairRef<const Eigen::VectorXd> edgeVelocity;
  SidePairRef<const Eigen::VectorXd> skinFrictionCoeff;
  SidePairRef<const Eigen::VectorXd> momentumThickness;
  SidePairRef<const Eigen::VectorXd> displacementThickness;
  SidePairRef<const Eigen::VectorXd> massFlux;
  SidePairRef<const Eigen::VectorXd> skinFrictionCoeffHistory;
  SidePairRef<const int> itran;

  static SetblInputView fromXFoil(const XFoil& xfoil);
};

struct SetblOutputView {
  bool& lblini;
  XFoil::BlCompressibilityParams& blCompressibility;
  XFoil::BlReynoldsParams& blReynolds;
  XFoil::BlTransitionParams& blTransition;
  SidePairRef<Eigen::VectorXd> edgeVelocity;
  SidePairRef<Eigen::VectorXd> skinFrictionCoeff;
  SidePairRef<Eigen::VectorXd> momentumThickness;
  SidePairRef<Eigen::VectorXd> displacementThickness;
  SidePairRef<Eigen::VectorXd> massFlux;
  SidePairRef<Eigen::VectorXd> skinFrictionCoeffHistory;
  SidePairRef<int> itran;
  std::vector<Eigen::Matrix<double, 3, 2>>& va;
  std::vector<Eigen::Matrix<double, 3, 2>>& vb;
  std::vector<Eigen::Matrix<double, 3, 2>>& vdel;
  double (&vm)[3][IZX][IZX];
  double (&vz)[3][2];
  FlowRegimeEnum& flowRegime;

  static SetblOutputView fromXFoil(XFoil& xfoil);
};
