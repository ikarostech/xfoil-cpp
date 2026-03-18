#pragma once

#include <vector>

#include "Eigen/Core"

#include "model/boundary_layer/runtime_types.hpp"
#include "numerics/side_pair.hpp"
#include "model/flow_regime.hpp"

using BoundaryLayerMatrix3x2d = Eigen::Matrix<double, 3, 2>;
using BoundaryLayerMatrix3x2dVector = std::vector<BoundaryLayerMatrix3x2d>;

enum class BoundaryLayerEdgeVelocityFallbackMode {
  UsePreviousStation,
  AverageNeighbors
};

struct BoundaryLayerEdgeVelocityDistribution {
  SidePair<Eigen::VectorXd> unew;
  SidePair<Eigen::VectorXd> u_ac;
};

struct BoundaryLayerClContributions {
  double cl = 0.0;
  double cl_a = 0.0;
  double cl_ms = 0.0;
  double cl_ac = 0.0;
};

struct BoundaryLayerDelta {
  Eigen::VectorXd dskinFrictionCoeff;
  Eigen::VectorXd dmomentumThickness;
  Eigen::VectorXd ddisplacementThickness;
  Eigen::VectorXd dedgeVelocity;
};

struct BoundaryLayerMetrics {
  double rmsContribution = 0.0;
  double maxChange = 0.0;
};

struct BoundaryLayerMixedModeStationContext {
  FlowRegimeEnum flowRegime = FlowRegimeEnum::Laminar;
  double xsi = 0.0;
  double uei = 0.0;
  double thi = 0.0;
  double dsi = 0.0;
  double cti = 0.0;
  double ami = 0.0;
  double dswaki = 0.0;
  double cte = 0.0;
  double dte = 0.0;
  double tte = 0.0;
  double dmax = 0.0;

  bool isSimilarity() const { return flowRegime == FlowRegimeEnum::Similarity; }
  bool isWake() const { return flowRegime == FlowRegimeEnum::Wake; }
};
