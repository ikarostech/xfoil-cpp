#pragma once

#include "model/boundary_layer/state.hpp"
#include "solver/boundary_layer/viscous_types.hpp"

struct BoundaryLayerAerodynamicContext {
  const Eigen::MatrixXd &dij;
  const Eigen::Matrix2Xd &foilPoints;
  double alpha = 0.0;
  double qinf = 0.0;
  double currentMach = 0.0;
  bool controlByAlpha = true;
};

class BoundaryLayerAerodynamicsOps {
 public:
  static BoundaryLayerEdgeVelocityDistribution computeNewUeDistribution(
      const SidePair<BoundaryLayerLattice> &lattice,
      const BoundaryLayerAerodynamicContext &context,
      const BoundaryLayerMatrix3x2dVector &vdel);
  static BoundaryLayerClContributions computeClFromEdgeVelocityDistribution(
      const SidePair<BoundaryLayerLattice> &lattice,
      const BoundaryLayerAerodynamicContext &context,
      const BoundaryLayerEdgeVelocityDistribution &distribution);
  static SidePair<Eigen::VectorXd> ueset(
      const SidePair<BoundaryLayerLattice> &lattice, const Eigen::MatrixXd &dij);
};
