#pragma once

#include "solver/boundary_layer/boundary_layer_state.hpp"
#include "solver/boundary_layer/viscous_types.hpp"

class XFoil;

class BoundaryLayerAerodynamicsOps {
 public:
  static BoundaryLayerEdgeVelocityDistribution computeNewUeDistribution(
      const SidePair<BoundaryLayerLattice> &lattice,
      const XFoil &xfoil,
      const BoundaryLayerMatrix3x2dVector &vdel);
  static BoundaryLayerClContributions computeClFromEdgeVelocityDistribution(
      const SidePair<BoundaryLayerLattice> &lattice,
      const XFoil &xfoil,
      const BoundaryLayerEdgeVelocityDistribution &distribution);
  static SidePair<Eigen::VectorXd> ueset(
      const SidePair<BoundaryLayerLattice> &lattice, const Eigen::MatrixXd &dij);
};
