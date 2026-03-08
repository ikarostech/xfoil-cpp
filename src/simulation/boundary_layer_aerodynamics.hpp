#pragma once

#include "simulation/BoundaryLayer.hpp"

class XFoil;

class BoundaryLayerAerodynamicsOps {
 public:
  static BoundaryLayerWorkflow::EdgeVelocityDistribution computeNewUeDistribution(
      const SidePair<BoundaryLayerLattice> &lattice,
      const XFoil &xfoil,
      const BoundaryLayerWorkflow::Matrix3x2dVector &vdel);
  static BoundaryLayerWorkflow::ClContributions computeClFromEdgeVelocityDistribution(
      const SidePair<BoundaryLayerLattice> &lattice,
      const XFoil &xfoil,
      const BoundaryLayerWorkflow::EdgeVelocityDistribution &distribution);
  static SidePair<Eigen::VectorXd> ueset(
      const SidePair<BoundaryLayerLattice> &lattice, const Eigen::MatrixXd &dij);
};
