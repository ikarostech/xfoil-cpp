#pragma once

#include "solver/boundary_layer/boundary_layer_workflow.hpp"
#include "solver/boundary_layer/boundary_layer_aerodynamics.hpp"

class XFoil;

class BoundaryLayerAerodynamicCoupling {
 public:
  static BoundaryLayerEdgeVelocityDistribution computeNewUeDistribution(
      const BoundaryLayerWorkflow &workflow, const XFoil &xfoil,
      const BoundaryLayerMatrix3x2dVector &vdel) {
    return BoundaryLayerAerodynamicsOps::computeNewUeDistribution(
        workflow.stateStore().lattice, xfoil, vdel);
  }

  static BoundaryLayerClContributions computeClFromEdgeVelocityDistribution(
      const BoundaryLayerWorkflow &workflow, const XFoil &xfoil,
      const BoundaryLayerEdgeVelocityDistribution &distribution) {
    return BoundaryLayerAerodynamicsOps::computeClFromEdgeVelocityDistribution(
        workflow.stateStore().lattice, xfoil, distribution);
  }
};
