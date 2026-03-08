#pragma once

#include "simulation/BoundaryLayer.hpp"

class BoundaryLayerRelaxationOps {
 public:
  static BoundaryLayerWorkflow::BoundaryLayerDelta buildBoundaryLayerDelta(
      const BoundaryLayerLattice &lattice_side, const Eigen::VectorXd &unew_side,
      const Eigen::VectorXd &u_ac_side, double dac,
      const BoundaryLayerWorkflow::Matrix3x2dVector &vdel);
  static BoundaryLayerWorkflow::BoundaryLayerMetrics evaluateSegmentRelaxation(
      const BoundaryLayerSideProfiles &profiles,
      const BoundaryLayerWorkflow::BoundaryLayerDelta &delta, double dhi,
      double dlo, double &relaxation);
  static BoundaryLayerSideProfiles applyBoundaryLayerDelta(
      const BoundaryLayerLattice &lattice_side, const Eigen::VectorXd &wgap,
      const BoundaryLayerWorkflow::BoundaryLayerDelta &delta, double relaxation,
      double hstinv, double gamm1);
};
