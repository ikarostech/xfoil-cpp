#pragma once

#include "model/boundary_layer/state.hpp"
#include "solver/boundary_layer/viscous_types.hpp"

class BoundaryLayerRelaxationOps {
 public:
  static BoundaryLayerDelta buildBoundaryLayerDelta(
      const BoundaryLayerLattice &lattice_side, const Eigen::VectorXd &unew_side,
      const Eigen::VectorXd &u_ac_side, double dac,
      const BoundaryLayerMatrix3x2dVector &vdel);
  static BoundaryLayerMetrics evaluateSegmentRelaxation(
      const BoundaryLayerSideProfiles &profiles,
      const BoundaryLayerDelta &delta, double dhi,
      double dlo, double &relaxation);
  static BoundaryLayerSideProfiles applyBoundaryLayerDelta(
      const BoundaryLayerLattice &lattice_side, const Eigen::VectorXd &wgap,
      const BoundaryLayerDelta &delta, double relaxation,
      double hstinv, double gamm1);
};
