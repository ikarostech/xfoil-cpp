#pragma once

#include "simulation/BoundaryLayer.hpp"

class BoundaryLayerRelaxationOps {
 public:
  explicit BoundaryLayerRelaxationOps(BoundaryLayerWorkflow &workflow)
      : workflow_(workflow) {}

  BoundaryLayerWorkflow::BoundaryLayerDelta buildBoundaryLayerDelta(
      int side, const Eigen::VectorXd &unew_side,
      const Eigen::VectorXd &u_ac_side, double dac,
      const BoundaryLayerWorkflow::Matrix3x2dVector &vdel) const;
  BoundaryLayerWorkflow::BoundaryLayerMetrics evaluateSegmentRelaxation(
      int side, const BoundaryLayerWorkflow::BoundaryLayerDelta &delta,
      double dhi, double dlo, double &relaxation) const;
  BoundaryLayerSideProfiles applyBoundaryLayerDelta(
      int side, const BoundaryLayerWorkflow::BoundaryLayerDelta &delta,
      double relaxation, double hstinv, double gamm1) const;

 private:
  BoundaryLayerWorkflow &workflow_;
};
