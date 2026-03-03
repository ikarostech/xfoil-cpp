#pragma once

#include <vector>

#include "Eigen/Core"
#include "simulation/BoundaryLayer.hpp"

class XFoil;

class Marcher {
 public:
  using Matrix3x2d = Eigen::Matrix<double, 3, 2>;
  using Matrix3x2dVector = std::vector<Matrix3x2d>;

  int resetSideState(BoundaryLayerWorkflow& workflow, int side, const Foil& foil,
                     const StagnationResult& stagnation);
  void storeStationStateCommon(BoundaryLayerWorkflow& workflow, int side,
                               int stationIndex, double ami, double cti,
                               double thi, double dsi, double uei, double xsi,
                               double dswaki);
  void syncStationRegimeStates(BoundaryLayerWorkflow& workflow, int side,
                               int stationIndex, bool wake);
  FlowRegimeEnum determineRegimeForStation(const BoundaryLayerWorkflow& workflow,
                                           int side, int stationIndex,
                                           bool similarity, bool wake) const;
  double fallbackEdgeVelocity(
      const BoundaryLayerWorkflow& workflow, int side, int stationIndex,
      BoundaryLayerWorkflow::EdgeVelocityFallbackMode edgeMode) const;
  template <typename StationContext>
  void resetStationKinematicsAfterFailure(
      BoundaryLayerWorkflow& workflow, int side, int stationIndex,
      StationContext& ctx,
      BoundaryLayerWorkflow::EdgeVelocityFallbackMode edgeMode);

  BoundaryLayerWorkflow::EdgeVelocityDistribution computeNewUeDistribution(
      const BoundaryLayerWorkflow& workflow, const XFoil& xfoil,
      const Matrix3x2dVector& vdel) const;
  BoundaryLayerWorkflow::QtanResult computeQtan(
      const BoundaryLayerWorkflow& workflow,
      const BoundaryLayerWorkflow::EdgeVelocityDistribution& distribution,
      int point_count) const;
  BoundaryLayerWorkflow::ClContributions computeClFromEdgeVelocityDistribution(
      const BoundaryLayerWorkflow& workflow, const XFoil& xfoil,
      const BoundaryLayerWorkflow::EdgeVelocityDistribution& distribution) const;
  BoundaryLayerWorkflow::BoundaryLayerDelta buildBoundaryLayerDelta(
      const BoundaryLayerWorkflow& workflow, int side,
      const Eigen::VectorXd& unew_side, const Eigen::VectorXd& u_ac_side,
      double dac, const XFoil& xfoil, const Matrix3x2dVector& vdel) const;
  BoundaryLayerWorkflow::BoundaryLayerMetrics evaluateSegmentRelaxation(
      const BoundaryLayerWorkflow& workflow, int side,
      const BoundaryLayerWorkflow::BoundaryLayerDelta& delta, double dhi,
      double dlo, double& relaxation) const;
  BoundaryLayerSideProfiles applyBoundaryLayerDelta(
      const BoundaryLayerWorkflow& workflow, int side,
      const BoundaryLayerWorkflow::BoundaryLayerDelta& delta, double relaxation,
      double hstinv, double gamm1) const;
};

template <typename StationContext>
inline void Marcher::resetStationKinematicsAfterFailure(
    BoundaryLayerWorkflow& workflow, int side, int stationIndex,
    StationContext& ctx,
    BoundaryLayerWorkflow::EdgeVelocityFallbackMode edgeMode) {
  if (ctx.dmax <= 0.1 || stationIndex < 2) {
    return;
  }

  if (stationIndex <= workflow.lattice.get(side).trailingEdgeIndex) {
    const double ratio =
        workflow.lattice.get(side).arcLengthCoordinates[stationIndex] /
        workflow.lattice.get(side).arcLengthCoordinates[stationIndex - 1];
    const double scale = std::sqrt(ratio);
    ctx.thi =
        workflow.lattice.get(side).profiles.momentumThickness[stationIndex - 1] *
        scale;
    ctx.dsi =
        workflow.lattice.get(side).profiles.displacementThickness[stationIndex - 1] *
        scale;
  } else {
    if (stationIndex == workflow.lattice.get(side).trailingEdgeIndex + 1) {
      ctx.cti = ctx.cte;
      ctx.thi = ctx.tte;
      ctx.dsi = ctx.dte;
    } else {
      ctx.thi = workflow.lattice.get(side).profiles.momentumThickness[stationIndex - 1];
      const double ratlen =
          (workflow.lattice.get(side).arcLengthCoordinates[stationIndex] -
           workflow.lattice.get(side).arcLengthCoordinates[stationIndex - 1]) /
          (10.0 * workflow.lattice.get(side).profiles.displacementThickness[stationIndex - 1]);
      ctx.dsi =
          (workflow.lattice.get(side).profiles.displacementThickness[stationIndex - 1] + ctx.thi * ratlen) /
          (1.0 + ratlen);
    }
  }

  ctx.uei = fallbackEdgeVelocity(workflow, side, stationIndex, edgeMode);

  if (stationIndex == workflow.lattice.get(side).profiles.transitionIndex) {
    ctx.cti = 0.05;
  }
  if (stationIndex > workflow.lattice.get(side).profiles.transitionIndex) {
    ctx.cti = workflow.lattice.get(side).profiles.skinFrictionCoeff[stationIndex - 1];
  }
}
