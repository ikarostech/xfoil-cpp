#pragma once

#include <sstream>
#include <vector>

#include "Eigen/Core"
#include "simulation/BoundaryLayer.hpp"

class XFoil;

class BoundaryLayerMarcher {
 public:
  using Matrix3x2d = Eigen::Matrix<double, 3, 2>;
  using Matrix3x2dVector = std::vector<Matrix3x2d>;

  struct MrchueStationContext {
    bool simi = false;
    bool wake = false;
    bool direct = true;
    double xsi = 0.0;
    double uei = 0.0;
    double thi = 0.0;
    double dsi = 0.0;
    double ami = 0.0;
    double cti = 0.0;
    double dswaki = 0.0;
    double cte = 0.0;
    double dte = 0.0;
    double tte = 0.0;
    double dmax = 0.0;
    double hmax = 0.0;
    double htarg = 0.0;
  };

  bool mrchdu(BoundaryLayerWorkflow& workflow, const Foil& foil,
              const StagnationResult& stagnation);
  bool mrchdu(BoundaryLayerWorkflow& workflow, BoundaryLayerState& state,
              const Foil& foil, const StagnationResult& stagnation);
  int resetSideState(BoundaryLayerWorkflow& workflow, int side, const Foil& foil,
                     const StagnationResult& stagnation);
  bool marchBoundaryLayerSide(BoundaryLayerWorkflow& workflow,
                              BoundaryLayerState& state, int side,
                              double deps, double senswt, double& sens,
                              double& sennew, double& ami, const Foil& foil,
                              const StagnationResult& stagnation);
  bool processBoundaryLayerStation(BoundaryLayerWorkflow& workflow,
                                   BoundaryLayerState& state, int side,
                                   int stationIndex, int previousTransition,
                                   double deps, double senswt, double& sens,
                                   double& sennew, double& ami,
                                   const Foil& foil);
  bool performMixedModeNewtonIteration(
      BoundaryLayerWorkflow& workflow, const Edge& edge, int side, int ibl,
      int itrold, BoundaryLayerWorkflow::MixedModeStationContext& ctx,
      double deps, double senswt, double& sens, double& sennew, double& ami);
  void handleMixedModeNonConvergence(
      BoundaryLayerWorkflow& workflow, int side, int ibl,
      BoundaryLayerWorkflow::MixedModeStationContext& ctx, double& ami);
  bool mrchue(BoundaryLayerWorkflow& workflow, const Foil& foil,
              const StagnationResult& stagnation);
  bool mrchue(BoundaryLayerWorkflow& workflow, BoundaryLayerState& state,
              const Foil& foil, const StagnationResult& stagnation);
  bool marchMrchueSide(BoundaryLayerWorkflow& workflow,
                       BoundaryLayerState& state, int side,
                       const Foil& foil, const StagnationResult& stagnation,
                       std::stringstream& ss);
  void initializeMrchueSide(BoundaryLayerWorkflow& workflow, int side,
                            double& thi, double& dsi, double& ami,
                            double& cti);
  void prepareMrchueStationContext(
      BoundaryLayerWorkflow& workflow, int side, int stationIndex,
      BoundaryLayerMarcher::MrchueStationContext& ctx, double thi,
      double dsi, double ami, double cti);
  bool performMrchueNewtonLoop(
      BoundaryLayerWorkflow& workflow, int side, int stationIndex,
      BoundaryLayerMarcher::MrchueStationContext& ctx, const Edge& edge,
      std::stringstream& ss);
  void handleMrchueStationFailure(
      BoundaryLayerWorkflow& workflow, int side, int stationIndex,
      BoundaryLayerMarcher::MrchueStationContext& ctx,
      std::stringstream& ss);
  void storeMrchueStationState(
      BoundaryLayerWorkflow& workflow, int side, int stationIndex,
      const BoundaryLayerMarcher::MrchueStationContext& ctx);
  void storeStationStateCommon(BoundaryLayerWorkflow& workflow, int side,
                               int stationIndex, double ami, double cti,
                               double thi, double dsi, double uei, double xsi,
                               double dswaki);
  template <typename StationContext>
  void resetStationKinematicsAfterFailure(
      BoundaryLayerWorkflow& workflow, int side, int stationIndex,
      StationContext& ctx,
      BoundaryLayerWorkflow::EdgeVelocityFallbackMode edgeMode);
  double fallbackEdgeVelocity(
      const BoundaryLayerWorkflow& workflow, int side, int stationIndex,
      BoundaryLayerWorkflow::EdgeVelocityFallbackMode edgeMode) const;
  BoundaryLayerWorkflow::EdgeVelocityDistribution computeNewUeDistribution(
      const BoundaryLayerWorkflow& workflow,
      const XFoil& xfoil,
      const Matrix3x2dVector& vdel) const;
  BoundaryLayerWorkflow::QtanResult computeQtan(
      const BoundaryLayerWorkflow& workflow,
      const BoundaryLayerWorkflow::EdgeVelocityDistribution& distribution,
      int point_count) const;
  BoundaryLayerWorkflow::ClContributions
  computeClFromEdgeVelocityDistribution(
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
      const BoundaryLayerWorkflow::BoundaryLayerDelta& delta,
      double relaxation, double hstinv, double gamm1) const;
  void syncStationRegimeStates(BoundaryLayerWorkflow& workflow, int side,
                               int stationIndex, bool wake);
  FlowRegimeEnum determineRegimeForStation(const BoundaryLayerWorkflow& workflow,
                                           int side, int stationIndex,
                                           bool similarity, bool wake) const;
  BoundaryLayerWorkflow::MixedModeStationContext prepareMixedModeStation(
      BoundaryLayerWorkflow& workflow, int side, int stationIndex,
      int previousTransition, double& ami);
};

template <typename StationContext>
inline void BoundaryLayerMarcher::resetStationKinematicsAfterFailure(
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
    ctx.thi = workflow.lattice.get(side).profiles.momentumThickness[stationIndex - 1] * scale;
    ctx.dsi = workflow.lattice.get(side).profiles.displacementThickness[stationIndex - 1] * scale;
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
