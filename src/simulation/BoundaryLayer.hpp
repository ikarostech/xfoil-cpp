#pragma once

#include <cmath>
#include <sstream>

#include "simulation/boundary_layer_state.hpp"
#include "domain/coefficient/bl_newton.hpp"
#include "domain/boundary_layer/boundary_layer_variables_solver.hpp"
#include "simulation/skin_friction_coefficients.hpp"

class XFoil;
enum class FlowRegimeEnum;

class BoundaryLayerWorkflow {
 public:
  BoundaryLayerVariablesSolver boundaryLayerVariablesSolver;
  enum class EdgeVelocityFallbackMode {
    UsePreviousStation,
    AverageNeighbors
  };
  struct EdgeVelocityDistribution {
    SidePair<Eigen::VectorXd> unew;
    SidePair<Eigen::VectorXd> u_ac;
  };
  struct QtanResult {
    Eigen::VectorXd qnew;
    Eigen::VectorXd q_ac;
  };
  struct ClContributions {
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
  struct MixedModeStationContext {
    bool simi = false;
    bool wake = false;
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
  };
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

  bool isStartOfWake(const XFoil& xfoil, int side, int stationIndex);
  void updateSystemMatricesForStation(XFoil& xfoil, int side,
                                      int stationIndex,
                                      MixedModeStationContext& ctx);
  void initializeFirstIterationState(XFoil& xfoil, int side,
                                     int stationIndex,
                                     int previousTransition,
                                     MixedModeStationContext& ctx,
                                     double& ueref, double& hkref,
                                     double& ami);
  void configureSimilarityRow(double ueref);
  void configureViscousRow(double hkref, double ueref,
                           double senswt, bool resetSensitivity,
                           bool averageSensitivity, double& sens,
                           double& sennew);
  bool applyMixedModeNewtonStep(XFoil& xfoil, int side, int stationIndex,
                                double deps, double& ami,
                                MixedModeStationContext& ctx);
  blData blvar(blData data, FlowRegimeEnum flowRegimeType);
  SkinFrictionCoefficients blmid(XFoil& xfoil,
                                 FlowRegimeEnum flowRegimeType);
  blData blprv(XFoil& xfoil, blData data, double xsi, double ami, double cti,
               double thi, double dsi, double dswaki, double uei) const;
  bool blsys(XFoil& xfoil);
  bool trdif(XFoil& xfoil);
  bool mrchdu(XFoil& xfoil);
  bool mrchdu(BoundaryLayerState& state, SidePair<BoundaryLayerLattice>& lattice,
              XFoil& xfoil);
  int resetSideState(int side, XFoil& xfoil);
  bool marchBoundaryLayerSide(BoundaryLayerState& state, int side,
                              double deps, double senswt, double& sens,
                              double& sennew, double& ami, XFoil& xfoil);
  bool processBoundaryLayerStation(BoundaryLayerState& state, int side,
                                   int stationIndex, int previousTransition,
                                   double deps, double senswt, double& sens,
                                   double& sennew, double& ami, XFoil& xfoil);
  bool mrchue(XFoil& xfoil);
  bool mrchue(BoundaryLayerState& state, SidePair<BoundaryLayerLattice>& lattice,
              XFoil& xfoil);
  bool marchMrchueSide(BoundaryLayerState& state, int side,
                       XFoil& xfoil, std::stringstream& ss);
  void initializeMrchueSide(int side, double& thi, double& dsi,
                            double& ami, double& cti, XFoil& xfoil);
  void prepareMrchueStationContext(int side, int stationIndex,
                                   MrchueStationContext& ctx, double thi,
                                   double dsi, double ami, double cti,
                                   XFoil& xfoil);
  bool performMrchueNewtonLoop(int side, int stationIndex,
                               MrchueStationContext& ctx, XFoil& xfoil,
                               std::stringstream& ss);
  void handleMrchueStationFailure(int side, int stationIndex,
                                  MrchueStationContext& ctx, XFoil& xfoil,
                                  std::stringstream& ss);
  void storeMrchueStationState(int side, int stationIndex,
                               const MrchueStationContext& ctx, XFoil& xfoil);
  void storeStationStateCommon(int side, int stationIndex, double ami,
                               double cti, double thi, double dsi, double uei,
                               double xsi, double dswaki, XFoil& xfoil);
  template <typename StationContext>
  void resetStationKinematicsAfterFailure(int side, int stationIndex,
                                          StationContext& ctx,
                                          EdgeVelocityFallbackMode edgeMode);
  double fallbackEdgeVelocity(int side, int stationIndex,
                              EdgeVelocityFallbackMode edgeMode) const;
  EdgeVelocityDistribution computeNewUeDistribution(const XFoil& xfoil) const;
  QtanResult computeQtan(const EdgeVelocityDistribution& distribution) const;
  ClContributions computeClFromEdgeVelocityDistribution(
      const XFoil& xfoil, const EdgeVelocityDistribution& distribution) const;
  BoundaryLayerDelta buildBoundaryLayerDelta(
      int side, const Eigen::VectorXd& unew_side,
      const Eigen::VectorXd& u_ac_side, double dac,
      const XFoil& xfoil) const;
  BoundaryLayerMetrics evaluateSegmentRelaxation(
      int side, const BoundaryLayerDelta& delta, double dhi, double dlo,
      double& relaxation) const;
  BoundaryLayerSideProfiles applyBoundaryLayerDelta(
      int side, const BoundaryLayerDelta& delta, double relaxation,
      XFoil& xfoil) const;
  void syncStationRegimeStates(int side, int stationIndex, bool wake,
                               XFoil& xfoil);
  FlowRegimeEnum determineRegimeForStation(int side, int stationIndex,
                                           bool similarity,
                                           bool wake) const;

  bool iblpan(XFoil& xfoil);
  bool iblsys(XFoil& xfoil);
  bool stfind(XFoil& xfoil);
  bool stmove(XFoil& xfoil);
  bool uicalc(XFoil& xfoil);
  bool tesys(XFoil& xfoil, double cte, double tte, double dte);
  bool trchek(XFoil& xfoil);

private:
  void copyStationState(int side, int destination, int source);

public:
  SidePair<BoundaryLayerLattice> lattice;
  BlSystemCoeffs blc;
  BoundaryLayerState state;
  int stagnationIndex = 0;
};

template <typename StationContext>
inline void BoundaryLayerWorkflow::resetStationKinematicsAfterFailure(
    int side, int stationIndex, StationContext& ctx,
    EdgeVelocityFallbackMode edgeMode) {
  if (ctx.dmax <= 0.1 || stationIndex < 2) {
    return;
  }

  if (stationIndex <= lattice.get(side).trailingEdgeIndex) {
    const double ratio =
        lattice.get(side).arcLengthCoordinates[stationIndex] /
        lattice.get(side).arcLengthCoordinates[stationIndex - 1];
    const double scale = std::sqrt(ratio);
    ctx.thi = lattice.get(side).profiles.momentumThickness[stationIndex - 1] * scale;
    ctx.dsi = lattice.get(side).profiles.displacementThickness[stationIndex - 1] * scale;
  } else {
    if (stationIndex == lattice.get(side).trailingEdgeIndex + 1) {
      ctx.cti = ctx.cte;
      ctx.thi = ctx.tte;
      ctx.dsi = ctx.dte;
    } else {
      ctx.thi = lattice.get(side).profiles.momentumThickness[stationIndex - 1];
      const double ratlen =
          (lattice.get(side).arcLengthCoordinates[stationIndex] -
           lattice.get(side).arcLengthCoordinates[stationIndex - 1]) /
          (10.0 * lattice.get(side).profiles.displacementThickness[stationIndex - 1]);
      ctx.dsi =
          (lattice.get(side).profiles.displacementThickness[stationIndex - 1] + ctx.thi * ratlen) /
          (1.0 + ratlen);
    }
  }

  ctx.uei = fallbackEdgeVelocity(side, stationIndex, edgeMode);

  if (stationIndex == lattice.get(side).transitionIndex) {
    ctx.cti = 0.05;
  }
  if (stationIndex > lattice.get(side).transitionIndex) {
    ctx.cti = lattice.get(side).profiles.skinFrictionCoeff[stationIndex - 1];
  }
}
