#pragma once

#include "simulation/boundary_layer_solver_ops.hpp"
#include "simulation/viscous_types.hpp"

class BoundaryLayerMixedModeOps {
 public:
  struct Context {
    SidePair<BoundaryLayerLattice> &lattice;
    BoundaryLayerState &state;
    FlowRegimeEnum &flowRegime;
    BlSystemCoeffs &blc;
    const BlCompressibilityParams &blCompressibility;
    BoundaryLayerVariablesSolver &boundaryLayerVariablesSolver;
    BoundaryLayerTransitionSolver &transitionSolver;
    BoundaryLayerSolverOps solverOps;
  };

  explicit BoundaryLayerMixedModeOps(Context context) : context_(context) {}

  void storeStationStateCommon(
      int side, int stationIndex,
      const BoundaryLayerMixedModeStationContext &ctx);
  double fallbackEdgeVelocity(
      int side, int stationIndex,
      BoundaryLayerEdgeVelocityFallbackMode edgeMode) const;
  void syncStationRegimeStates(int side, int stationIndex,
                               FlowRegimeEnum stationRegime);
  FlowRegimeEnum determineRegimeForStation(int side, int stationIndex) const;
  void updateSystemMatricesForStation(
      const Edge &edge, int side, int stationIndex,
      BoundaryLayerMixedModeStationContext &ctx);
  void initializeFirstIterationState(
      int side, int stationIndex, int previousTransition,
      BoundaryLayerMixedModeStationContext &ctx, double &ueref,
      double &hkref);
  void configureSimilarityRow(double ueref);
  void configureViscousRow(double hkref, double ueref, double senswt,
                           bool resetSensitivity, bool averageSensitivity,
                           double &sens, double &sennew);
  bool applyMixedModeNewtonStep(
      int side, int stationIndex, double &ami,
      BoundaryLayerMixedModeStationContext &ctx);
  void checkTransitionIfNeeded(int side, int stationIndex, bool skipCheck,
                               int laminarAdvance, double &ami);
  void resetStationKinematicsAfterFailure(
      int side, int stationIndex,
      BoundaryLayerMixedModeStationContext &ctx,
      BoundaryLayerEdgeVelocityFallbackMode edgeMode);
  void recoverStationAfterFailure(
      int side, int stationIndex,
      BoundaryLayerMixedModeStationContext &ctx, double &ami,
      BoundaryLayerEdgeVelocityFallbackMode edgeMode,
      int laminarAdvance);

 private:
  Context context_;
};
