#pragma once

#include "simulation/BoundaryLayer.hpp"
#include "simulation/boundary_layer_solver_ops.hpp"

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
      const BoundaryLayerWorkflow::MixedModeStationContext &ctx);
  double fallbackEdgeVelocity(
      int side, int stationIndex,
      BoundaryLayerWorkflow::EdgeVelocityFallbackMode edgeMode) const;
  void syncStationRegimeStates(int side, int stationIndex,
                               FlowRegimeEnum stationRegime);
  FlowRegimeEnum determineRegimeForStation(int side, int stationIndex) const;
  void updateSystemMatricesForStation(
      const Edge &edge, int side, int stationIndex,
      BoundaryLayerWorkflow::MixedModeStationContext &ctx);
  void initializeFirstIterationState(
      int side, int stationIndex, int previousTransition,
      BoundaryLayerWorkflow::MixedModeStationContext &ctx, double &ueref,
      double &hkref);
  void configureSimilarityRow(double ueref);
  void configureViscousRow(double hkref, double ueref, double senswt,
                           bool resetSensitivity, bool averageSensitivity,
                           double &sens, double &sennew);
  bool applyMixedModeNewtonStep(
      int side, int stationIndex, double &ami,
      BoundaryLayerWorkflow::MixedModeStationContext &ctx);
  void checkTransitionIfNeeded(int side, int stationIndex, bool skipCheck,
                               int laminarAdvance, double &ami);
  void resetStationKinematicsAfterFailure(
      int side, int stationIndex,
      BoundaryLayerWorkflow::MixedModeStationContext &ctx,
      BoundaryLayerWorkflow::EdgeVelocityFallbackMode edgeMode);
  void recoverStationAfterFailure(
      int side, int stationIndex,
      BoundaryLayerWorkflow::MixedModeStationContext &ctx, double &ami,
      BoundaryLayerWorkflow::EdgeVelocityFallbackMode edgeMode,
      int laminarAdvance);

 private:
  Context context_;
};
