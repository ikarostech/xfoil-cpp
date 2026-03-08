#pragma once

#include "simulation/BoundaryLayer.hpp"

class BoundaryLayerMixedModeOps {
 public:
  explicit BoundaryLayerMixedModeOps(BoundaryLayerWorkflow &workflow)
      : workflow_(workflow) {}

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
  BoundaryLayerWorkflow &workflow_;
};
