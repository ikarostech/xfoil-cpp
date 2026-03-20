#pragma once

#include <string_view>

#include "infrastructure/logger.hpp"
#include "solver/boundary_layer/runtime/state.hpp"
#include "solver/march/mrchue_linear_system.hpp"
#include "solver/boundary_layer/initialization/setbl_access.hpp"

class BoundaryLayerMarchAccess {
 public:
  explicit BoundaryLayerMarchAccess(BoundaryLayerWorkflow &workflow)
      : workflow_(workflow) {}

  BoundaryLayerState &state() const { return workflow_.state(); }
  Eigen::VectorXd &wgap() const { return workflow_.wakeGap(); }
  BlSystemCoeffs &blc() const { return workflow_.systemCoefficients(); }
  FlowRegimeEnum &flowRegime() const { return workflow_.flowRegime(); }
  BlCompressibilityParams &blCompressibility() const {
    return workflow_.compressibility();
  }
  BlReynoldsParams &blReynolds() const { return workflow_.reynolds(); }
  blDiff &xt() const { return workflow_.transitionSensitivity(); }

  int readSideStationCount(int side) const {
    return workflow_.readSideStationCount(side);
  }
  BoundaryLayerStationReadModel readStationModel(
      int side, int stationIndex) const {
    return workflow_.readStationModel(side, stationIndex);
  }
  BoundaryLayerTrailingEdgeReadModel readTrailingEdgeModel() const {
    return workflow_.readTrailingEdgeModel();
  }

  void emitMarchInfoLog(std::string_view message) const {
    Logger::instance().write(std::string(message));
  }

  double readNewtonRhs(int row) const {
    return workflow_.readNewtonRhs(row);
  }
  void solveMrchueDirectNewtonSystem() const { workflow_.solveDirectNewtonSystem(); }
  void solveMrchueInverseNewtonSystem(double htarg) const {
    workflow_.solveInverseNewtonSystem(htarg);
  }
  void runTransitionCheckForMrchue(int side, int stationIndex, double &ami,
                                   double &cti) const {
    workflow_.runTransitionCheckForMrchue(side, stationIndex, ami, cti);
  }
  double calcHtarg(int stationIndex, int side, bool wake) const {
    return workflow_.calcHtarg(stationIndex, side, wake);
  }

  bool isStartOfWake(int side, int stationIndex) const {
    return workflow_.isStartOfWake(side, stationIndex);
  }
  FlowRegimeEnum applyFlowRegimeCandidate(FlowRegimeEnum candidate) const {
    workflow_.flowRegime() = candidate;
    return workflow_.flowRegime();
  }
  FlowRegimeEnum currentFlowRegime() const { return workflow_.flowRegime(); }
  FlowRegimeEnum determineRegimeForStation(int side, int stationIndex) const {
    return workflow_.determineRegimeForStation(side, stationIndex);
  }
  int resetSideState(int side, const Foil &foil,
                     const StagnationResult &stagnation) const {
    return workflow_.resetSideState(side, foil, stagnation);
  }

  bool tesys(const BoundaryLayerSideProfiles &top_profiles,
             const BoundaryLayerSideProfiles &bottom_profiles,
             const Edge &edge) const {
    return workflow_.tesys(top_profiles, bottom_profiles, edge);
  }
  bool blsys() const { return workflow_.blsys(); }
  bool blkin(BoundaryLayerState &state) const { return workflow_.blkin(state); }
  blData blprv(blData data, double xsi, double ami, double cti, double thi,
               double dsi, double dswaki, double uei) const {
    return workflow_.blprv(data, xsi, ami, cti, thi, dsi, dswaki, uei);
  }

  void updateSystemMatricesForStation(
      const Edge &edge, int side, int stationIndex,
      BoundaryLayerMixedModeStationContext &ctx) const {
    workflow_.updateSystemMatricesForStation(edge, side, stationIndex, ctx);
  }
  void initializeFirstIterationState(
      int side, int stationIndex, int previousTransition,
      BoundaryLayerMixedModeStationContext &ctx, double &ueref,
      double &hkref) const {
    workflow_.initializeFirstIterationState(side, stationIndex,
                                            previousTransition, ctx, ueref,
                                            hkref);
  }
  void configureSimilarityRow(double ueref) const { workflow_.configureSimilarityRow(ueref); }
  void configureViscousRow(double hkref, double ueref, double senswt,
                           bool resetSensitivity, bool averageSensitivity,
                           double &sens, double &sennew) const {
    workflow_.configureViscousRow(hkref, ueref, senswt, resetSensitivity,
                                  averageSensitivity, sens, sennew);
  }
  bool applyMixedModeNewtonStep(
      int side, int stationIndex, double &ami,
      BoundaryLayerMixedModeStationContext &ctx) const {
    return workflow_.applyMixedModeNewtonStep(side, stationIndex, ami, ctx);
  }
  void checkTransitionIfNeeded(int side, int stationIndex, bool skipCheck,
                               int laminarAdvance, double &ami) const {
    workflow_.checkTransitionIfNeeded(side, stationIndex, skipCheck,
                                      laminarAdvance, ami);
  }
  void storeStationStateCommon(
      int side, int stationIndex,
      const BoundaryLayerMixedModeStationContext &ctx) const {
    workflow_.storeStationStateCommon(side, stationIndex, ctx);
  }
  void recoverStationAfterFailure(
      int side, int stationIndex,
      BoundaryLayerMixedModeStationContext &ctx, double &ami,
      BoundaryLayerEdgeVelocityFallbackMode edgeMode,
      int laminarAdvance) const {
    workflow_.recoverStationAfterFailure(side, stationIndex, ctx, ami, edgeMode,
                                         laminarAdvance);
  }
  bool solveTeSystemForCurrentProfiles(const Edge &edge) const {
    return workflow_.solveTeSystemForCurrentProfiles(edge);
  }

 private:
  BoundaryLayerWorkflow &workflow_;
};

inline BoundaryLayerMarchAccess makeBoundaryLayerMarchAccess(
    BoundaryLayerWorkflow &workflow) {
  return BoundaryLayerMarchAccess(workflow);
}
