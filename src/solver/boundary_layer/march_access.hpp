#pragma once

#include <string_view>

#include "infrastructure/logger.hpp"
#include "solver/boundary_layer/runtime/state.hpp"
#include "solver/march/mrchue_linear_system.hpp"
#include "solver/boundary_layer/initialization/setbl_access.hpp"

class BoundaryLayerMarchAccess {
 public:
  BoundaryLayerMarchAccess(BoundaryLayerWorkflow &workflow,
                           BoundaryLayerSolverOps solver_ops)
      : workflow_(workflow), solverOps_(solver_ops) {}

  BoundaryLayerState &state() const { return workflow_.state(); }
  SidePair<BoundaryLayerLattice> &lattice() const { return workflow_.lattice(); }
  Eigen::VectorXd &wgap() const { return workflow_.wakeGap(); }
  BlSystemCoeffs &blc() const { return workflow_.systemCoefficients(); }
  BoundaryLayerTransitionSolver &transitionSolver() const {
    return workflow_.transitionSolver;
  }
  FlowRegimeEnum &flowRegime() const { return workflow_.flowRegime(); }
  BlCompressibilityParams &blCompressibility() const {
    return workflow_.compressibility();
  }
  BlReynoldsParams &blReynolds() const { return workflow_.reynolds(); }
  blDiff &xt() const { return workflow_.transitionSensitivity(); }

  int readSideStationCount(int side) const {
    return BoundaryLayerRuntimeStateOps::readSideStationCount(workflow_.lattice(),
                                                              side);
  }
  BoundaryLayerStationReadModel readStationModel(
      int side, int stationIndex) const {
    return BoundaryLayerRuntimeStateOps::readStationModel(
        workflow_.lattice(), workflow_.wakeGap(), side, stationIndex);
  }
  BoundaryLayerTrailingEdgeReadModel readTrailingEdgeModel() const {
    return BoundaryLayerRuntimeStateOps::readTrailingEdgeModel(
        workflow_.lattice());
  }

  void emitMarchInfoLog(std::string_view message) const {
    Logger::instance().write(std::string(message));
  }

  double readNewtonRhs(int row) const {
    return MrchueLinearSystemOps::readNewtonRhs(workflow_.systemCoefficients(),
                                                row);
  }
  void solveMrchueDirectNewtonSystem() const {
    MrchueLinearSystemOps::solveDirect(workflow_.systemCoefficients());
  }
  void solveMrchueInverseNewtonSystem(double htarg) const {
    MrchueLinearSystemOps::solveInverse(workflow_.systemCoefficients(),
                                        workflow_.state(), htarg);
  }
  void runTransitionCheckForMrchue(int side, int stationIndex, double &ami,
                                   double &cti) const {
    workflow_.runTransitionCheckForMrchue(side, stationIndex, ami, cti);
  }
  double calcHtarg(int stationIndex, int side, bool wake) const {
    return workflow_.calcHtarg(stationIndex, side, wake);
  }

  bool isStartOfWake(int side, int stationIndex) const {
    return stationIndex == workflow_.lattice(side).trailingEdgeIndex + 1;
  }
  FlowRegimeEnum applyFlowRegimeCandidate(FlowRegimeEnum candidate) const {
    workflow_.flowRegime() = candidate;
    return workflow_.flowRegime();
  }
  FlowRegimeEnum currentFlowRegime() const { return workflow_.flowRegime(); }
  FlowRegimeEnum determineRegimeForStation(int side, int stationIndex) const {
    return makeMixedModeOps().determineRegimeForStation(side, stationIndex);
  }
  int resetSideState(int side, const Foil &foil,
                     const StagnationResult &stagnation) const {
    return workflow_.resetSideState(side, foil, stagnation);
  }

  bool tesys(const BoundaryLayerSideProfiles &top_profiles,
             const BoundaryLayerSideProfiles &bottom_profiles,
             const Edge &edge) const {
    return solverOps_.tesys(top_profiles, bottom_profiles, edge);
  }
  bool blsys() const { return solverOps_.blsys(); }
  bool blkin(BoundaryLayerState &state) const {
    return solverOps_.blkin(state);
  }
  blData blprv(blData data, double xsi, double ami, double cti, double thi,
               double dsi, double dswaki, double uei) const {
    return solverOps_.blprv(data, xsi, ami, cti, thi, dsi, dswaki, uei);
  }

  void updateSystemMatricesForStation(
      const Edge &edge, int side, int stationIndex,
      BoundaryLayerMixedModeStationContext &ctx) const {
    makeMixedModeOps().updateSystemMatricesForStation(edge, side, stationIndex,
                                                      ctx);
  }
  void initializeFirstIterationState(
      int side, int stationIndex, int previousTransition,
      BoundaryLayerMixedModeStationContext &ctx, double &ueref,
      double &hkref) const {
    makeMixedModeOps().initializeFirstIterationState(side, stationIndex,
                                                     previousTransition, ctx,
                                                     ueref, hkref);
  }
  void configureSimilarityRow(double ueref) const {
    makeMixedModeOps().configureSimilarityRow(ueref);
  }
  void configureViscousRow(double hkref, double ueref, double senswt,
                           bool resetSensitivity, bool averageSensitivity,
                           double &sens, double &sennew) const {
    makeMixedModeOps().configureViscousRow(hkref, ueref, senswt,
                                           resetSensitivity,
                                           averageSensitivity, sens, sennew);
  }
  bool applyMixedModeNewtonStep(
      int side, int stationIndex, double &ami,
      BoundaryLayerMixedModeStationContext &ctx) const {
    return makeMixedModeOps().applyMixedModeNewtonStep(side, stationIndex, ami,
                                                       ctx);
  }
  void checkTransitionIfNeeded(int side, int stationIndex, bool skipCheck,
                               int laminarAdvance, double &ami) const {
    makeMixedModeOps().checkTransitionIfNeeded(side, stationIndex, skipCheck,
                                               laminarAdvance, ami);
  }
  void storeStationStateCommon(
      int side, int stationIndex,
      const BoundaryLayerMixedModeStationContext &ctx) const {
    makeMixedModeOps().storeStationStateCommon(side, stationIndex, ctx);
  }
  void recoverStationAfterFailure(
      int side, int stationIndex,
      BoundaryLayerMixedModeStationContext &ctx, double &ami,
      BoundaryLayerEdgeVelocityFallbackMode edgeMode,
      int laminarAdvance) const {
    makeMixedModeOps().recoverStationAfterFailure(side, stationIndex, ctx, ami,
                                                  edgeMode, laminarAdvance);
  }

 private:
  BoundaryLayerMixedModeOps makeMixedModeOps() const {
    return BoundaryLayerMixedModeOps({workflow_.lattice(),
                                      workflow_.state(),
                                      workflow_.flowRegime(),
                                      workflow_.systemCoefficients(),
                                      workflow_.compressibility(),
                                      workflow_.boundaryLayerVariablesSolver,
                                      workflow_.transitionSolver,
                                      solverOps_});
  }

  BoundaryLayerWorkflow &workflow_;
  BoundaryLayerSolverOps solverOps_;
};

inline BoundaryLayerMarchAccess makeBoundaryLayerMarchAccess(
    BoundaryLayerWorkflow &workflow) {
  return BoundaryLayerMarchAccess(workflow,
                                  makeBoundaryLayerSetblSolverOps(workflow));
}
