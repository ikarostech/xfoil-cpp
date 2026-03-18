#pragma once

#include <string_view>

#include "infrastructure/logger.hpp"
#include "solver/boundary_layer/runtime/state.hpp"
#include "solver/march/mrchue_linear_system.hpp"
#include "solver/boundary_layer/initialization/setbl_access.hpp"

struct BoundaryLayerMarchContextData {
  BoundaryLayerState &state;
  SidePair<BoundaryLayerLattice> &lattice;
  Eigen::VectorXd &wgap;
  BlSystemCoeffs &blc;
  blDiff &xt;
};

class BoundaryLayerMarchAccess {
 public:
  BoundaryLayerMarchAccess(BoundaryLayerMarchContextData context,
                           BoundaryLayerVariablesSolver &variable_solver,
                           BoundaryLayerTransitionSolver &transition_solver,
                           FlowRegimeEnum &flow_regime,
                           BlCompressibilityParams &bl_compressibility,
                           BlReynoldsParams &bl_reynolds,
                           BlTransitionParams &bl_transition,
                           BoundaryLayerSolverOps solver_ops)
      : context_(context),
        variableSolver_(variable_solver),
        transitionSolver_(transition_solver),
        flowRegime_(flow_regime),
        blCompressibility_(bl_compressibility),
        blReynolds_(bl_reynolds),
        blTransition_(bl_transition),
        solverOps_(solver_ops) {}

  BoundaryLayerState &state() const { return context_.state; }
  SidePair<BoundaryLayerLattice> &lattice() const { return context_.lattice; }
  Eigen::VectorXd &wgap() const { return context_.wgap; }
  BlSystemCoeffs &blc() const { return context_.blc; }
  BoundaryLayerTransitionSolver &transitionSolver() const {
    return transitionSolver_;
  }
  FlowRegimeEnum &flowRegime() const { return flowRegime_; }
  BlCompressibilityParams &blCompressibility() const {
    return blCompressibility_;
  }
  BlReynoldsParams &blReynolds() const { return blReynolds_; }
  blDiff &xt() const { return context_.xt; }

  int readSideStationCount(int side) const {
    return BoundaryLayerRuntimeStateOps::readSideStationCount(context_.lattice,
                                                              side);
  }
  BoundaryLayerStationReadModel readStationModel(
      int side, int stationIndex) const {
    return BoundaryLayerRuntimeStateOps::readStationModel(
        context_.lattice, context_.wgap, side, stationIndex);
  }
  BoundaryLayerTrailingEdgeReadModel readTrailingEdgeModel() const {
    return BoundaryLayerRuntimeStateOps::readTrailingEdgeModel(context_.lattice);
  }

  void emitMarchInfoLog(std::string_view message) const {
    Logger::instance().write(std::string(message));
  }

  double readNewtonRhs(int row) const {
    return MrchueLinearSystemOps::readNewtonRhs(context_.blc, row);
  }
  void solveMrchueDirectNewtonSystem() const {
    MrchueLinearSystemOps::solveDirect(context_.blc);
  }
  void solveMrchueInverseNewtonSystem(double htarg) const {
    MrchueLinearSystemOps::solveInverse(context_.blc, context_.state, htarg);
  }

  bool isStartOfWake(int side, int stationIndex) const {
    return stationIndex == context_.lattice.get(side).trailingEdgeIndex + 1;
  }
  FlowRegimeEnum applyFlowRegimeCandidate(FlowRegimeEnum candidate) const {
    flowRegime_ = candidate;
    return flowRegime_;
  }
  FlowRegimeEnum currentFlowRegime() const { return flowRegime_; }
  FlowRegimeEnum determineRegimeForStation(int side, int stationIndex) const {
    return makeMixedModeOps().determineRegimeForStation(side, stationIndex);
  }
  int resetSideState(int side, const Foil &foil,
                     const StagnationResult &stagnation) const {
    return BoundaryLayerRuntimeStateOps::resetSideState(
        context_.lattice, blTransition_, flowRegime_, side, foil, stagnation);
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
    return BoundaryLayerMixedModeOps({context_.lattice,
                                      context_.state,
                                      flowRegime_,
                                      context_.blc,
                                      blCompressibility_,
                                      variableSolver_,
                                      transitionSolver_,
                                      solverOps_});
  }

  BoundaryLayerMarchContextData context_;
  BoundaryLayerVariablesSolver &variableSolver_;
  BoundaryLayerTransitionSolver &transitionSolver_;
  FlowRegimeEnum &flowRegime_;
  BlCompressibilityParams &blCompressibility_;
  BlReynoldsParams &blReynolds_;
  BlTransitionParams &blTransition_;
  BoundaryLayerSolverOps solverOps_;
};

inline BoundaryLayerMarchContextData makeBoundaryLayerMarchContext(
    BoundaryLayerWorkflow &workflow) {
  auto &state_store = workflow.stateStore();
  auto &workspace = workflow.workspace();
  return {workspace.state,
          state_store.lattice,
          state_store.wgap,
          workspace.blc,
          workspace.xt};
}

inline BoundaryLayerMarchAccess makeBoundaryLayerMarchAccess(
    BoundaryLayerWorkflow &workflow) {
  return BoundaryLayerMarchAccess(makeBoundaryLayerMarchContext(workflow),
                                  workflow.boundaryLayerVariablesSolver,
                                  workflow.transitionSolver,
                                  workflow.stateStore().flowRegime,
                                  workflow.stateStore().blCompressibility,
                                  workflow.stateStore().blReynolds,
                                  workflow.stateStore().blTransition,
                                  makeBoundaryLayerSetblSolverOps(workflow));
}
