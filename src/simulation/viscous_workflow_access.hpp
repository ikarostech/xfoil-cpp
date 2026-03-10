#pragma once

#include <string_view>

#include "infrastructure/logger.hpp"
#include "simulation/BoundaryLayer.hpp"
#include "simulation/boundary_layer_aerodynamics.hpp"
#include "simulation/boundary_layer_mixed_mode.hpp"
#include "simulation/boundary_layer_runtime_state.hpp"
#include "simulation/boundary_layer_solver_ops.hpp"
#include "simulation/march/mrchue_linear_system.hpp"

struct BoundaryLayerSetblContext {
  BoundaryLayerState &state;
  SidePair<BoundaryLayerLattice> &lattice;
  Eigen::VectorXd &wgap;
  BlSystemCoeffs &blc;
  int &nsys;
};

struct BoundaryLayerMarchContextData {
  BoundaryLayerState &state;
  SidePair<BoundaryLayerLattice> &lattice;
  Eigen::VectorXd &wgap;
  BlSystemCoeffs &blc;
  BoundaryLayerVariablesSolver &variableSolver;
  BlDiffSolver &blDiffSolver;
  BoundaryLayerTransitionSolver &transitionSolver;
  FlowRegimeEnum &flowRegime;
  BlCompressibilityParams &blCompressibility;
  BlReynoldsParams &blReynolds;
  BlTransitionParams &blTransition;
  int &nsys;
  blDiff &xt;
};

class BoundaryLayerSetblAccess {
 public:
  BoundaryLayerSetblAccess(BoundaryLayerSetblContext context,
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

  BoundaryLayerSetblContext context() const { return context_; }

  BoundaryLayerState &state() const { return context_.state; }
  SidePair<BoundaryLayerLattice> &lattice() const { return context_.lattice; }
  Eigen::VectorXd &wgap() const { return context_.wgap; }
  BlSystemCoeffs &blc() const { return context_.blc; }
  FlowRegimeEnum &flowRegime() const { return flowRegime_; }
  BlCompressibilityParams &blCompressibility() const {
    return blCompressibility_;
  }
  BlReynoldsParams &blReynolds() const { return blReynolds_; }
  BlTransitionParams &blTransition() const { return blTransition_; }
  int systemSize() const { return context_.nsys; }
  FlowRegimeEnum determineRegimeForStation(int side, int stationIndex) const {
    return makeMixedModeOps().determineRegimeForStation(side, stationIndex);
  }

  bool tesys(const BoundaryLayerSideProfiles &top_profiles,
             const BoundaryLayerSideProfiles &bottom_profiles,
             const Edge &edge) const {
    return solverOps_.tesys(top_profiles, bottom_profiles, edge);
  }
  bool blsys() const { return solverOps_.blsys(); }
  bool blkin(BoundaryLayerState &state) const { return solverOps_.blkin(state); }
  blData blprv(blData data, double xsi, double ami, double cti, double thi,
               double dsi, double dswaki, double uei) const {
    return solverOps_.blprv(data, xsi, ami, cti, thi, dsi, dswaki, uei);
  }
  SkinFrictionCoefficients blmid(FlowRegimeEnum flow_regime) const {
    return solverOps_.blmid(flow_regime);
  }
  void runTransitionCheck() const { transitionSolver_.trchek(); }
  void solveWakeState() const {
    context_.state.station2 =
        variableSolver_.solve(context_.state.station2, FlowRegimeEnum::Wake);
  }
  double xifset(const Foil &foil, const StagnationResult &stagnation,
                int side) const {
    return BoundaryLayerRuntimeStateOps::xifset(context_.lattice, foil,
                                                stagnation, side);
  }
  SidePair<Eigen::VectorXd> ueset(const Eigen::MatrixXd &dij) const {
    return BoundaryLayerAerodynamicsOps::ueset(context_.lattice, dij);
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

  BoundaryLayerSetblContext context_;
  BoundaryLayerVariablesSolver &variableSolver_;
  BoundaryLayerTransitionSolver &transitionSolver_;
  FlowRegimeEnum &flowRegime_;
  BlCompressibilityParams &blCompressibility_;
  BlReynoldsParams &blReynolds_;
  BlTransitionParams &blTransition_;
  BoundaryLayerSolverOps solverOps_;
};

class BoundaryLayerMarchAccess {
 public:
  explicit BoundaryLayerMarchAccess(BoundaryLayerMarchContextData context)
      : context_(context) {}

  BoundaryLayerState &state() const { return context_.state; }
  SidePair<BoundaryLayerLattice> &lattice() const { return context_.lattice; }
  Eigen::VectorXd &wgap() const { return context_.wgap; }
  BlSystemCoeffs &blc() const { return context_.blc; }
  BoundaryLayerTransitionSolver &transitionSolver() const {
    return context_.transitionSolver;
  }
  FlowRegimeEnum &flowRegime() const { return context_.flowRegime; }
  BlCompressibilityParams &blCompressibility() const {
    return context_.blCompressibility;
  }
  BlReynoldsParams &blReynolds() const { return context_.blReynolds; }
  int systemSize() const { return context_.nsys; }
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
    context_.flowRegime = candidate;
    return context_.flowRegime;
  }
  FlowRegimeEnum currentFlowRegime() const { return context_.flowRegime; }
  FlowRegimeEnum determineRegimeForStation(int side, int stationIndex) const {
    return makeMixedModeOps().determineRegimeForStation(side, stationIndex);
  }
  int resetSideState(int side, const Foil &foil,
                     const StagnationResult &stagnation) const {
    return BoundaryLayerRuntimeStateOps::resetSideState(
        context_.lattice, context_.blTransition, context_.flowRegime, side,
        foil, stagnation);
  }

  bool tesys(const BoundaryLayerSideProfiles &top_profiles,
             const BoundaryLayerSideProfiles &bottom_profiles,
             const Edge &edge) const {
    return makeSolverOps().tesys(top_profiles, bottom_profiles, edge);
  }
  bool blsys() const { return makeSolverOps().blsys(); }
  bool blkin(BoundaryLayerState &state) const {
    return makeSolverOps().blkin(state);
  }
  blData blprv(blData data, double xsi, double ami, double cti, double thi,
               double dsi, double dswaki, double uei) const {
    return makeSolverOps().blprv(data, xsi, ami, cti, thi, dsi, dswaki, uei);
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
  BoundaryLayerSolverOps makeSolverOps() const {
    return BoundaryLayerSolverOps({context_.variableSolver,
                                   context_.blDiffSolver,
                                   context_.transitionSolver,
                                   context_.flowRegime,
                                   context_.blCompressibility,
                                   context_.blReynolds,
                                   context_.blTransition,
                                   context_.state,
                                   context_.blc,
                                   context_.lattice});
  }

  BoundaryLayerMixedModeOps makeMixedModeOps() const {
    return BoundaryLayerMixedModeOps({context_.lattice,
                                      context_.state,
                                      context_.flowRegime,
                                      context_.blc,
                                      context_.blCompressibility,
                                      context_.variableSolver,
                                      context_.transitionSolver,
                                      makeSolverOps()});
  }

  BoundaryLayerMarchContextData context_;
};

inline BoundaryLayerSetblContext makeBoundaryLayerSetblContext(
    BoundaryLayerWorkflow &workflow) {
  auto &state_store = workflow.stateStore();
  auto &workspace = workflow.workspace();
  return {workspace.state,
          state_store.lattice,
          state_store.wgap,
          workspace.blc,
          workspace.nsys};
}

inline BoundaryLayerSolverOps makeBoundaryLayerSetblSolverOps(
    BoundaryLayerWorkflow &workflow) {
  auto &state_store = workflow.stateStore();
  auto &workspace = workflow.workspace();
  return BoundaryLayerSolverOps({workflow.boundaryLayerVariablesSolver,
                                 workflow.blDiffSolver,
                                 workflow.transitionSolver,
                                 state_store.flowRegime,
                                 state_store.blCompressibility,
                                 state_store.blReynolds,
                                 state_store.blTransition,
                                 workspace.state,
                                 workspace.blc,
                                 state_store.lattice});
}

inline BoundaryLayerMarchContextData makeBoundaryLayerMarchContext(
    BoundaryLayerWorkflow &workflow) {
  auto &state_store = workflow.stateStore();
  auto &workspace = workflow.workspace();
  return {workspace.state,
          state_store.lattice,
          state_store.wgap,
          workspace.blc,
          workflow.boundaryLayerVariablesSolver,
          workflow.blDiffSolver,
          workflow.transitionSolver,
          state_store.flowRegime,
          state_store.blCompressibility,
          state_store.blReynolds,
          state_store.blTransition,
          workspace.nsys,
          workspace.xt};
}

inline BoundaryLayerSetblAccess makeBoundaryLayerSetblAccess(
    BoundaryLayerWorkflow &workflow) {
  return BoundaryLayerSetblAccess(makeBoundaryLayerSetblContext(workflow),
                                  workflow.boundaryLayerVariablesSolver,
                                  workflow.transitionSolver,
                                  workflow.stateStore().flowRegime,
                                  workflow.stateStore().blCompressibility,
                                  workflow.stateStore().blReynolds,
                                  workflow.stateStore().blTransition,
                                  makeBoundaryLayerSetblSolverOps(workflow));
}

inline BoundaryLayerMarchAccess makeBoundaryLayerMarchAccess(
    BoundaryLayerWorkflow &workflow) {
  return BoundaryLayerMarchAccess(makeBoundaryLayerMarchContext(workflow));
}
