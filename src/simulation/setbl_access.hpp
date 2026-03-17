#pragma once

#include "simulation/BoundaryLayer.hpp"
#include "simulation/boundary_layer_aerodynamics.hpp"
#include "simulation/boundary_layer_mixed_mode.hpp"
#include "simulation/boundary_layer_runtime_state.hpp"
#include "simulation/boundary_layer_solver_ops.hpp"

struct BoundaryLayerSetblContext {
  BoundaryLayerState &state;
  SidePair<BoundaryLayerLattice> &lattice;
  Eigen::VectorXd &wgap;
  BlSystemCoeffs &blc;
  int &nsys;
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
