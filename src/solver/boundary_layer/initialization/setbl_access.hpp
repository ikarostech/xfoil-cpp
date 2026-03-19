#pragma once

#include "solver/boundary_layer/workflow/workflow.hpp"
#include "solver/boundary_layer/boundary_layer_aerodynamics.hpp"
#include "solver/boundary_layer/workflow/mixed_mode.hpp"
#include "solver/boundary_layer/runtime/state.hpp"
#include "solver/boundary_layer/workflow/solver_ops.hpp"

class BoundaryLayerSetblAccess {
 public:
  BoundaryLayerSetblAccess(BoundaryLayerWorkflow &workflow,
                           BoundaryLayerSolverOps solver_ops)
      : workflow_(workflow), solverOps_(solver_ops) {}

  BoundaryLayerState &state() const { return workflow_.state(); }
  SidePair<BoundaryLayerLattice> &lattice() const { return workflow_.lattice(); }
  Eigen::VectorXd &wgap() const { return workflow_.wakeGap(); }
  BlSystemCoeffs &blc() const { return workflow_.systemCoefficients(); }
  FlowRegimeEnum &flowRegime() const { return workflow_.flowRegime(); }
  BlCompressibilityParams &blCompressibility() const {
    return workflow_.compressibility();
  }
  BlReynoldsParams &blReynolds() const { return workflow_.reynolds(); }
  BlTransitionParams &blTransition() const { return workflow_.transition(); }
  int systemSize() const { return workflow_.systemSize(); }
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
  void runTransitionCheck() const { workflow_.transitionSolver.trchek(); }
  void solveWakeState() const {
    workflow_.state().station2 = workflow_.boundaryLayerVariablesSolver.solve(
        workflow_.state().station2, FlowRegimeEnum::Wake);
  }
  double xifset(const Foil &foil, const StagnationResult &stagnation,
                int side) const {
    return workflow_.computeForcedTransitionArcLength(foil, stagnation, side);
  }
  SidePair<Eigen::VectorXd> ueset(const Eigen::MatrixXd &dij) const {
    return BoundaryLayerAerodynamicsOps::ueset(workflow_.lattice(), dij);
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

inline BoundaryLayerSolverOps makeBoundaryLayerSetblSolverOps(
    BoundaryLayerWorkflow &workflow) {
  return BoundaryLayerSolverOps({workflow.boundaryLayerVariablesSolver,
                                 workflow.blDiffSolver,
                                 workflow.transitionSolver,
                                 workflow.flowRegime(),
                                 workflow.compressibility(),
                                 workflow.reynolds(),
                                 workflow.transition(),
                                 workflow.state(),
                                 workflow.systemCoefficients(),
                                 workflow.lattice()});
}

inline BoundaryLayerSetblAccess makeBoundaryLayerSetblAccess(
    BoundaryLayerWorkflow &workflow) {
  return BoundaryLayerSetblAccess(workflow,
                                  makeBoundaryLayerSetblSolverOps(workflow));
}
