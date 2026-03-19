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
  const BlSystemCoeffs &blcConst() const {
    return workflow_.systemCoefficients();
  }
  int stationCount(int side) const { return workflow_.lattice(side).stationCount; }
  int trailingEdgeIndex(int side) const {
    return workflow_.lattice(side).trailingEdgeIndex;
  }
  int topTrailingEdgeIndex() const { return trailingEdgeIndex(1); }
  int bottomTrailingEdgeIndex() const { return trailingEdgeIndex(2); }
  int transitionIndex(int side) const {
    return workflow_.lattice(side).profiles.transitionIndex;
  }
  void setTransitionIndex(int side, int stationIndex) const {
    workflow_.lattice(side).profiles.transitionIndex = stationIndex;
  }
  int stationToPanel(int side, int stationIndex) const {
    return workflow_.lattice(side).stationToPanel[stationIndex];
  }
  int leadingEdgePanelIndex(int side) const { return stationToPanel(side, 0); }
  int trailingEdgePanelIndex(int side) const {
    return stationToPanel(side, trailingEdgeIndex(side));
  }
  int stationToSystem(int side, int stationIndex) const {
    return workflow_.lattice(side).stationToSystem[stationIndex];
  }
  int trailingEdgeSystemIndex(int side) const {
    return stationToSystem(side, trailingEdgeIndex(side));
  }
  double arcLengthCoordinate(int side, int stationIndex) const {
    return workflow_.lattice(side).arcLengthCoordinates[stationIndex];
  }
  double panelInfluenceFactor(int side, int stationIndex) const {
    return workflow_.lattice(side).panelInfluenceFactor[stationIndex];
  }
  double leadingEdgePanelInfluenceFactor(int side) const {
    return panelInfluenceFactor(side, 0);
  }
  double trailingEdgePanelInfluenceFactor(int side) const {
    return panelInfluenceFactor(side, trailingEdgeIndex(side));
  }
  double inviscidEdgeVelocitySensitivityToAlpha(int side,
                                                int stationIndex) const {
    return workflow_.lattice(side).inviscidEdgeVelocityMatrix(1, stationIndex);
  }
  double wakeGapAt(int side, int stationIndex) const {
    const int wake_index = stationIndex - trailingEdgeIndex(side);
    return workflow_.wakeGap()[wake_index - 1];
  }
  BoundaryLayerSideProfiles &profiles(int side) {
    return workflow_.lattice(side).profiles;
  }
  const BoundaryLayerSideProfiles &profiles(int side) const {
    return workflow_.lattice(side).profiles;
  }
  void copyProfilesTo(SidePair<BoundaryLayerSideProfiles> &profiles) const {
    profiles.top = workflow_.lattice(1).profiles;
    profiles.bottom = workflow_.lattice(2).profiles;
  }
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
  double currentAmplification() const {
    return workflow_.state().station2.param.amplz;
  }
  double currentSkinFrictionHistory() const {
    return workflow_.state().station2.cqz.scalar;
  }
  void replaceState(const BoundaryLayerState &state) const {
    workflow_.state() = state;
  }
  void advanceState() const { workflow_.state().stepbl(); }
  bool solveTeSystemForCurrentProfiles(const Edge &edge) const {
    return solverOps_.tesys(workflow_.lattice(1).profiles,
                            workflow_.lattice(2).profiles, edge);
  }
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
