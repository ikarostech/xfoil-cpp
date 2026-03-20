#pragma once

#include "solver/boundary_layer/workflow/workflow.hpp"
#include "solver/boundary_layer/boundary_layer_aerodynamics.hpp"
#include "solver/boundary_layer/workflow/mixed_mode.hpp"
#include "solver/boundary_layer/runtime/state.hpp"
#include "solver/boundary_layer/workflow/solver_ops.hpp"

class BoundaryLayerSetblAccess {
 public:
  explicit BoundaryLayerSetblAccess(BoundaryLayerWorkflow &workflow)
      : workflow_(workflow) {}

  Eigen::VectorXd &wgap() const { return workflow_.wakeGap(); }
  BlSystemCoeffs &blc() const { return workflow_.systemCoefficients(); }
  const BlSystemCoeffs &blcConst() const {
    return workflow_.systemCoefficients();
  }
  int stationCount(int side) const { return workflow_.stationCount(side); }
  int trailingEdgeIndex(int side) const { return workflow_.trailingEdgeIndex(side); }
  int topTrailingEdgeIndex() const { return trailingEdgeIndex(1); }
  int bottomTrailingEdgeIndex() const { return trailingEdgeIndex(2); }
  int transitionIndex(int side) const { return workflow_.transitionIndex(side); }
  void setTransitionIndex(int side, int stationIndex) const {
    workflow_.setTransitionIndex(side, stationIndex);
  }
  int stationToPanel(int side, int stationIndex) const {
    return workflow_.stationToPanel(side, stationIndex);
  }
  int leadingEdgePanelIndex(int side) const { return stationToPanel(side, 0); }
  int trailingEdgePanelIndex(int side) const {
    return stationToPanel(side, trailingEdgeIndex(side));
  }
  int stationToSystem(int side, int stationIndex) const {
    return workflow_.stationToSystem(side, stationIndex);
  }
  int trailingEdgeSystemIndex(int side) const {
    return workflow_.trailingEdgeSystemIndex(side);
  }
  double arcLengthCoordinate(int side, int stationIndex) const {
    return workflow_.readStationModel(side, stationIndex).arcLength;
  }
  double panelInfluenceFactor(int side, int stationIndex) const {
    return workflow_.panelInfluenceFactor(side, stationIndex);
  }
  double leadingEdgePanelInfluenceFactor(int side) const {
    return panelInfluenceFactor(side, 0);
  }
  double trailingEdgePanelInfluenceFactor(int side) const {
    return panelInfluenceFactor(side, trailingEdgeIndex(side));
  }
  double inviscidEdgeVelocitySensitivityToAlpha(int side,
                                                int stationIndex) const {
    return workflow_.inviscidEdgeVelocitySensitivityToAlpha(side, stationIndex);
  }
  double wakeGapAt(int side, int stationIndex) const {
    const int wake_index = stationIndex - trailingEdgeIndex(side);
    return workflow_.wakeGap()[wake_index - 1];
  }
  BoundaryLayerSideProfiles &profiles(int side) { return workflow_.profiles(side); }
  const BoundaryLayerSideProfiles &profiles(int side) const { return workflow_.profiles(side); }
  void copyProfilesTo(SidePair<BoundaryLayerSideProfiles> &profiles) const {
    workflow_.copyProfilesTo(profiles);
  }
  FlowRegimeEnum &flowRegime() const { return workflow_.flowRegime(); }
  BlCompressibilityParams &blCompressibility() const {
    return workflow_.compressibility();
  }
  BlReynoldsParams &blReynolds() const { return workflow_.reynolds(); }
  BlTransitionParams &blTransition() const { return workflow_.transition(); }
  int systemSize() const { return workflow_.systemSize(); }
  FlowRegimeEnum determineRegimeForStation(int side, int stationIndex) const {
    return workflow_.determineRegimeForStation(side, stationIndex);
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
  SkinFrictionCoefficients blmid(FlowRegimeEnum flow_regime) const {
    return workflow_.blmid(flow_regime);
  }
  void runTransitionCheck() const { workflow_.runTransitionCheck(); }
  double currentAmplification() const { return workflow_.currentAmplification(); }
  double previousAmplification() const { return workflow_.previousAmplification(); }
  double currentSkinFrictionHistory() const { return workflow_.currentSkinFrictionHistory(); }
  BoundaryLayerState snapshotState() const { return workflow_.snapshotState(); }
  void replaceState(const BoundaryLayerState &state) const { workflow_.replaceState(state); }
  void advanceState() const { workflow_.advanceState(); }
  bool solveTeSystemForCurrentProfiles(const Edge &edge) const {
    return workflow_.solveTeSystemForCurrentProfiles(edge);
  }
  void solveWakeState() const { workflow_.solveWakeState(); }
  double xifset(const Foil &foil, const StagnationResult &stagnation,
                int side) const {
    return workflow_.computeForcedTransitionArcLength(foil, stagnation, side);
  }
  SidePair<Eigen::VectorXd> ueset(const Eigen::MatrixXd &dij) const {
    return workflow_.computeInviscidEdgeVelocitySensitivity(dij);
  }

 private:
  BoundaryLayerWorkflow &workflow_;
};

inline BoundaryLayerSetblAccess makeBoundaryLayerSetblAccess(
    BoundaryLayerWorkflow &workflow) {
  return BoundaryLayerSetblAccess(workflow);
}
