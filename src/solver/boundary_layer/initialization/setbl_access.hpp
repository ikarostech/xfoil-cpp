#pragma once

#include "solver/boundary_layer/boundary_layer.hpp"
#include "solver/boundary_layer/boundary_layer_aerodynamics.hpp"
#include "solver/boundary_layer/workflow/mixed_mode.hpp"
#include "solver/boundary_layer/runtime/state.hpp"
#include "solver/boundary_layer/workflow/solver_ops.hpp"

class BoundaryLayerSetblAccess {
 public:
  explicit BoundaryLayerSetblAccess(BoundaryLayer &boundaryLayer)
      : boundaryLayer_(boundaryLayer) {}

  Eigen::VectorXd &wgap() const { return boundaryLayer_.wakeGap(); }
  BlSystemCoeffs &blc() const { return boundaryLayer_.systemCoefficients(); }
  const BlSystemCoeffs &blcConst() const {
    return boundaryLayer_.systemCoefficients();
  }
  int stationCount(int side) const { return boundaryLayer_.stationCount(side); }
  int trailingEdgeIndex(int side) const { return boundaryLayer_.trailingEdgeIndex(side); }
  int topTrailingEdgeIndex() const { return trailingEdgeIndex(1); }
  int bottomTrailingEdgeIndex() const { return trailingEdgeIndex(2); }
  int transitionIndex(int side) const { return boundaryLayer_.transitionIndex(side); }
  void setTransitionIndex(int side, int stationIndex) const {
    boundaryLayer_.setTransitionIndex(side, stationIndex);
  }
  int stationToPanel(int side, int stationIndex) const {
    return boundaryLayer_.stationToPanel(side, stationIndex);
  }
  int leadingEdgePanelIndex(int side) const { return stationToPanel(side, 0); }
  int trailingEdgePanelIndex(int side) const {
    return stationToPanel(side, trailingEdgeIndex(side));
  }
  int stationToSystem(int side, int stationIndex) const {
    return boundaryLayer_.stationToSystem(side, stationIndex);
  }
  int trailingEdgeSystemIndex(int side) const {
    return boundaryLayer_.trailingEdgeSystemIndex(side);
  }
  double arcLengthCoordinate(int side, int stationIndex) const {
    return boundaryLayer_.readStationModel(side, stationIndex).arcLength;
  }
  double panelInfluenceFactor(int side, int stationIndex) const {
    return boundaryLayer_.panelInfluenceFactor(side, stationIndex);
  }
  double leadingEdgePanelInfluenceFactor(int side) const {
    return panelInfluenceFactor(side, 0);
  }
  double trailingEdgePanelInfluenceFactor(int side) const {
    return panelInfluenceFactor(side, trailingEdgeIndex(side));
  }
  double inviscidEdgeVelocitySensitivityToAlpha(int side,
                                                int stationIndex) const {
    return boundaryLayer_.inviscidEdgeVelocitySensitivityToAlpha(side, stationIndex);
  }
  double wakeGapAt(int side, int stationIndex) const {
    const int wake_index = stationIndex - trailingEdgeIndex(side);
    return boundaryLayer_.wakeGap()[wake_index - 1];
  }
  BoundaryLayerSideProfiles &profiles(int side) { return boundaryLayer_.profiles(side); }
  const BoundaryLayerSideProfiles &profiles(int side) const { return boundaryLayer_.profiles(side); }
  void copyProfilesTo(SidePair<BoundaryLayerSideProfiles> &profiles) const {
    boundaryLayer_.copyProfilesTo(profiles);
  }
  FlowRegimeEnum &flowRegime() const { return boundaryLayer_.flowRegime(); }
  BlCompressibilityParams &blCompressibility() const {
    return boundaryLayer_.compressibility();
  }
  BlReynoldsParams &blReynolds() const { return boundaryLayer_.reynolds(); }
  BlTransitionParams &blTransition() const { return boundaryLayer_.transition(); }
  int systemSize() const { return boundaryLayer_.systemSize(); }
  FlowRegimeEnum determineRegimeForStation(int side, int stationIndex) const {
    return boundaryLayer_.determineRegimeForStation(side, stationIndex);
  }

  bool tesys(const BoundaryLayerSideProfiles &top_profiles,
             const BoundaryLayerSideProfiles &bottom_profiles,
             const Edge &edge) const {
    return boundaryLayer_.tesys(top_profiles, bottom_profiles, edge);
  }
  bool blsys() const { return boundaryLayer_.blsys(); }
  bool blkin(BoundaryLayerState &state) const { return boundaryLayer_.blkin(state); }
  blData blprv(blData data, double xsi, double ami, double cti, double thi,
               double dsi, double dswaki, double uei) const {
    return boundaryLayer_.blprv(data, xsi, ami, cti, thi, dsi, dswaki, uei);
  }
  SkinFrictionCoefficients blmid(FlowRegimeEnum flow_regime) const {
    return boundaryLayer_.blmid(flow_regime);
  }
  void runTransitionCheck() const { boundaryLayer_.runTransitionCheck(); }
  double currentAmplification() const { return boundaryLayer_.currentAmplification(); }
  double previousAmplification() const { return boundaryLayer_.previousAmplification(); }
  double currentSkinFrictionHistory() const { return boundaryLayer_.currentSkinFrictionHistory(); }
  BoundaryLayerState snapshotState() const { return boundaryLayer_.snapshotState(); }
  void replaceState(const BoundaryLayerState &state) const { boundaryLayer_.replaceState(state); }
  void advanceState() const { boundaryLayer_.advanceState(); }
  bool solveTeSystemForCurrentProfiles(const Edge &edge) const {
    return boundaryLayer_.solveTeSystemForCurrentProfiles(edge);
  }
  void solveWakeState() const { boundaryLayer_.solveWakeState(); }
  double xifset(const Foil &foil, const StagnationResult &stagnation,
                int side) const {
    return boundaryLayer_.computeForcedTransitionArcLength(foil, stagnation, side);
  }
  SidePair<Eigen::VectorXd> ueset(const Eigen::MatrixXd &dij) const {
    return boundaryLayer_.computeInviscidEdgeVelocitySensitivity(dij);
  }

 private:
  BoundaryLayer &boundaryLayer_;
};
