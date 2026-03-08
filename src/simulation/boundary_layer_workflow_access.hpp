#pragma once

#include "simulation/BoundaryLayer.hpp"

class BoundaryLayerWorkflowAccess {
 public:
  explicit BoundaryLayerWorkflowAccess(BoundaryLayerWorkflow &workflow)
      : workflow_(workflow) {}

  BoundaryLayerState &state() const { return workflow_.state; }
  SidePair<BoundaryLayerLattice> &lattice() const { return workflow_.lattice; }
  Eigen::VectorXd &wgap() const { return workflow_.wgap; }
  BlSystemCoeffs &blc() const { return workflow_.blc; }
  BoundaryLayerVariablesSolver &variableSolver() const {
    return workflow_.boundaryLayerVariablesSolver;
  }
  BoundaryLayerTransitionSolver &transitionSolver() const {
    return workflow_.transitionSolver;
  }
  FlowRegimeEnum &flowRegime() const { return workflow_.flowRegime; }
  BlCompressibilityParams &blCompressibility() const {
    return workflow_.blCompressibility;
  }
  BlReynoldsParams &blReynolds() const { return workflow_.blReynolds; }
  BlTransitionParams &blTransition() const { return workflow_.blTransition; }
  int systemSize() const { return workflow_.nsys; }

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
  double xifset(const Foil &foil, const StagnationResult &stagnation,
                int side) const {
    return workflow_.xifset(foil, stagnation, side);
  }
  SidePair<Eigen::VectorXd> ueset(const Eigen::MatrixXd &dij) const {
    return workflow_.ueset(dij);
  }
  FlowRegimeEnum determineRegimeForStation(int side, int station) const {
    return workflow_.determineRegimeForStation(side, station);
  }

 private:
  BoundaryLayerWorkflow &workflow_;
};
