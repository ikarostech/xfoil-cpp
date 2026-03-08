#pragma once

#include <string_view>

#include "infrastructure/logger.hpp"
#include "simulation/BoundaryLayer.hpp"
#include "simulation/boundary_layer_aerodynamics.hpp"
#include "simulation/boundary_layer_mixed_mode.hpp"
#include "simulation/boundary_layer_runtime_state.hpp"
#include "simulation/boundary_layer_solver_ops.hpp"

class BoundaryLayerWorkflowAccess {
 public:
  struct Context {
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

  explicit BoundaryLayerWorkflowAccess(Context context) : context_(context) {}

  BoundaryLayerState &state() const { return context_.state; }
  SidePair<BoundaryLayerLattice> &lattice() const { return context_.lattice; }
  Eigen::VectorXd &wgap() const { return context_.wgap; }
  BlSystemCoeffs &blc() const { return context_.blc; }
  BoundaryLayerVariablesSolver &variableSolver() const {
    return context_.variableSolver;
  }
  BoundaryLayerTransitionSolver &transitionSolver() const {
    return context_.transitionSolver;
  }
  FlowRegimeEnum &flowRegime() const { return context_.flowRegime; }
  BlCompressibilityParams &blCompressibility() const {
    return context_.blCompressibility;
  }
  BlReynoldsParams &blReynolds() const { return context_.blReynolds; }
  BlTransitionParams &blTransition() const { return context_.blTransition; }
  int systemSize() const { return context_.nsys; }
  blDiff &xt() const { return context_.xt; }

  int readSideStationCount(int side) const {
    return BoundaryLayerRuntimeStateOps::readSideStationCount(context_.lattice,
                                                              side);
  }
  BoundaryLayerWorkflow::StationReadModel readStationModel(
      int side, int stationIndex) const {
    return BoundaryLayerRuntimeStateOps::readStationModel(
        context_.lattice, context_.wgap, side, stationIndex);
  }
  BoundaryLayerWorkflow::TrailingEdgeReadModel readTrailingEdgeModel() const {
    return BoundaryLayerRuntimeStateOps::readTrailingEdgeModel(context_.lattice);
  }

  void emitMarchInfoLog(std::string_view message) const {
    Logger::instance().write(std::string(message));
  }

  double readNewtonRhs(int row) const { return context_.blc.rhs[row]; }
  void solveMrchueDirectNewtonSystem() const {
    context_.blc.a2(3, 0) = 0.0;
    context_.blc.a2(3, 1) = 0.0;
    context_.blc.a2(3, 2) = 0.0;
    context_.blc.a2(3, 3) = 1.0;
    context_.blc.rhs[3] = 0.0;
    context_.blc.rhs =
        context_.blc.a2.block(0, 0, 4, 4).fullPivLu().solve(context_.blc.rhs);
  }
  void solveMrchueInverseNewtonSystem(double htarg) const {
    context_.blc.a2(3, 0) = 0.0;
    context_.blc.a2(3, 1) = context_.state.station2.hkz.t();
    context_.blc.a2(3, 2) = context_.state.station2.hkz.d();
    context_.blc.a2(3, 3) = context_.state.station2.hkz.u();
    context_.blc.rhs[3] = htarg - context_.state.station2.hkz.scalar;
    context_.blc.rhs =
        context_.blc.a2.block(0, 0, 4, 4).fullPivLu().solve(context_.blc.rhs);
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
  bool blkin(BoundaryLayerState &state) const { return makeSolverOps().blkin(state); }
  blData blprv(blData data, double xsi, double ami, double cti, double thi,
               double dsi, double dswaki, double uei) const {
    return makeSolverOps().blprv(data, xsi, ami, cti, thi, dsi, dswaki, uei);
  }
  SkinFrictionCoefficients blmid(FlowRegimeEnum flow_regime) const {
    return makeSolverOps().blmid(flow_regime);
  }
  double xifset(const Foil &foil, const StagnationResult &stagnation,
                int side) const {
    return BoundaryLayerRuntimeStateOps::xifset(context_.lattice, foil,
                                                stagnation, side);
  }
  SidePair<Eigen::VectorXd> ueset(const Eigen::MatrixXd &dij) const {
    return BoundaryLayerAerodynamicsOps::ueset(context_.lattice, dij);
  }

  void updateSystemMatricesForStation(
      const Edge &edge, int side, int stationIndex,
      BoundaryLayerWorkflow::MixedModeStationContext &ctx) const {
    makeMixedModeOps().updateSystemMatricesForStation(edge, side, stationIndex,
                                                      ctx);
  }
  void initializeFirstIterationState(
      int side, int stationIndex, int previousTransition,
      BoundaryLayerWorkflow::MixedModeStationContext &ctx, double &ueref,
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
      BoundaryLayerWorkflow::MixedModeStationContext &ctx) const {
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
      const BoundaryLayerWorkflow::MixedModeStationContext &ctx) const {
    makeMixedModeOps().storeStationStateCommon(side, stationIndex, ctx);
  }
  void recoverStationAfterFailure(
      int side, int stationIndex,
      BoundaryLayerWorkflow::MixedModeStationContext &ctx, double &ami,
      BoundaryLayerWorkflow::EdgeVelocityFallbackMode edgeMode,
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

  Context context_;
};

inline BoundaryLayerWorkflowAccess makeBoundaryLayerWorkflowAccess(
    BoundaryLayerWorkflow &workflow) {
  return BoundaryLayerWorkflowAccess({workflow.state,
                                      workflow.lattice,
                                      workflow.wgap,
                                      workflow.blc,
                                      workflow.boundaryLayerVariablesSolver,
                                      workflow.blDiffSolver,
                                      workflow.transitionSolver,
                                      workflow.flowRegime,
                                      workflow.blCompressibility,
                                      workflow.blReynolds,
                                      workflow.blTransition,
                                      workflow.nsys,
                                      workflow.xt});
}
