#include "solver/boundary_layer/workflow/workflow.hpp"

#include <algorithm>

#include "model/boundary_layer.hpp"
#include "solver/boundary_layer/workflow/mixed_mode.hpp"
#include "solver/boundary_layer/workflow/relaxation.hpp"
#include "solver/boundary_layer/workflow/solver_ops.hpp"

using BoundaryContext = BoundaryLayerMixedModeStationContext;

namespace {
BoundaryLayerSolverOps makeSolverOps(BoundaryLayerWorkflow &workflow) {
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

BoundaryLayerMixedModeOps makeMixedModeOps(BoundaryLayerWorkflow &workflow) {
  auto &state_store = workflow.stateStore();
  auto &workspace = workflow.workspace();
  return BoundaryLayerMixedModeOps({state_store.lattice,
                                    workspace.state,
                                    state_store.flowRegime,
                                    workspace.blc,
                                    state_store.blCompressibility,
                                    workflow.boundaryLayerVariablesSolver,
                                    workflow.transitionSolver,
                                    makeSolverOps(workflow)});
}
}

void BoundaryLayerWorkflow::storeStationStateCommon(
    int side, int stationIndex,
    const BoundaryLayerMixedModeStationContext &ctx) {
  makeMixedModeOps(*this).storeStationStateCommon(side, stationIndex, ctx);
}

double BoundaryLayerWorkflow::fallbackEdgeVelocity(
    int side, int stationIndex,
    BoundaryLayerEdgeVelocityFallbackMode edgeMode) const {
  return makeMixedModeOps(const_cast<BoundaryLayerWorkflow &>(*this))
      .fallbackEdgeVelocity(side, stationIndex, edgeMode);
}

BoundaryLayerDelta
BoundaryLayerWorkflow::buildBoundaryLayerDelta(
    int side, const Eigen::VectorXd &unew_side, const Eigen::VectorXd &u_ac_side,
    double dac, const BoundaryLayerMatrix3x2dVector &vdel) const {
  const auto &state_store = stateStore();
  return BoundaryLayerRelaxationOps::buildBoundaryLayerDelta(
      state_store.lattice.get(side), unew_side, u_ac_side, dac, vdel);
}

BoundaryLayerMetrics
BoundaryLayerWorkflow::evaluateSegmentRelaxation(
    int side, const BoundaryLayerDelta &delta, double dhi, double dlo,
    double &relaxation) const {
  const auto &state_store = stateStore();
  return BoundaryLayerRelaxationOps::evaluateSegmentRelaxation(
      state_store.lattice.get(side).profiles, delta, dhi, dlo, relaxation);
}

BoundaryLayerSideProfiles BoundaryLayerWorkflow::applyBoundaryLayerDelta(
    int side, const BoundaryLayerDelta &delta, double relaxation, double hstinv,
    double gamm1) const {
  const auto &state_store = stateStore();
  return BoundaryLayerRelaxationOps::applyBoundaryLayerDelta(
      state_store.lattice.get(side), state_store.wgap, delta, relaxation,
      hstinv, gamm1);
}

void BoundaryLayerWorkflow::syncStationRegimeStates(int side, int stationIndex,
                                                    FlowRegimeEnum stationRegime) {
  makeMixedModeOps(*this).syncStationRegimeStates(side, stationIndex,
                                                  stationRegime);
}

FlowRegimeEnum BoundaryLayerWorkflow::determineRegimeForStation(
    int side, int stationIndex) const {
  return makeMixedModeOps(const_cast<BoundaryLayerWorkflow &>(*this))
      .determineRegimeForStation(side, stationIndex);
}

bool BoundaryLayerWorkflow::blkin(BoundaryLayerState &state) {
  return makeSolverOps(*this).blkin(state);
}

void BoundaryLayerWorkflow::updateSystemMatricesForStation(
    const Edge &edge, int side, int stationIndex, BoundaryContext &ctx) {
  makeMixedModeOps(*this).updateSystemMatricesForStation(edge, side,
                                                         stationIndex, ctx);
}

void BoundaryLayerWorkflow::initializeFirstIterationState(
    int side, int stationIndex, int previousTransition, BoundaryContext &ctx,
    double &ueref, double &hkref) {
  makeMixedModeOps(*this).initializeFirstIterationState(
      side, stationIndex, previousTransition, ctx, ueref, hkref);
}

void BoundaryLayerWorkflow::configureSimilarityRow(double ueref) {
  makeMixedModeOps(*this).configureSimilarityRow(ueref);
}

void BoundaryLayerWorkflow::configureViscousRow(double hkref, double ueref,
                                                double senswt,
                                                bool resetSensitivity,
                                                bool averageSensitivity,
                                                double &sens, double &sennew) {
  makeMixedModeOps(*this).configureViscousRow(
      hkref, ueref, senswt, resetSensitivity, averageSensitivity, sens,
      sennew);
}

bool BoundaryLayerWorkflow::applyMixedModeNewtonStep(int side, int stationIndex,
                                                     double &ami,
                                                     BoundaryContext &ctx) {
  return makeMixedModeOps(*this).applyMixedModeNewtonStep(
      side, stationIndex, ami, ctx);
}

SkinFrictionCoefficients
BoundaryLayerWorkflow::blmid(FlowRegimeEnum flowRegimeType) {
  return makeSolverOps(*this).blmid(flowRegimeType);
}

blData BoundaryLayerWorkflow::blprv(blData data, double xsi, double ami,
                                    double cti, double thi, double dsi,
                                    double dswaki, double uei) const {
  return makeSolverOps(const_cast<BoundaryLayerWorkflow &>(*this))
      .blprv(data, xsi, ami, cti, thi, dsi, dswaki, uei);
}

bool BoundaryLayerWorkflow::blsys() {
  return makeSolverOps(*this).blsys();
}

bool BoundaryLayerWorkflow::tesys(
    const BoundaryLayerSideProfiles &top_profiles,
    const BoundaryLayerSideProfiles &bottom_profiles, const Edge &edge) {
  return makeSolverOps(*this).tesys(top_profiles, bottom_profiles, edge);
}

void BoundaryLayerWorkflow::checkTransitionIfNeeded(int side, int stationIndex,
                                                    bool skipCheck,
                                                    int laminarAdvance,
                                                    double &ami) {
  makeMixedModeOps(*this).checkTransitionIfNeeded(
      side, stationIndex, skipCheck, laminarAdvance, ami);
}

void BoundaryLayerWorkflow::resetStationKinematicsAfterFailure(
    int side, int stationIndex, BoundaryLayerMixedModeStationContext &ctx,
    BoundaryLayerEdgeVelocityFallbackMode edgeMode) {
  makeMixedModeOps(*this).resetStationKinematicsAfterFailure(
      side, stationIndex, ctx, edgeMode);
}

void BoundaryLayerWorkflow::recoverStationAfterFailure(
    int side, int stationIndex, BoundaryLayerMixedModeStationContext &ctx,
    double &ami, BoundaryLayerEdgeVelocityFallbackMode edgeMode,
    int laminarAdvance) {
  makeMixedModeOps(*this).recoverStationAfterFailure(
      side, stationIndex, ctx, ami, edgeMode, laminarAdvance);
}
