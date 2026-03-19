#include "solver/boundary_layer/workflow/workflow.hpp"

#include <algorithm>

#include "model/boundary_layer.hpp"
#include "solver/boundary_layer/workflow/mixed_mode.hpp"
#include "solver/boundary_layer/workflow/relaxation.hpp"
#include "solver/boundary_layer/workflow/solver_ops.hpp"
#include "solver/boundary_layer/runtime/state.hpp"

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

void BoundaryLayerWorkflow::setTransitionLocations(double top, double bottom) {
  state_store_.lattice.top.transitionLocation = top;
  state_store_.lattice.bottom.transitionLocation = bottom;
}

void BoundaryLayerWorkflow::initializeLattices(int size) {
  state_store_.lattice.top = BoundaryLayerLattice(size);
  state_store_.lattice.bottom = BoundaryLayerLattice(size);
}

void BoundaryLayerWorkflow::initializeWakeGap(int wake_nodes) {
  state_store_.wgap = Eigen::VectorXd::Zero(wake_nodes);
}

void BoundaryLayerWorkflow::assignInviscidEdgeVelocity(
    const SidePair<Eigen::Matrix2Xd> &velocity) {
  state_store_.lattice.top.inviscidEdgeVelocityMatrix = velocity.top;
  state_store_.lattice.bottom.inviscidEdgeVelocityMatrix = velocity.bottom;
}

void BoundaryLayerWorkflow::seedEdgeVelocityFromInviscid() {
  for (int side = 1; side <= 2; ++side) {
    auto &side_lattice = state_store_.lattice.get(side);
    for (int ibl = 0; ibl < side_lattice.stationCount - 1; ++ibl) {
      side_lattice.profiles.edgeVelocity[ibl] =
          side_lattice.inviscidEdgeVelocityMatrix(0, ibl);
    }
  }
}

void BoundaryLayerWorkflow::applyProfiles(
    const SidePair<BoundaryLayerSideProfiles> &profiles) {
  state_store_.lattice.top.profiles = profiles.top;
  state_store_.lattice.bottom.profiles = profiles.bottom;
}

void BoundaryLayerWorkflow::applyProfiles(
    SidePair<BoundaryLayerSideProfiles> &&profiles) {
  state_store_.lattice.top.profiles = std::move(profiles.top);
  state_store_.lattice.bottom.profiles = std::move(profiles.bottom);
}

StagnationResult BoundaryLayerWorkflow::findStagnation(
    const Eigen::Matrix2Xd &surface_vortex,
    const Eigen::VectorXd &spline_length) const {
  return geometry_.stfind(surface_vortex, spline_length);
}

bool BoundaryLayerWorkflow::buildPanelMap(int point_count, int wake_point_count) {
  return geometry_.iblpan(point_count, wake_point_count);
}

bool BoundaryLayerWorkflow::rebuildArcLengthCoordinates(const Foil &foil) {
  return geometry_.xicalc(foil);
}

bool BoundaryLayerWorkflow::buildSystemMapping() {
  return geometry_.iblsys(workspace_.nsys);
}

SidePair<Eigen::Matrix2Xd>
BoundaryLayerWorkflow::computeInviscidEdgeVelocity(
    const Eigen::Matrix2Xd &qinv_matrix) const {
  return geometry_.uicalc(qinv_matrix);
}

bool BoundaryLayerWorkflow::moveStagnation(
    const Eigen::Matrix2Xd &surface_vortex,
    const Eigen::VectorXd &spline_length, const Foil &foil,
    const Eigen::Matrix2Xd &qinv_matrix, StagnationResult &stagnation) {
  return geometry_.stmove(surface_vortex, spline_length, foil, qinv_matrix,
                          stagnation, workspace_.nsys);
}

int BoundaryLayerWorkflow::resetSideState(int side, const Foil &foil,
                                          const StagnationResult &stagnation) {
  return BoundaryLayerRuntimeStateOps::resetSideState(
      state_store_.lattice, state_store_.blTransition, state_store_.flowRegime,
      side, foil, stagnation);
}

double BoundaryLayerWorkflow::computeForcedTransitionArcLength(
    const Foil &foil, const StagnationResult &stagnation, int side) const {
  return BoundaryLayerRuntimeStateOps::xifset(state_store_.lattice, foil,
                                              stagnation, side);
}

void BoundaryLayerWorkflow::runTransitionCheckForMrchue(
    int side, int stationIndex, double &ami, double &cti, int laminarAdvance) {
  transitionSolver.trchek();
  ami = workspace_.state.station2.param.amplz;
  if (state_store_.flowRegime == FlowRegimeEnum::Transition) {
    state_store_.lattice.get(side).profiles.transitionIndex = stationIndex;
    if (cti <= 0.0) {
      cti = 0.03;
      workspace_.state.station2.param.sz = cti;
    }
    return;
  }

  state_store_.lattice.get(side).profiles.transitionIndex =
      stationIndex + laminarAdvance;
}

double BoundaryLayerWorkflow::calcHtarg(int stationIndex, int side, bool wake) {
  const int transition_index =
      state_store_.lattice.get(side).profiles.transitionIndex;
  const auto &station1 = workspace_.state.station1;
  auto &station2 = workspace_.state.station2;

  if (stationIndex < transition_index) {
    return station1.hkz.scalar +
           0.03 * (station2.param.xz - station1.param.xz) / station1.param.tz;
  }

  if (stationIndex == transition_index) {
    return station1.hkz.scalar +
           (0.03 * (workspace_.xt.scalar - station1.param.xz) -
            0.15 * (station2.param.xz - workspace_.xt.scalar)) /
               station1.param.tz;
  }

  if (wake) {
    const double cst =
        0.03 * (station2.param.xz - station1.param.xz) / station1.param.tz;
    auto euler = [](double hk2, double hk1, double cst_local) {
      return hk2 - (hk2 + cst_local * std::pow(hk2 - 1, 3) - hk1) /
                       (1 + 3 * cst_local * std::pow(hk2 - 1, 2));
    };

    station2.hkz.scalar = station1.hkz.scalar;
    for (int i = 0; i < 3; ++i) {
      station2.hkz.scalar = euler(station2.hkz.scalar, station1.hkz.scalar, cst);
    }
    return station2.hkz.scalar;
  }

  return station1.hkz.scalar -
         0.15 * (station2.param.xz - station1.param.xz) / station1.param.tz;
}
