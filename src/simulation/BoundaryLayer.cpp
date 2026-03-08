#include "BoundaryLayer.hpp"

#include <algorithm>

#include "domain/boundary_layer.hpp"
#include "simulation/boundary_layer_mixed_mode.hpp"
#include "simulation/boundary_layer_relaxation.hpp"
#include "simulation/boundary_layer_solver_ops.hpp"

using BoundaryContext = BoundaryLayerMixedModeStationContext;

namespace {
BoundaryLayerSolverOps makeSolverOps(BoundaryLayerWorkflow &workflow) {
  return BoundaryLayerSolverOps({workflow.boundaryLayerVariablesSolver,
                                 workflow.blDiffSolver,
                                 workflow.transitionSolver,
                                 workflow.flowRegime,
                                 workflow.blCompressibility,
                                 workflow.blReynolds,
                                 workflow.blTransition,
                                 workflow.state,
                                 workflow.blc,
                                 workflow.lattice});
}

BoundaryLayerMixedModeOps makeMixedModeOps(BoundaryLayerWorkflow &workflow) {
  return BoundaryLayerMixedModeOps({workflow.lattice,
                                    workflow.state,
                                    workflow.flowRegime,
                                    workflow.blc,
                                    workflow.blCompressibility,
                                    workflow.boundaryLayerVariablesSolver,
                                    workflow.transitionSolver,
                                    makeSolverOps(workflow)});
}
} // namespace

double BoundaryLayerWorkflow::adjustDisplacementForHkLimit(
    double displacementThickness, double momentumThickness, double msq,
    double hklim) {
  const double h = displacementThickness / momentumThickness;

  boundary_layer::KineticShapeParameterResult hkin_result =
      boundary_layer::hkin(h, msq);

  const double dh = std::max(0.0, hklim - hkin_result.hk) / hkin_result.hk_h;
  return displacementThickness + dh * momentumThickness;
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
  return BoundaryLayerRelaxationOps::buildBoundaryLayerDelta(
      lattice.get(side), unew_side, u_ac_side, dac, vdel);
}

BoundaryLayerMetrics
BoundaryLayerWorkflow::evaluateSegmentRelaxation(
    int side, const BoundaryLayerDelta &delta, double dhi, double dlo,
    double &relaxation) const {
  return BoundaryLayerRelaxationOps::evaluateSegmentRelaxation(
      lattice.get(side).profiles, delta, dhi, dlo, relaxation);
}

BoundaryLayerSideProfiles BoundaryLayerWorkflow::applyBoundaryLayerDelta(
    int side, const BoundaryLayerDelta &delta, double relaxation, double hstinv,
    double gamm1) const {
  return BoundaryLayerRelaxationOps::applyBoundaryLayerDelta(
      lattice.get(side), wgap, delta, relaxation, hstinv, gamm1);
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
