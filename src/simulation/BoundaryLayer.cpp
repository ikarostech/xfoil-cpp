#include "BoundaryLayer.hpp"

#include <algorithm>
#include <sstream>
#include <string>

#include "domain/boundary_layer.hpp"
#include "domain/coefficient/aero_coefficients.hpp"
#include "domain/flow_state.hpp"
#include "infrastructure/logger.hpp"
#include "simulation/boundary_layer_aerodynamics.hpp"
#include "simulation/boundary_layer_mixed_mode.hpp"
#include "simulation/boundary_layer_relaxation.hpp"
#include "simulation/boundary_layer_runtime_state.hpp"
#include "simulation/boundary_layer_setbl.hpp"
#include "simulation/boundary_layer_solver_ops.hpp"

using BoundaryContext = BoundaryLayerWorkflow::MixedModeStationContext;

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

BoundaryLayerWorkflow::EdgeVelocityDistribution
BoundaryLayerWorkflow::computeNewUeDistribution(
    const XFoil &xfoil, const Matrix3x2dVector &vdel) const {
  return BoundaryLayerAerodynamicsOps::computeNewUeDistribution(lattice, xfoil,
                                                                vdel);
}

BoundaryLayerWorkflow::ClContributions
BoundaryLayerWorkflow::computeClFromEdgeVelocityDistribution(
    const XFoil &xfoil, const EdgeVelocityDistribution &distribution) const {
  return BoundaryLayerAerodynamicsOps::computeClFromEdgeVelocityDistribution(
      lattice, xfoil, distribution);
}

int BoundaryLayerWorkflow::resetSideState(int side, const Foil &foil,
                                          const StagnationResult &stagnation) {
  return BoundaryLayerRuntimeStateOps::resetSideState(
      lattice, blTransition, flowRegime, side, foil, stagnation);
}

BoundaryLayerWorkflow::StationReadModel
BoundaryLayerWorkflow::readStationModel(int side, int stationIndex) const {
  return BoundaryLayerRuntimeStateOps::readStationModel(lattice, wgap, side,
                                                        stationIndex);
}

int BoundaryLayerWorkflow::readSideStationCount(int side) const {
  return BoundaryLayerRuntimeStateOps::readSideStationCount(lattice, side);
}

BoundaryLayerWorkflow::TrailingEdgeReadModel
BoundaryLayerWorkflow::readTrailingEdgeModel() const {
  return BoundaryLayerRuntimeStateOps::readTrailingEdgeModel(lattice);
}

FlowRegimeEnum
BoundaryLayerWorkflow::applyFlowRegimeCandidate(FlowRegimeEnum candidate) {
  flowRegime = candidate;
  return flowRegime;
}

FlowRegimeEnum BoundaryLayerWorkflow::currentFlowRegime() const {
  return flowRegime;
}

void BoundaryLayerWorkflow::emitMarchInfoLog(std::string_view message) const {
  Logger::instance().write(std::string(message));
}

double BoundaryLayerWorkflow::readNewtonRhs(int row) const {
  return blc.rhs[row];
}

void BoundaryLayerWorkflow::solveMrchueDirectNewtonSystem() {
  blc.a2(3, 0) = 0.0;
  blc.a2(3, 1) = 0.0;
  blc.a2(3, 2) = 0.0;
  blc.a2(3, 3) = 1.0;
  blc.rhs[3] = 0.0;
  blc.rhs = blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);
}

void BoundaryLayerWorkflow::solveMrchueInverseNewtonSystem(double htarg) {
  blc.a2(3, 0) = 0.0;
  blc.a2(3, 1) = state.station2.hkz.t();
  blc.a2(3, 2) = state.station2.hkz.d();
  blc.a2(3, 3) = state.station2.hkz.u();
  blc.rhs[3] = htarg - state.station2.hkz.scalar;
  blc.rhs = blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);
}

void BoundaryLayerWorkflow::storeStationStateCommon(
    int side, int stationIndex, const MixedModeStationContext &ctx) {
  makeMixedModeOps(*this).storeStationStateCommon(side, stationIndex, ctx);
}

double BoundaryLayerWorkflow::fallbackEdgeVelocity(
    int side, int stationIndex, EdgeVelocityFallbackMode edgeMode) const {
  return makeMixedModeOps(const_cast<BoundaryLayerWorkflow &>(*this))
      .fallbackEdgeVelocity(side, stationIndex, edgeMode);
}

BoundaryLayerWorkflow::BoundaryLayerDelta
BoundaryLayerWorkflow::buildBoundaryLayerDelta(
    int side, const Eigen::VectorXd &unew_side, const Eigen::VectorXd &u_ac_side,
    double dac, const Matrix3x2dVector &vdel) const {
  return BoundaryLayerRelaxationOps::buildBoundaryLayerDelta(
      lattice.get(side), unew_side, u_ac_side, dac, vdel);
}

BoundaryLayerWorkflow::BoundaryLayerMetrics
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

bool BoundaryLayerWorkflow::isStartOfWake(int side, int stationIndex) {
  return stationIndex == lattice.get(side).trailingEdgeIndex + 1;
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

void BoundaryLayerWorkflow::applySetblOutput(SetblOutputView &output) {
  blCompressibility = output.blCompressibility;
  blReynolds = output.blReynolds;
  lattice.top.profiles = std::move(output.profiles.top);
  lattice.bottom.profiles = std::move(output.profiles.bottom);
  flowRegime = output.flowRegime;
  blTransition = output.blTransition;
}

void BoundaryLayerWorkflow::checkTransitionIfNeeded(int side, int stationIndex,
                                                    bool skipCheck,
                                                    int laminarAdvance,
                                                    double &ami) {
  makeMixedModeOps(*this).checkTransitionIfNeeded(
      side, stationIndex, skipCheck, laminarAdvance, ami);
}

void BoundaryLayerWorkflow::resetStationKinematicsAfterFailure(
    int side, int stationIndex, MixedModeStationContext &ctx,
    EdgeVelocityFallbackMode edgeMode) {
  makeMixedModeOps(*this).resetStationKinematicsAfterFailure(
      side, stationIndex, ctx, edgeMode);
}

void BoundaryLayerWorkflow::recoverStationAfterFailure(
    int side, int stationIndex, MixedModeStationContext &ctx, double &ami,
    EdgeVelocityFallbackMode edgeMode, int laminarAdvance) {
  makeMixedModeOps(*this).recoverStationAfterFailure(
      side, stationIndex, ctx, ami, edgeMode, laminarAdvance);
}

SetblOutputView BoundaryLayerWorkflow::setbl(
    SidePairRef<const BoundaryLayerSideProfiles> profiles,
    const FlowState &analysis_state, const AeroCoefficients &aero_coeffs,
    double acrit, const Foil &foil, const StagnationResult &stagnation,
    const Eigen::MatrixXd &dij, bool bl_initialized) {
  return runBoundaryLayerSetbl(*this, profiles, analysis_state, aero_coeffs,
                               acrit, foil, stagnation, dij, bl_initialized);
}

SidePair<Eigen::VectorXd>
BoundaryLayerWorkflow::ueset(const Eigen::MatrixXd &dij) const {
  return BoundaryLayerAerodynamicsOps::ueset(lattice, dij);
}

/** -----------------------------------------------------
 * 	   sets forced-transition bl coordinate locations.
 * ----------------------------------------------------- */
double BoundaryLayerWorkflow::xifset(const Foil &foil,
                                     const StagnationResult &stagnation,
                                     int is) const {
  return BoundaryLayerRuntimeStateOps::xifset(lattice, foil, stagnation, is);
}
