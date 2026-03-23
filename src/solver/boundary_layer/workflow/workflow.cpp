#include "solver/boundary_layer/workflow/workflow.hpp"

#include <algorithm>

#include "model/boundary_layer.hpp"
#include "solver/boundary_layer/boundary_layer_aerodynamics.hpp"
#include "solver/boundary_layer/runtime/state.hpp"
#include "solver/boundary_layer/workflow/mixed_mode.hpp"
#include "solver/boundary_layer/workflow/relaxation.hpp"
#include "solver/boundary_layer/workflow/solver_ops.hpp"
#include "solver/march/mrchue_linear_system.hpp"

using BoundaryContext = BoundaryLayerMixedModeStationContext;

BoundaryLayerSolverOps BoundaryLayer::makeSolverOps() {
    return BoundaryLayerSolverOps({boundaryLayerVariablesSolver, blDiffSolver, transitionSolver,
                                   state_store_.flowRegime, state_store_.blCompressibility, state_store_.blReynolds,
                                   state_store_.blTransition, workspace_.state, workspace_.blc, state_store_.lattice});
}

BoundaryLayerMixedModeOps BoundaryLayer::makeMixedModeOps() {
    return BoundaryLayerMixedModeOps({state_store_.lattice, workspace_.state, state_store_.flowRegime, workspace_.blc,
                                      state_store_.blCompressibility, state_store_.blReynolds,
                                      boundaryLayerVariablesSolver, transitionSolver, makeSolverOps()});
}

void BoundaryLayer::storeStationStateCommon(int side, int stationIndex,
                                            const BoundaryLayerMixedModeStationContext &ctx) {
    makeMixedModeOps().storeStationStateCommon(side, stationIndex, ctx);
    refreshTrailingEdgeFeature();
}

double BoundaryLayer::fallbackEdgeVelocity(int side, int stationIndex,
                                           BoundaryLayerEdgeVelocityFallbackMode edgeMode) const {
    return const_cast<BoundaryLayer *>(this)->makeMixedModeOps().fallbackEdgeVelocity(side, stationIndex, edgeMode);
}

BoundaryLayerDelta BoundaryLayer::buildBoundaryLayerDelta(int side, const Eigen::VectorXd &unew_side,
                                                          const Eigen::VectorXd &u_ac_side, double dac,
                                                          const BoundaryLayerMatrix3x2dVector &vdel) const {
    return BoundaryLayerRelaxationOps::buildBoundaryLayerDelta(state_store_.lattice.get(side), unew_side, u_ac_side,
                                                               dac, vdel);
}

BoundaryLayerMetrics BoundaryLayer::evaluateSegmentRelaxation(int side, const BoundaryLayerDelta &delta, double dhi,
                                                              double dlo, double &relaxation) const {
    return BoundaryLayerRelaxationOps::evaluateSegmentRelaxation(state_store_.lattice.get(side).profiles, delta, dhi,
                                                                 dlo, relaxation);
}

BoundaryLayerSideState BoundaryLayer::applyBoundaryLayerDelta(int side, const BoundaryLayerDelta &delta,
                                                              double relaxation, double hstinv, double gamm1) const {
    return BoundaryLayerRelaxationOps::applyBoundaryLayerDelta(state_store_.lattice.get(side), state_store_.wgap, delta,
                                                               relaxation, hstinv, gamm1);
}

void BoundaryLayer::syncStationRegimeStates(int side, int stationIndex, FlowRegimeEnum stationRegime) {
    makeMixedModeOps().syncStationRegimeStates(side, stationIndex, stationRegime);
}

FlowRegimeEnum BoundaryLayer::determineRegimeForStation(int side, int stationIndex) const {
    return const_cast<BoundaryLayer *>(this)->makeMixedModeOps().determineRegimeForStation(side, stationIndex);
}

void BoundaryLayer::updateSystemMatricesForStation(const Edge &edge, int side, int stationIndex, BoundaryContext &ctx) {
    makeMixedModeOps().updateSystemMatricesForStation(edge, side, stationIndex, ctx);
}

void BoundaryLayer::initializeFirstIterationState(int side, int stationIndex, int previousTransition,
                                                  BoundaryContext &ctx, double &ueref, double &hkref) {
    makeMixedModeOps().initializeFirstIterationState(side, stationIndex, previousTransition, ctx, ueref, hkref);
}

void BoundaryLayer::configureSimilarityRow(double ueref) {
    makeMixedModeOps().configureSimilarityRow(ueref);
}

void BoundaryLayer::configureViscousRow(double hkref, double ueref, double senswt, bool resetSensitivity,
                                        bool averageSensitivity, double &sens, double &sennew) {
    makeMixedModeOps().configureViscousRow(hkref, ueref, senswt, resetSensitivity, averageSensitivity, sens, sennew);
}

bool BoundaryLayer::applyMixedModeNewtonStep(int side, int stationIndex, double &ami, BoundaryContext &ctx) {
    return makeMixedModeOps().applyMixedModeNewtonStep(side, stationIndex, ami, ctx);
}

SkinFrictionCoefficients BoundaryLayer::blmid(FlowRegimeEnum flowRegimeType) {
    return makeSolverOps().blmid(flowRegimeType);
}

bool BoundaryLayer::blsys() {
    return makeSolverOps().blsys();
}

bool BoundaryLayer::tesys(const BoundaryLayerSideState &top_profiles, const BoundaryLayerSideState &bottom_profiles,
                          const Edge &edge) {
    return makeSolverOps().tesys(top_profiles, bottom_profiles, edge);
}

void BoundaryLayer::checkTransitionIfNeeded(int side, int stationIndex, bool skipCheck, int laminarAdvance,
                                            double &ami) {
    makeMixedModeOps().checkTransitionIfNeeded(side, stationIndex, skipCheck, laminarAdvance, ami);
}

void BoundaryLayer::resetStationKinematicsAfterFailure(int side, int stationIndex,
                                                       BoundaryLayerMixedModeStationContext &ctx,
                                                       BoundaryLayerEdgeVelocityFallbackMode edgeMode) {
    makeMixedModeOps().resetStationKinematicsAfterFailure(side, stationIndex, ctx, edgeMode);
}

void BoundaryLayer::recoverStationAfterFailure(int side, int stationIndex, BoundaryLayerMixedModeStationContext &ctx,
                                               double &ami, BoundaryLayerEdgeVelocityFallbackMode edgeMode,
                                               int laminarAdvance) {
    makeMixedModeOps().recoverStationAfterFailure(side, stationIndex, ctx, ami, edgeMode, laminarAdvance);
}

void BoundaryLayer::setTransitionLocations(double top, double bottom) {
    state_store_.lattice.top.transitionLocation    = top;
    state_store_.lattice.bottom.transitionLocation = bottom;
}

double BoundaryLayer::transitionLocation(int side) const {
    return state_store_.lattice.get(side).transitionLocation;
}

void BoundaryLayer::setFlowRegime(FlowRegimeEnum flowRegime) {
    state_store_.flowRegime = flowRegime;
}

void BoundaryLayer::clearState() {
    workspace_.state = BoundaryLayerStationWindow{};
}

void BoundaryLayer::resetSideMetadata() {
    for (int side = 1; side <= 2; ++side) {
        setTransitionIndex(side, 0);
        setStationCount(side, 0);
        setTrailingEdgeIndex(side, 0);
    }
}

void BoundaryLayer::resetPhysicsState() {
    state_store_.stagnation        = BoundaryLayerStagnationFeature{};
    state_store_.trailingEdge      = BoundaryLayerTrailingEdgeFeature{};
    state_store_.flowRegime        = FlowRegimeEnum::Laminar;
    state_store_.blCompressibility = {};
    state_store_.blReynolds        = {};
    state_store_.blTransition      = {};
}

void BoundaryLayer::resetTransportState() {
    state_store_.blCompressibility.qinfbl    = 0.0;
    state_store_.blCompressibility.tkbl      = 0.0;
    state_store_.blCompressibility.tkbl_ms   = 0.0;
    state_store_.blCompressibility.rstbl     = 0.0;
    state_store_.blCompressibility.rstbl_ms  = 0.0;
    state_store_.blCompressibility.hstinv    = 0.0;
    state_store_.blCompressibility.hstinv_ms = 0.0;
    state_store_.blCompressibility.gm1bl     = 0.0;

    state_store_.blReynolds.reybl    = 0.0;
    state_store_.blReynolds.reybl_ms = 0.0;
    state_store_.blReynolds.reybl_re = 0.0;

    state_store_.blTransition.xiforc = 0.0;
    state_store_.blTransition.amcrit = 0.0;
}

void BoundaryLayer::resetForReinitialization() {
    clearState();
    resetSideMetadata();
    resetPhysicsState();
    resetTransportState();
    zeroProfiles();
}

void BoundaryLayer::setStagnationState(const StagnationResult &stagnation) {
    state_store_.stagnation.index  = stagnation.stagnationIndex;
    state_store_.stagnation.sst    = stagnation.sst;
    state_store_.stagnation.sst_go = stagnation.sst_go;
    state_store_.stagnation.sst_gp = stagnation.sst_gp;
    state_store_.stagnation.found  = stagnation.found;
}

const BoundaryLayerStagnationFeature &BoundaryLayer::stagnationFeature() const {
    return state_store_.stagnation;
}

const BoundaryLayerTrailingEdgeFeature &BoundaryLayer::trailingEdgeFeature() const {
    return state_store_.trailingEdge;
}

void BoundaryLayer::refreshTrailingEdgeFeature() {
    auto &feature = state_store_.trailingEdge;
    feature       = BoundaryLayerTrailingEdgeFeature{};

    const auto &top    = state_store_.lattice.top;
    const auto &bottom = state_store_.lattice.bottom;
    if (top.stationCount <= 0 || bottom.stationCount <= 0) {
        return;
    }
    if (top.profiles.momentumThickness.size() <= top.trailingEdgeIndex ||
        bottom.profiles.momentumThickness.size() <= bottom.trailingEdgeIndex ||
        top.profiles.displacementThickness.size() <= top.trailingEdgeIndex ||
        bottom.profiles.displacementThickness.size() <= bottom.trailingEdgeIndex ||
        top.profiles.skinFrictionCoeff.size() <= top.trailingEdgeIndex ||
        bottom.profiles.skinFrictionCoeff.size() <= bottom.trailingEdgeIndex) {
        return;
    }

    feature.topIndex                    = top.trailingEdgeIndex;
    feature.bottomIndex                 = bottom.trailingEdgeIndex;
    feature.topMomentumThickness        = top.profiles.momentumThickness[feature.topIndex];
    feature.bottomMomentumThickness     = bottom.profiles.momentumThickness[feature.bottomIndex];
    feature.topDisplacementThickness    = top.profiles.displacementThickness[feature.topIndex];
    feature.bottomDisplacementThickness = bottom.profiles.displacementThickness[feature.bottomIndex];
    feature.topSkinFrictionCoeff        = top.profiles.skinFrictionCoeff[feature.topIndex];
    feature.bottomSkinFrictionCoeff     = bottom.profiles.skinFrictionCoeff[feature.bottomIndex];
}

void BoundaryLayer::clearPanelMap() {
    setStationCount(1, 0);
    setStationCount(2, 0);
}

void BoundaryLayer::zeroProfiles() {
    for (int side = 1; side <= 2; ++side) {
        state_store_.lattice.get(side).profiles.zeroStateVectors();
    }
    refreshTrailingEdgeFeature();
}

bool BoundaryLayer::hasValidPanelMap(int total_nodes) const {
    const auto &top    = state_store_.lattice.top;
    const auto &bottom = state_store_.lattice.bottom;
    if (!top.hasValidPanelMap(total_nodes) || !bottom.hasValidPanelMap(total_nodes)) {
        return false;
    }

    return true;
}

int BoundaryLayer::trailingEdgeSystemIndex(int side) const {
    return state_store_.lattice.get(side).trailingEdgeSystemIndex();
}

BoundaryLayerEdgeVelocityDistribution
BoundaryLayer::computeNewUeDistribution(const BoundaryLayerAerodynamicContext &context,
                                        const BoundaryLayerMatrix3x2dVector &vdel) const {
    return BoundaryLayerAerodynamicsOps::computeNewUeDistribution(state_store_.lattice, context, vdel);
}

BoundaryLayerClContributions
BoundaryLayer::computeClFromEdgeVelocityDistribution(const BoundaryLayerAerodynamicContext &context,
                                                     const BoundaryLayerEdgeVelocityDistribution &distribution) const {
    return BoundaryLayerAerodynamicsOps::computeClFromEdgeVelocityDistribution(state_store_.lattice, context,
                                                                               distribution);
}

void BoundaryLayer::initializeLattices(int size) {
    const double top_transition_location           = state_store_.lattice.top.transitionLocation;
    const double bottom_transition_location        = state_store_.lattice.bottom.transitionLocation;
    state_store_.lattice.top                       = BoundaryLayerLattice(size);
    state_store_.lattice.bottom                    = BoundaryLayerLattice(size);
    state_store_.lattice.top.transitionLocation    = top_transition_location;
    state_store_.lattice.bottom.transitionLocation = bottom_transition_location;
}

void BoundaryLayer::initializeWakeGap(int wake_nodes) {
    state_store_.wgap = Eigen::VectorXd::Zero(wake_nodes);
}

void BoundaryLayer::assignInviscidEdgeVelocity(const SidePair<Eigen::Matrix2Xd> &velocity) {
    state_store_.lattice.top.inviscidEdgeVelocityMatrix    = velocity.top;
    state_store_.lattice.bottom.inviscidEdgeVelocityMatrix = velocity.bottom;
}

void BoundaryLayer::seedEdgeVelocityFromInviscid() {
    for (int side = 1; side <= 2; ++side) {
        auto &side_lattice = state_store_.lattice.get(side);
        for (int ibl = 0; ibl < side_lattice.stationCount - 1; ++ibl) {
            side_lattice.profiles.edgeVelocity[ibl] = side_lattice.inviscidEdgeVelocityMatrix(0, ibl);
        }
    }
}

void BoundaryLayer::applyProfiles(const SidePair<BoundaryLayerSideState> &profiles) {
    state_store_.lattice.top.profiles    = profiles.top;
    state_store_.lattice.bottom.profiles = profiles.bottom;
    refreshTrailingEdgeFeature();
}

void BoundaryLayer::applyProfiles(SidePair<BoundaryLayerSideState> &&profiles) {
    state_store_.lattice.top.profiles    = std::move(profiles.top);
    state_store_.lattice.bottom.profiles = std::move(profiles.bottom);
    refreshTrailingEdgeFeature();
}

StagnationResult BoundaryLayer::findStagnation(const Eigen::Matrix2Xd &surface_vortex,
                                               const Eigen::VectorXd &spline_length) const {
    return geometry_.stfind(surface_vortex, spline_length);
}

bool BoundaryLayer::buildPanelMap(int point_count, int wake_point_count) {
    return geometry_.iblpan(point_count, wake_point_count);
}

bool BoundaryLayer::rebuildArcLengthCoordinates(const Foil &foil) {
    return geometry_.xicalc(foil);
}

bool BoundaryLayer::buildSystemMapping() {
    return geometry_.iblsys(workspace_.nsys);
}

SidePair<Eigen::Matrix2Xd> BoundaryLayer::computeInviscidEdgeVelocity(const Eigen::Matrix2Xd &qinv_matrix) const {
    return geometry_.uicalc(qinv_matrix);
}

bool BoundaryLayer::moveStagnation(const Eigen::Matrix2Xd &surface_vortex, const Eigen::VectorXd &spline_length,
                                   const Foil &foil, const Eigen::Matrix2Xd &qinv_matrix,
                                   StagnationResult &stagnation) {
    const bool moved = geometry_.stmove(surface_vortex, spline_length, foil, qinv_matrix, stagnation, workspace_.nsys);
    if (moved) {
        refreshTrailingEdgeFeature();
    }
    return moved;
}

int BoundaryLayer::resetSideState(int side, const Foil &foil, const StagnationResult &stagnation) {
    return BoundaryLayerRuntimeStateOps::resetSideState(state_store_.lattice, state_store_.blTransition,
                                                        state_store_.flowRegime, side, foil, stagnation);
}

double BoundaryLayer::computeForcedTransitionArcLength(const Foil &foil, const StagnationResult &stagnation,
                                                       int side) const {
    return BoundaryLayerRuntimeStateOps::xifset(state_store_.lattice, foil, stagnation, side);
}

int BoundaryLayer::readSideStationCount(int side) const {
    return BoundaryLayerRuntimeStateOps::readSideStationCount(state_store_.lattice, side);
}

BoundaryLayerStationReadModel BoundaryLayer::readStationModel(int side, int stationIndex) const {
    return BoundaryLayerRuntimeStateOps::readStationModel(state_store_.lattice, state_store_.wgap, side, stationIndex);
}

BoundaryLayerSideReadModel BoundaryLayer::readSideModel(int side) const {
    return BoundaryLayerRuntimeStateOps::readSideModel(state_store_.lattice, side);
}

BoundaryLayerTrailingEdgeReadModel BoundaryLayer::readTrailingEdgeModel() const {
    BoundaryLayerTrailingEdgeReadModel model;
    const auto &feature               = state_store_.trailingEdge;
    model.topTrailingEdgeIndex        = feature.topIndex;
    model.bottomTrailingEdgeIndex     = feature.bottomIndex;
    model.topMomentumThickness        = feature.topMomentumThickness;
    model.bottomMomentumThickness     = feature.bottomMomentumThickness;
    model.topDisplacementThickness    = feature.topDisplacementThickness;
    model.bottomDisplacementThickness = feature.bottomDisplacementThickness;
    model.topSkinFrictionCoeff        = feature.topSkinFrictionCoeff;
    model.bottomSkinFrictionCoeff     = feature.bottomSkinFrictionCoeff;
    return model;
}

bool BoundaryLayer::isStartOfWake(int side, int stationIndex) const {
    return state_store_.lattice.get(side).isStartOfWake(stationIndex);
}

void BoundaryLayer::copyProfilesTo(SidePair<BoundaryLayerSideState> &profiles) const {
    profiles.top    = state_store_.lattice.top.profiles;
    profiles.bottom = state_store_.lattice.bottom.profiles;
}

double BoundaryLayer::inviscidEdgeVelocitySensitivityToAlpha(int side, int stationIndex) const {
    return state_store_.lattice.get(side).inviscidEdgeVelocityMatrix(1, stationIndex);
}

double BoundaryLayer::currentAmplification() const {
    return workspace_.state.station2.param.amplz;
}

double BoundaryLayer::previousAmplification() const {
    return workspace_.state.station1.param.amplz;
}

double BoundaryLayer::currentSkinFrictionHistory() const {
    return workspace_.state.station2.cqz.scalar;
}

BoundaryLayerStationWindow BoundaryLayer::snapshotState() const {
    return workspace_.state;
}

double BoundaryLayer::readCurrentShapeFactor() const {
    return workspace_.state.current().hkz.scalar;
}

void BoundaryLayer::refreshCurrentStationState(double xsi, double ami, double cti, double thi, double dsi,
                                               double dswaki, double uei) {
    BoundaryLayerPhysics::refreshCurrentStation(workspace_.state, state_store_.blCompressibility,
                                                state_store_.blReynolds, xsi, ami, cti, thi, dsi, dswaki, uei);
}

void BoundaryLayer::updateCurrentStationKinematics() {
    BoundaryLayerPhysics::blkin(workspace_.state, state_store_.blCompressibility, state_store_.blReynolds);
}

void BoundaryLayer::replaceState(const BoundaryLayerStationWindow &state) {
    workspace_.state = state;
}

void BoundaryLayer::advanceState() {
    workspace_.state.stepbl();
}

void BoundaryLayer::runTransitionCheck() {
    transitionSolver.trchek();
}

bool BoundaryLayer::solveTeSystemForCurrentProfiles(const Edge &edge) {
    return makeSolverOps().tesys(state_store_.lattice.top.profiles, state_store_.lattice.bottom.profiles, edge);
}

void BoundaryLayer::solveWakeState() {
    workspace_.state.station2 = boundaryLayerVariablesSolver.solve(workspace_.state.station2, FlowRegimeEnum::Wake);
}

SidePair<Eigen::VectorXd> BoundaryLayer::computeInviscidEdgeVelocitySensitivity(const Eigen::MatrixXd &dij) const {
    return BoundaryLayerAerodynamicsOps::ueset(state_store_.lattice, dij);
}

double BoundaryLayer::readNewtonRhs(int row) const {
    return MrchueLinearSystemOps::readNewtonRhs(workspace_.blc, row);
}

void BoundaryLayer::solveDirectNewtonSystem() {
    MrchueLinearSystemOps::solveDirect(workspace_.blc);
}

void BoundaryLayer::solveInverseNewtonSystem(double htarg) {
    MrchueLinearSystemOps::solveInverse(workspace_.blc, workspace_.state, htarg);
}

void BoundaryLayer::applyInitializationState(const BlCompressibilityParams &compressibility,
                                             const BlReynoldsParams &reynolds, const BlTransitionParams &transition,
                                             FlowRegimeEnum flowRegime, SidePair<BoundaryLayerSideState> profiles) {
    state_store_.blCompressibility = compressibility;
    state_store_.blReynolds        = reynolds;
    state_store_.blTransition      = transition;
    state_store_.flowRegime        = flowRegime;
    applyProfiles(std::move(profiles));
}

void BoundaryLayer::runTransitionCheckForMrchue(int side, int stationIndex, double &ami, double &cti,
                                                int laminarAdvance) {
    transitionSolver.trchek();
    ami = workspace_.state.station2.param.amplz;
    if (state_store_.flowRegime == FlowRegimeEnum::Transition) {
        state_store_.lattice.get(side).profiles.transitionIndex = stationIndex;
        if (cti <= 0.0) {
            cti                                = 0.03;
            workspace_.state.station2.param.sz = cti;
        }
        return;
    }

    state_store_.lattice.get(side).profiles.transitionIndex = stationIndex + laminarAdvance;
}

double BoundaryLayer::calcHtarg(int stationIndex, int side, bool wake) {
    const int transition_index = state_store_.lattice.get(side).profiles.transitionIndex;
    const auto &station1       = workspace_.state.station1;
    auto &station2             = workspace_.state.station2;

    if (stationIndex < transition_index) {
        return station1.hkz.scalar + 0.03 * (station2.param.xz - station1.param.xz) / station1.param.tz;
    }

    if (stationIndex == transition_index) {
        return station1.hkz.scalar +
               (0.03 * (workspace_.xt.scalar - station1.param.xz) - 0.15 * (station2.param.xz - workspace_.xt.scalar)) /
                   station1.param.tz;
    }

    if (wake) {
        const double cst = 0.03 * (station2.param.xz - station1.param.xz) / station1.param.tz;
        auto euler       = [](double hk2, double hk1, double cst_local) {
            return hk2 - (hk2 + cst_local * std::pow(hk2 - 1, 3) - hk1) / (1 + 3 * cst_local * std::pow(hk2 - 1, 2));
        };

        station2.hkz.scalar = station1.hkz.scalar;
        for (int i = 0; i < 3; ++i) {
            station2.hkz.scalar = euler(station2.hkz.scalar, station1.hkz.scalar, cst);
        }
        return station2.hkz.scalar;
    }

    return station1.hkz.scalar - 0.15 * (station2.param.xz - station1.param.xz) / station1.param.tz;
}
