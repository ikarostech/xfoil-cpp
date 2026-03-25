#pragma once

#include <array>
#include <cmath>

#include "model/analysis/stagnation_feature.hpp"
#include "model/boundary_layer/reference/bl_compressibility_params.hpp"
#include "model/boundary_layer/reference/bl_reynolds_params.hpp"
#include "model/boundary_layer/reference/bl_transition_params.hpp"
#include "model/boundary_layer/features.hpp"
#include "model/boundary_layer/physics.hpp"
#include "model/boundary_layer/skin_friction_coefficients.hpp"
#include "model/boundary_layer/state.hpp"
#include "model/flow_regime.hpp"
#include "numerics/boundary_layer/diff_system.hpp"
#include "numerics/boundary_layer/variables.hpp"
#include "numerics/coefficient/bl_newton.hpp"
#include "solver/boundary_layer/boundary_layer_aerodynamics.hpp"
#include "solver/boundary_layer/boundary_layer_geometry.hpp"
#include "solver/boundary_layer/runtime/state.hpp"
#include "solver/boundary_layer/viscous_types.hpp"
#include "solver/boundary_layer/workflow/transition.hpp"

class Edge;
class BoundaryLayer;
class BoundaryLayerMarchAccess;
class BoundaryLayerSetblAccess;
class BoundaryLayerSolverOps;
class BoundaryLayerMixedModeOps;

struct BoundaryLayerStationWindowStore {
    Eigen::VectorXd wgap;
    SidePair<BoundaryLayerLattice> lattice;
    StagnationFeature stagnation;
    BoundaryLayerTrailingEdgeFeature trailingEdge;
    FlowRegimeEnum flowRegime = FlowRegimeEnum::Laminar;
    BlCompressibilityParams blCompressibility{};
    BlReynoldsParams blReynolds{};
    BlTransitionParams blTransition{};
};

struct BoundaryLayerWorkspace {
    BlSystemCoeffs blc;
    BoundaryLayerStationWindow state;
    blDiff xt;
    int nsys = 0;
};

class BoundaryLayer {
  public:
    BoundaryLayer();

    BoundaryLayerVariablesSolver boundaryLayerVariablesSolver;
    BlDiffSolver blDiffSolver;
    BoundaryLayerTransitionSolver transitionSolver;

    struct QtanResult {
        Eigen::VectorXd qnew;
        Eigen::VectorXd q_ac;
    };
    void updateSystemMatricesForStation(const Edge &edge, int side, int stationIndex,
                                        BoundaryLayerMixedModeStationContext &ctx);
    void initializeFirstIterationState(int side, int stationIndex, int previousTransition,
                                       BoundaryLayerMixedModeStationContext &ctx, double &ueref, double &hkref);
    void configureSimilarityRow(double ueref);
    void configureViscousRow(double hkref, double ueref, double senswt, bool resetSensitivity, bool averageSensitivity,
                             double &sens, double &sennew);
    bool applyMixedModeNewtonStep(int side, int stationIndex, double &ami, BoundaryLayerMixedModeStationContext &ctx);
    void checkTransitionIfNeeded(int side, int stationIndex, bool skipCheck, int laminarAdvance, double &ami);
    bool blsys();
    SkinFrictionCoefficients blmid(FlowRegimeEnum flowRegimeType);
    void storeStationStateCommon(int side, int stationIndex, const BoundaryLayerMixedModeStationContext &ctx);
    void syncStationRegimeStates(int side, int stationIndex, FlowRegimeEnum stationRegime);
    FlowRegimeEnum determineRegimeForStation(int side, int stationIndex) const;
    double fallbackEdgeVelocity(int side, int stationIndex, BoundaryLayerEdgeVelocityFallbackMode edgeMode) const;
    void resetStationKinematicsAfterFailure(int side, int stationIndex, BoundaryLayerMixedModeStationContext &ctx,
                                            BoundaryLayerEdgeVelocityFallbackMode edgeMode);
    void recoverStationAfterFailure(int side, int stationIndex, BoundaryLayerMixedModeStationContext &ctx, double &ami,
                                    BoundaryLayerEdgeVelocityFallbackMode edgeMode, int laminarAdvance);
    BoundaryLayerDelta buildBoundaryLayerDelta(int side, const Eigen::VectorXd &unew_side,
                                               const Eigen::VectorXd &u_ac_side, double dac,
                                               const BoundaryLayerMatrix3x2dVector &vdel) const;
    BoundaryLayerMetrics evaluateSegmentRelaxation(int side, const BoundaryLayerDelta &delta, double dhi, double dlo,
                                                   double &relaxation) const;
    BoundaryLayerSideState applyBoundaryLayerDelta(int side, const BoundaryLayerDelta &delta, double relaxation,
                                                   double hstinv, double gamm1) const;

    bool tesys(const BoundaryLayerSideState &top_profiles, const BoundaryLayerSideState &bottom_profiles,
               const Edge &edge);

    int trailingEdgeIndex(int side) const {
        return state_store_.lattice.get(side).trailingEdgeIndex;
    }
    void setTrailingEdgeIndex(int side, int index) {
        state_store_.lattice.get(side).trailingEdgeIndex = index;
    }
    int stationToSystem(int side, int stationIndex) const {
        return state_store_.lattice.get(side).stationToSystem[stationIndex];
    }
    int stationToPanel(int side, int stationIndex) const {
        return state_store_.lattice.get(side).stationToPanel[stationIndex];
    }
    double panelInfluenceFactor(int side, int stationIndex) const {
        return state_store_.lattice.get(side).panelInfluenceFactor[stationIndex];
    }
    SidePairRef<const BoundaryLayerSideState> currentProfiles() const {
        return {state_store_.lattice.top.profiles, state_store_.lattice.bottom.profiles};
    }
    void clearState();
    void resetSideMetadata();
    void resetPhysicsState();
    void resetTransportState();
    void resetForReinitialization();
    void setStagnationState(const StagnationFeature &stagnation);
    const StagnationFeature &stagnationFeature() const;
    const BoundaryLayerTrailingEdgeFeature &trailingEdgeFeature() const;
    void refreshTrailingEdgeFeature();
    void clearPanelMap();
    void zeroProfiles();
    bool hasValidPanelMap(int total_nodes) const;
    int trailingEdgeSystemIndex(int side) const;
    BoundaryLayerEdgeVelocityDistribution computeNewUeDistribution(const BoundaryLayerAerodynamicContext &context,
                                                                   const BoundaryLayerMatrix3x2dVector &vdel) const;
    BoundaryLayerClContributions
    computeClFromEdgeVelocityDistribution(const BoundaryLayerAerodynamicContext &context,
                                          const BoundaryLayerEdgeVelocityDistribution &distribution) const;

    Eigen::VectorXd &wakeGap() {
        return state_store_.wgap;
    }
    const Eigen::VectorXd &wakeGap() const {
        return state_store_.wgap;
    }
    BlSystemCoeffs &systemCoefficients() {
        return workspace_.blc;
    }
    const BlSystemCoeffs &systemCoefficients() const {
        return workspace_.blc;
    }
    int &systemSize() {
        return workspace_.nsys;
    }
    int systemSize() const {
        return workspace_.nsys;
    }
    blDiff &transitionSensitivity() {
        return workspace_.xt;
    }
    const blDiff &transitionSensitivity() const {
        return workspace_.xt;
    }
    FlowRegimeEnum &flowRegime() {
        return state_store_.flowRegime;
    }
    const FlowRegimeEnum &flowRegime() const {
        return state_store_.flowRegime;
    }
    void setFlowRegime(FlowRegimeEnum flowRegime);
    BlCompressibilityParams &compressibility() {
        return state_store_.blCompressibility;
    }
    const BlCompressibilityParams &compressibility() const {
        return state_store_.blCompressibility;
    }
    BlReynoldsParams &reynolds() {
        return state_store_.blReynolds;
    }
    const BlReynoldsParams &reynolds() const {
        return state_store_.blReynolds;
    }
    BlTransitionParams &transition() {
        return state_store_.blTransition;
    }
    const BlTransitionParams &transition() const {
        return state_store_.blTransition;
    }
    void setTransitionLocations(double top, double bottom);
    double transitionLocation(int side) const;
    void initializeLattices(int size);
    void initializeWakeGap(int wake_nodes);
    void assignInviscidEdgeVelocity(const SidePair<Eigen::Matrix2Xd> &velocity);
    void seedEdgeVelocityFromInviscid();
    void applyProfiles(const SidePair<BoundaryLayerSideState> &profiles);
    void applyProfiles(SidePair<BoundaryLayerSideState> &&profiles);
    StagnationFeature findStagnation(const Eigen::Matrix2Xd &surface_vortex, const Eigen::VectorXd &spline_length) const;
    bool buildPanelMap(int point_count, int wake_point_count);
    bool rebuildArcLengthCoordinates(const Foil &foil);
    bool buildSystemMapping();
    SidePair<Eigen::Matrix2Xd> computeInviscidEdgeVelocity(const Eigen::Matrix2Xd &qinv_matrix) const;
    bool moveStagnation(const Eigen::Matrix2Xd &surface_vortex, const Eigen::VectorXd &spline_length, const Foil &foil,
                        const Eigen::Matrix2Xd &qinv_matrix, StagnationFeature &stagnation);
    int resetSideState(int side, const Foil &foil, const StagnationFeature &stagnation);
    double computeForcedTransitionArcLength(const Foil &foil, const StagnationFeature &stagnation, int side) const;
    int readSideStationCount(int side) const;
    BoundaryLayerStationReadModel readStationModel(int side, int stationIndex) const;
    BoundaryLayerSideReadModel readSideModel(int side) const;
    BoundaryLayerTrailingEdgeReadModel readTrailingEdgeModel() const;
    bool isStartOfWake(int side, int stationIndex) const;
    void copyProfilesTo(SidePair<BoundaryLayerSideState> &profiles) const;
    double inviscidEdgeVelocitySensitivityToAlpha(int side, int stationIndex) const;
    double currentAmplification() const;
    double previousAmplification() const;
    double currentSkinFrictionHistory() const;
    BoundaryLayerStationWindow snapshotState() const;
    double readCurrentShapeFactor() const;
    void refreshCurrentStationState(double xsi, double ami, double cti, double thi, double dsi, double dswaki,
                                    double uei);
    void updateCurrentStationKinematics();
    void replaceState(const BoundaryLayerStationWindow &state);
    void advanceState();
    void runTransitionCheck();
    bool solveTeSystemForCurrentProfiles(const Edge &edge);
    void solveWakeState();
    SidePair<Eigen::VectorXd> computeInviscidEdgeVelocitySensitivity(const Eigen::MatrixXd &dij) const;
    double readNewtonRhs(int row) const;
    void solveDirectNewtonSystem();
    void solveInverseNewtonSystem(double htarg);
    void applyInitializationState(const BlCompressibilityParams &compressibility, const BlReynoldsParams &reynolds,
                                  const BlTransitionParams &transition, FlowRegimeEnum flowRegime,
                                  SidePair<BoundaryLayerSideState> profiles);
    void runTransitionCheckForMrchue(int side, int stationIndex, double &ami, double &cti, int laminarAdvance = 2);
    double calcHtarg(int stationIndex, int side, bool wake);

  private:
    friend class BoundaryLayerMarchAccess;
    friend class BoundaryLayerSetblAccess;

    int stationCount(int side) const {
        return state_store_.lattice.get(side).stationCount;
    }
    void setStationCount(int side, int count) {
        state_store_.lattice.get(side).stationCount = count;
    }
    int transitionIndex(int side) const {
        return state_store_.lattice.get(side).profiles.transitionIndex;
    }
    void setTransitionIndex(int side, int index) {
        state_store_.lattice.get(side).profiles.transitionIndex = index;
    }
    const BoundaryLayerSideState &profiles(int side) const {
        return state_store_.lattice.get(side).profiles;
    }
    BoundaryLayerSideState &profiles(int side) {
        return state_store_.lattice.get(side).profiles;
    }

    SidePair<BoundaryLayerLattice> &lattice() {
        return state_store_.lattice;
    }
    const SidePair<BoundaryLayerLattice> &lattice() const {
        return state_store_.lattice;
    }
    BoundaryLayerLattice &lattice(int side) {
        return state_store_.lattice.get(side);
    }
    const BoundaryLayerLattice &lattice(int side) const {
        return state_store_.lattice.get(side);
    }
    BoundaryLayerStationWindow &state() {
        return workspace_.state;
    }
    const BoundaryLayerStationWindow &state() const {
        return workspace_.state;
    }
    BoundaryLayerSolverOps makeSolverOps();
    BoundaryLayerMixedModeOps makeMixedModeOps();

    BoundaryLayerStationWindowStore state_store_;
    BoundaryLayerWorkspace workspace_;
    BoundaryLayerGeometry geometry_;
};
