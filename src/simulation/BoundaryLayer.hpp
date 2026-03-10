#pragma once

#include <array>
#include <cmath>

#include "domain/boundary_layer/bl_compressibility_params.hpp"
#include "domain/boundary_layer/bl_reynolds_params.hpp"
#include "domain/boundary_layer/bl_transition_params.hpp"
#include "domain/boundary_layer/boundary_layer_diff_solver.hpp"
#include "domain/boundary_layer/boundary_layer_variables_solver.hpp"
#include "domain/coefficient/bl_newton.hpp"
#include "domain/flow_regime.hpp"
#include "simulation/BoundaryLayer_transition.hpp"
#include "simulation/boundary_layer_geometry.hpp"
#include "simulation/boundary_layer_state.hpp"
#include "simulation/viscous_types.hpp"
#include "simulation/skin_friction_coefficients.hpp"

class Edge;

struct BoundaryLayerStateStore {
    Eigen::VectorXd wgap;
    SidePair<BoundaryLayerLattice> lattice;
    FlowRegimeEnum flowRegime = FlowRegimeEnum::Laminar;
    BlCompressibilityParams blCompressibility{};
    BlReynoldsParams blReynolds{};
    BlTransitionParams blTransition{};
    int stagnationIndex  = 0;
    double stagnationSst = 0.0;
};

struct BoundaryLayerWorkspace {
    BlSystemCoeffs blc;
    BoundaryLayerState state;
    blDiff xt;
    int nsys = 0;
};

class BoundaryLayerWorkflow {
  public:
    BoundaryLayerWorkflow();

    BoundaryLayerVariablesSolver boundaryLayerVariablesSolver;
    BlDiffSolver blDiffSolver;
    BoundaryLayerTransitionSolver transitionSolver;

    // Sutherland's const./T0 (assumes stagnation conditions are at STP).
    static constexpr double kHvrat = 0.35;

    struct QtanResult {
        Eigen::VectorXd qnew;
        Eigen::VectorXd q_ac;
    };
    void updateSystemMatricesForStation(const Edge &edge, int side, int stationIndex, BoundaryLayerMixedModeStationContext &ctx);
    void initializeFirstIterationState(int side, int stationIndex, int previousTransition, BoundaryLayerMixedModeStationContext &ctx,
                                       double &ueref, double &hkref);
    void configureSimilarityRow(double ueref);
    void configureViscousRow(double hkref, double ueref, double senswt, bool resetSensitivity, bool averageSensitivity,
                             double &sens, double &sennew);
    bool applyMixedModeNewtonStep(int side, int stationIndex, double &ami, BoundaryLayerMixedModeStationContext &ctx);
    void checkTransitionIfNeeded(int side, int stationIndex, bool skipCheck, int laminarAdvance, double &ami);
    SkinFrictionCoefficients blmid(FlowRegimeEnum flowRegimeType);
    blData blprv(blData data, double xsi, double ami, double cti, double thi, double dsi, double dswaki,
                 double uei) const;
    bool blsys();
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
    BoundaryLayerSideProfiles applyBoundaryLayerDelta(int side, const BoundaryLayerDelta &delta, double relaxation,
                                                      double hstinv, double gamm1) const;

    bool blkin(BoundaryLayerState &state);
    bool tesys(const BoundaryLayerSideProfiles &top_profiles, const BoundaryLayerSideProfiles &bottom_profiles,
               const Edge &edge);

    BoundaryLayerStateStore &stateStore() {
        return state_store_;
    }
    const BoundaryLayerStateStore &stateStore() const {
        return state_store_;
    }
    BoundaryLayerWorkspace &workspace() {
        return workspace_;
    }
    const BoundaryLayerWorkspace &workspace() const {
        return workspace_;
    }
    BoundaryLayerGeometry &geometryService() {
        return geometry_;
    }
    const BoundaryLayerGeometry &geometryService() const {
        return geometry_;
    }

  private:
    static double adjustDisplacementForHkLimit(double displacementThickness, double momentumThickness, double msq,
                                               double hklim);

    BoundaryLayerStateStore state_store_;
    BoundaryLayerWorkspace workspace_;
    BoundaryLayerGeometry geometry_;
};
