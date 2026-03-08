#pragma once

#include <array>
#include <cmath>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

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
#include "simulation/skin_friction_coefficients.hpp"

struct FlowState;
struct AeroCoefficients;
struct SetblOutputView;
class Foil;
class Edge;
class XFoil;
class BoundaryLayerWorkflow {
  public:
    BoundaryLayerWorkflow();
    using Matrix3x2d       = Eigen::Matrix<double, 3, 2>;
    using Matrix3x2dVector = std::vector<Matrix3x2d>;

    BoundaryLayerVariablesSolver boundaryLayerVariablesSolver;
    BlDiffSolver blDiffSolver;
    BoundaryLayerTransitionSolver transitionSolver;
    FlowRegimeEnum flowRegime = FlowRegimeEnum::Laminar;
    BlCompressibilityParams blCompressibility{};
    BlReynoldsParams blReynolds{};
    BlTransitionParams blTransition{};

    // Sutherland's const./T0 (assumes stagnation conditions are at STP).
    static constexpr double kHvrat = 0.35;

    enum class EdgeVelocityFallbackMode { UsePreviousStation, AverageNeighbors };
    struct EdgeVelocityDistribution {
        SidePair<Eigen::VectorXd> unew;
        SidePair<Eigen::VectorXd> u_ac;
    };
    struct QtanResult {
        Eigen::VectorXd qnew;
        Eigen::VectorXd q_ac;
    };
    struct ClContributions {
        double cl    = 0.0;
        double cl_a  = 0.0;
        double cl_ms = 0.0;
        double cl_ac = 0.0;
    };
    struct BoundaryLayerDelta {
        Eigen::VectorXd dskinFrictionCoeff;
        Eigen::VectorXd dmomentumThickness;
        Eigen::VectorXd ddisplacementThickness;
        Eigen::VectorXd dedgeVelocity;
    };
    struct BoundaryLayerMetrics {
        double rmsContribution = 0.0;
        double maxChange       = 0.0;
    };
    struct MixedModeStationContext {
        FlowRegimeEnum flowRegime = FlowRegimeEnum::Laminar;
        double xsi                = 0.0;
        double uei                = 0.0;
        double thi                = 0.0;
        double dsi                = 0.0;
        double cti                = 0.0;
        double ami                = 0.0;
        double dswaki             = 0.0;
        double cte                = 0.0;
        double dte                = 0.0;
        double tte                = 0.0;
        double dmax               = 0.0;

        bool isSimilarity() const {
            return flowRegime == FlowRegimeEnum::Similarity;
        }

        bool isWake() const {
            return flowRegime == FlowRegimeEnum::Wake;
        }
    };
    struct StationReadModel {
        int stationCount      = 0;
        int trailingEdgeIndex = 0;
        int transitionIndex   = 0;
        double arcLength      = 0.0;
        double edgeVelocity   = 0.0;
        double momentumThickness = 0.0;
        double displacementThickness = 0.0;
        double skinFrictionCoeff    = 0.0;
        double wakeGap              = 0.0;
    };
    struct TrailingEdgeReadModel {
        int topTrailingEdgeIndex    = 0;
        int bottomTrailingEdgeIndex = 0;
        double topMomentumThickness    = 0.0;
        double bottomMomentumThickness = 0.0;
        double topDisplacementThickness    = 0.0;
        double bottomDisplacementThickness = 0.0;
        double topSkinFrictionCoeff    = 0.0;
        double bottomSkinFrictionCoeff = 0.0;
    };

    int readSideStationCount(int side) const;
    StationReadModel readStationModel(int side, int stationIndex) const;
    TrailingEdgeReadModel readTrailingEdgeModel() const;
    void emitMarchInfoLog(std::string_view message) const;
    void emitMarchFailureLog(std::string_view phase, int side, int stationIndex, double residual) const;
    double readNewtonRhs(int row) const;
    void solveMrchueDirectNewtonSystem();
    void solveMrchueInverseNewtonSystem(double htarg);
    bool isStartOfWake(int side, int stationIndex);
    FlowRegimeEnum applyFlowRegimeCandidate(FlowRegimeEnum candidate);
    FlowRegimeEnum currentFlowRegime() const;
    void updateSystemMatricesForStation(const Edge &edge, int side, int stationIndex, MixedModeStationContext &ctx);
    void initializeFirstIterationState(int side, int stationIndex, int previousTransition, MixedModeStationContext &ctx,
                                       double &ueref, double &hkref);
    void configureSimilarityRow(double ueref);
    void configureViscousRow(double hkref, double ueref, double senswt, bool resetSensitivity, bool averageSensitivity,
                             double &sens, double &sennew);
    bool applyMixedModeNewtonStep(int side, int stationIndex, double &ami, MixedModeStationContext &ctx);
    void checkTransitionIfNeeded(int side, int stationIndex, bool skipCheck, int laminarAdvance, double &ami);
    SetblOutputView setbl(SidePairRef<const BoundaryLayerSideProfiles> profiles, const FlowState &analysis_state,
                          const AeroCoefficients &aero_coeffs, double acrit, const Foil &foil,
                          const StagnationResult &stagnation, const Eigen::MatrixXd &dij, bool bl_initialized);
    void applySetblOutput(SetblOutputView &output);
    SkinFrictionCoefficients blmid(FlowRegimeEnum flowRegimeType);
    blData blprv(blData data, double xsi, double ami, double cti, double thi, double dsi, double dswaki,
                 double uei) const;
    bool blsys();
    SidePair<Eigen::VectorXd> ueset(const Eigen::MatrixXd &dij) const;
    EdgeVelocityDistribution computeNewUeDistribution(const XFoil &xfoil, const Matrix3x2dVector &vdel) const;
    ClContributions computeClFromEdgeVelocityDistribution(const XFoil &xfoil,
                                                          const EdgeVelocityDistribution &distribution) const;
    int resetSideState(int side, const Foil &foil, const StagnationResult &stagnation);
    void storeStationStateCommon(int side, int stationIndex, const MixedModeStationContext &ctx);
    void syncStationRegimeStates(int side, int stationIndex, FlowRegimeEnum stationRegime);
    FlowRegimeEnum determineRegimeForStation(int side, int stationIndex) const;
    double fallbackEdgeVelocity(int side, int stationIndex, EdgeVelocityFallbackMode edgeMode) const;
    void resetStationKinematicsAfterFailure(int side, int stationIndex, MixedModeStationContext &ctx,
                                            EdgeVelocityFallbackMode edgeMode);
    void recoverStationAfterFailure(int side, int stationIndex, MixedModeStationContext &ctx, double &ami,
                                    EdgeVelocityFallbackMode edgeMode, int laminarAdvance);
    BoundaryLayerDelta buildBoundaryLayerDelta(int side, const Eigen::VectorXd &unew_side,
                                               const Eigen::VectorXd &u_ac_side, double dac,
                                               const Matrix3x2dVector &vdel) const;
    BoundaryLayerMetrics evaluateSegmentRelaxation(int side, const BoundaryLayerDelta &delta, double dhi, double dlo,
                                                   double &relaxation) const;
    BoundaryLayerSideProfiles applyBoundaryLayerDelta(int side, const BoundaryLayerDelta &delta, double relaxation,
                                                      double hstinv, double gamm1) const;

    bool blkin(BoundaryLayerState &state);
    bool tesys(const BoundaryLayerSideProfiles &top_profiles, const BoundaryLayerSideProfiles &bottom_profiles,
               const Edge &edge);
    double xifset(const Foil &foil, const StagnationResult &stagnation, int is) const;

  private:
    static double adjustDisplacementForHkLimit(double displacementThickness, double momentumThickness, double msq,
                                               double hklim);

  public:
    Eigen::VectorXd wgap;
    SidePair<BoundaryLayerLattice> lattice;
    BlSystemCoeffs blc;
    BoundaryLayerState state;
    blDiff xt;
    int nsys             = 0;
    int stagnationIndex  = 0;
    double stagnationSst = 0.0;
    BoundaryLayerGeometry geometry;
};
