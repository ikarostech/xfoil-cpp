#pragma once

#include <array>
#include <cmath>
#include <sstream>
#include <string>

#include "simulation/boundary_layer_state.hpp"
#include "domain/coefficient/bl_newton.hpp"
#include "domain/boundary_layer/boundary_layer_variables_solver.hpp"
#include "domain/boundary_layer/bl_compressibility_params.hpp"
#include "domain/boundary_layer/bl_reynolds_params.hpp"
#include "domain/boundary_layer/bl_transition_params.hpp"
#include "simulation/skin_friction_coefficients.hpp"
#include "domain/boundary_layer/boundary_layer_diff_solver.hpp"
#include "simulation/BoundaryLayer_transition.hpp"
#include "simulation/boundary_layer_geometry.hpp"
#include "domain/flow_regime.hpp"

class XFoil;
struct SetblOutputView;
class Foil;
class Edge;
class BoundaryLayerWorkflow {
 public:
  BoundaryLayerWorkflow();

  BoundaryLayerVariablesSolver boundaryLayerVariablesSolver;
  BlDiffSolver blDiffSolver;
  BoundaryLayerTransitionSolver transitionSolver;
  FlowRegimeEnum flowRegime = FlowRegimeEnum::Laminar;
  BlCompressibilityParams blCompressibility{};
  BlReynoldsParams blReynolds{};
  BlTransitionParams blTransition{};

  // Sutherland's const./T0 (assumes stagnation conditions are at STP).
  static constexpr double kHvrat = 0.35;

  enum class EdgeVelocityFallbackMode {
    UsePreviousStation,
    AverageNeighbors
  };
  struct EdgeVelocityDistribution {
    SidePair<Eigen::VectorXd> unew;
    SidePair<Eigen::VectorXd> u_ac;
  };
  struct LeTeSensitivities {
    SidePair<Eigen::VectorXd> ule_m;
    SidePair<Eigen::VectorXd> ute_m;
  };
  struct EdgeVelocitySensitivityResult {
    SidePair<Eigen::VectorXd> edgeVelocity;
    SidePair<Eigen::VectorXd> outputEdgeVelocity;
    SidePair<int> jvte;
    SidePair<double> dule;
    SidePair<Eigen::VectorXd> ule_m;
    SidePair<Eigen::VectorXd> ute_m;
    SidePair<double> ule_a;
  };
  struct QtanResult {
    Eigen::VectorXd qnew;
    Eigen::VectorXd q_ac;
  };
  struct ClContributions {
    double cl = 0.0;
    double cl_a = 0.0;
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
    double maxChange = 0.0;
  };
  struct MixedModeStationContext {
    bool simi = false;
    bool wake = false;
    double xsi = 0.0;
    double uei = 0.0;
    double thi = 0.0;
    double dsi = 0.0;
    double cti = 0.0;
    double ami = 0.0;
    double dswaki = 0.0;
    double cte = 0.0;
    double dte = 0.0;
    double tte = 0.0;
    double dmax = 0.0;
  };
  struct BlInitializationPlan {
    bool needsInitialization = false;
    std::string message;
  };
  struct StationPrimaryVars {
    double xsi = 0.0;
    double uei = 0.0;
    double thi = 0.0;
    double mdi = 0.0;
    double dsi = 0.0;
    double dswaki = 0.0;
    double ami = 0.0;
    double cti = 0.0;
  };
  struct TeWakeCoefficients {
    double tte = 0.0;
    double cte = 0.0;
    double dte = 0.0;
    double tte_tte1 = 0.0;
    double tte_tte2 = 0.0;
    double dte_mte1 = 0.0;
    double dte_ute1 = 0.0;
    double dte_mte2 = 0.0;
    double dte_ute2 = 0.0;
    double cte_cte1 = 0.0;
    double cte_cte2 = 0.0;
    double cte_tte1 = 0.0;
    double cte_tte2 = 0.0;
  };
  struct SimilarityStationCoefficients {
    Eigen::VectorXd u_m1;
    Eigen::VectorXd d_m1;
  };
  struct SideSweepInitResult {
    double u_a1 = 0.0;
    double d_a1 = 0.0;
    double due1 = 0.0;
    double dds1 = 0.0;
    double xiforc = 0.0;
  };
  struct StationUpdateResult {
    BoundaryLayerState state;
    Eigen::VectorXd u_m2;
    Eigen::VectorXd d_m2;
    double u_a2 = 0.0;
    double d_a2 = 0.0;
    double due2 = 0.0;
    double dds2 = 0.0;
  };
  struct TransitionLogResult {
    std::string message;
  };
  struct TeWakeUpdateResult {
    bool isStartOfWake = false;
    Eigen::VectorXd d_m1;
    double due1 = 0.0;
    double dds1 = 0.0;
    TeWakeCoefficients coeffs;
  };
  struct TeWakeJacobianAdjustments {
    double vz[3][2] = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
    Eigen::Matrix<double, 3, 2> vb = Eigen::Matrix<double, 3, 2>::Zero();
  };
  struct StationArraysAdvanceResult {
    Eigen::VectorXd u_m1;
    Eigen::VectorXd d_m1;
    double u_a1 = 0.0;
    double d_a1 = 0.0;
    double due1 = 0.0;
    double dds1 = 0.0;
  };

  bool isStartOfWake(int side, int stationIndex);
  void updateSystemMatricesForStation(XFoil& xfoil, int side,
                                      int stationIndex,
                                      MixedModeStationContext& ctx);
  void initializeFirstIterationState(int side, int stationIndex,
                                     int previousTransition,
                                     MixedModeStationContext& ctx,
                                     double& ueref, double& hkref,
                                     double& ami);
  void configureSimilarityRow(double ueref);
  void configureViscousRow(double hkref, double ueref,
                           double senswt, bool resetSensitivity,
                           bool averageSensitivity, double& sens,
                           double& sennew);
  bool applyMixedModeNewtonStep(XFoil& xfoil, int side, int stationIndex,
                                double deps, double& ami,
                                MixedModeStationContext& ctx);
  void checkTransitionIfNeeded(XFoil& xfoil, int side, int stationIndex,
                               bool skipCheck, int laminarAdvance,
                               double& ami);
  BlInitializationPlan computeBlInitializationPlan(bool lblini) const;
  SimilarityStationCoefficients resetSimilarityStationCoefficients(
      const Eigen::VectorXd& u_m1, const Eigen::VectorXd& d_m1) const;
  SideSweepInitResult initializeSideSweepState(
      const Foil& foil, const StagnationResult& stagnation, int is) const;
  StationPrimaryVars loadStationPrimaryVars(int is, int ibl,
                                            bool stationIsWake,
                                            const SetblOutputView& output,
                                            double ami,
                                            double cti) const;
  StationUpdateResult updateStationMatricesAndState(
      int is, int ibl, int iv, const StationPrimaryVars& vars,
      const SidePair<Eigen::VectorXd>& usav, const SetblOutputView& output,
      const BoundaryLayerState& base_state, int system_size,
      const Eigen::MatrixXd& dij);
  TransitionLogResult buildTransitionLog(
      bool stationIsTransitionCandidate, FlowRegimeEnum flowRegime) const;
  TeWakeUpdateResult computeTeWakeCoefficients(
      int is, int ibl, const SidePair<Eigen::VectorXd>& usav,
      const SidePair<Eigen::VectorXd>& ute_m, const SidePair<int>& jvte,
      const Eigen::VectorXd& d_m1_template, const SetblOutputView& output,
      const Edge& edge) const;
  TeWakeJacobianAdjustments computeTeWakeJacobianAdjustments(
      const TeWakeCoefficients& coeffs) const;
  EdgeVelocitySensitivityResult prepareEdgeVelocityAndSensitivities(
      SidePairRef<const BoundaryLayerSideProfiles> profiles,
      const Eigen::MatrixXd& dij, int nsys) const;
  LeTeSensitivities computeLeTeSensitivities(int ile1, int ile2, int ite1,
                                             int ite2, int nsys,
                                             const Eigen::MatrixXd& dij) const;
  StationArraysAdvanceResult advanceStationArrays(
      const Eigen::VectorXd& u_m2, const Eigen::VectorXd& d_m2, double u_a2,
      double d_a2, double due2, double dds2) const;
  void assembleBlJacobianForStation(
      int is, int iv, int nsys, const SidePairRef<const Eigen::VectorXd>& d_m,
      const SidePairRef<const Eigen::VectorXd>& u_m,
      const SidePairRef<const double>& xi_ule,
      const SidePairRef<const Eigen::VectorXd>& ule_m,
      const SidePairRef<const double>& ule_a,
      const SidePairRef<const double>& u_a,
      const SidePairRef<const double>& d_a,
      const SidePairRef<const double>& due,
      const SidePairRef<const double>& dds,
      const SidePairRef<const double>& dule, bool controlByAlpha,
      double re_clmr, double msq_clmr, SetblOutputView& output);
  SkinFrictionCoefficients blmid(FlowRegimeEnum flowRegimeType);
  blData blprv(blData data, double xsi, double ami, double cti,
               double thi, double dsi, double dswaki, double uei) const;
  bool blsys();
  SidePair<Eigen::VectorXd> ueset(const Eigen::MatrixXd& dij) const;

  bool iblpan(int point_count, int wake_point_count,
              std::string* error_message);
  bool iblsys(XFoil& xfoil);
  StagnationResult stfind(const Eigen::Matrix2Xd& surface_vortex,
                          const Eigen::VectorXd& spline_length) const;
  bool stmove(XFoil& xfoil);
  bool xicalc(const Foil& foil);
  SidePair<Eigen::Matrix2Xd> uicalc(const Eigen::Matrix2Xd& qinv_matrix) const;
  bool blkin(BoundaryLayerState& state);
  bool tesys(const BoundaryLayerSideProfiles& top_profiles,
             const BoundaryLayerSideProfiles& bottom_profiles,
             const Edge& edge);
  double calcHtarg(int ibl, int is, bool wake);
  double xifset(const Foil& foil, const StagnationResult& stagnation,
                int is) const;

private:
  static double adjustDisplacementForHkLimit(double displacementThickness,
                                             double momentumThickness,
                                             double msq, double hklim);

public:
  Eigen::VectorXd wgap;
  SidePair<BoundaryLayerLattice> lattice;
  BlSystemCoeffs blc;
  BoundaryLayerState state;
  blDiff xt;
  int stagnationIndex = 0;
  double stagnationSst = 0.0;
  BoundaryLayerGeometry geometry;
};
