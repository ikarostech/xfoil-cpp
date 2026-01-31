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
#include "simulation/BoundaryLayerStore.hpp"
#include "domain/flow_regime.hpp"

class XFoil;
struct SetblOutputView;
class Foil;
class Edge;

class BoundaryLayerWorkflow {
 public:
  BoundaryLayerVariablesSolver boundaryLayerVariablesSolver;
  BoundaryLayerStore boundaryLayerStore;
  BlDiffSolver blDiffSolver;
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
  struct StagnationResult {
    int stagnationIndex = 0;
    double sst = 0.0;
    double sst_go = 0.0;
    double sst_gp = 0.0;
    bool found = true;
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
  struct MrchueStationContext {
    bool simi = false;
    bool wake = false;
    bool direct = true;
    double xsi = 0.0;
    double uei = 0.0;
    double thi = 0.0;
    double dsi = 0.0;
    double ami = 0.0;
    double cti = 0.0;
    double dswaki = 0.0;
    double cte = 0.0;
    double dte = 0.0;
    double tte = 0.0;
    double dmax = 0.0;
    double hmax = 0.0;
    double htarg = 0.0;
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
  blData blvar(blData data, FlowRegimeEnum flowRegimeType);
  SkinFrictionCoefficients blmid(FlowRegimeEnum flowRegimeType);
  blData blprv(blData data, double xsi, double ami, double cti,
               double thi, double dsi, double dswaki, double uei) const;
  bool blsys();
  bool trdif();
  bool mrchdu(XFoil& xfoil);
  bool mrchdu(BoundaryLayerState& state, XFoil& xfoil);
  int resetSideState(int side, XFoil& xfoil);
  bool marchBoundaryLayerSide(BoundaryLayerState& state, int side,
                              double deps, double senswt, double& sens,
                              double& sennew, double& ami, XFoil& xfoil);
  bool processBoundaryLayerStation(BoundaryLayerState& state, int side,
                                   int stationIndex, int previousTransition,
                                   double deps, double senswt, double& sens,
                                   double& sennew, double& ami, XFoil& xfoil);
  bool mrchue(XFoil& xfoil);
  bool mrchue(BoundaryLayerState& state, XFoil& xfoil);
  bool marchMrchueSide(BoundaryLayerState& state, int side,
                       XFoil& xfoil, std::stringstream& ss);
  void initializeMrchueSide(int side, double& thi, double& dsi,
                            double& ami, double& cti);
  void prepareMrchueStationContext(int side, int stationIndex,
                                   MrchueStationContext& ctx, double thi,
                                   double dsi, double ami, double cti);
  bool performMrchueNewtonLoop(int side, int stationIndex,
                               MrchueStationContext& ctx, XFoil& xfoil,
                               std::stringstream& ss);
  void handleMrchueStationFailure(int side, int stationIndex,
                                  MrchueStationContext& ctx, XFoil& xfoil,
                                  std::stringstream& ss);
  void storeMrchueStationState(int side, int stationIndex,
                               const MrchueStationContext& ctx);
  void storeStationStateCommon(int side, int stationIndex, double ami,
                               double cti, double thi, double dsi, double uei,
                               double xsi, double dswaki);
  template <typename StationContext>
  void resetStationKinematicsAfterFailure(int side, int stationIndex,
                                          StationContext& ctx,
                                          EdgeVelocityFallbackMode edgeMode);
  double fallbackEdgeVelocity(int side, int stationIndex,
                              EdgeVelocityFallbackMode edgeMode) const;
  EdgeVelocityDistribution computeNewUeDistribution(const XFoil& xfoil) const;
  QtanResult computeQtan(const EdgeVelocityDistribution& distribution,
                         int point_count) const;
  ClContributions computeClFromEdgeVelocityDistribution(
      const XFoil& xfoil, const EdgeVelocityDistribution& distribution) const;
  BoundaryLayerDelta buildBoundaryLayerDelta(
      int side, const Eigen::VectorXd& unew_side,
      const Eigen::VectorXd& u_ac_side, double dac,
      const XFoil& xfoil) const;
  BoundaryLayerMetrics evaluateSegmentRelaxation(
      int side, const BoundaryLayerDelta& delta, double dhi, double dlo,
      double& relaxation) const;
  BoundaryLayerSideProfiles applyBoundaryLayerDelta(
      int side, const BoundaryLayerDelta& delta, double relaxation,
      double hstinv, double gamm1) const;
  SidePair<Eigen::VectorXd> ueset(const Eigen::MatrixXd& dij) const;
  void syncStationRegimeStates(int side, int stationIndex, bool wake);
  FlowRegimeEnum determineRegimeForStation(int side, int stationIndex,
                                           bool similarity,
                                           bool wake) const;
  MixedModeStationContext prepareMixedModeStation(int side, int stationIndex,
                                                  int previousTransition,
                                                  double& ami);

  bool iblpan(int point_count, int wake_point_count,
              std::string* error_message);
  bool iblsys(XFoil& xfoil);
  StagnationResult stfind(const Eigen::Matrix2Xd& surface_vortex,
                          const Eigen::VectorXd& spline_length) const;
  bool stmove(XFoil& xfoil);
  bool xicalc(const Foil& foil);
  static SidePair<Eigen::VectorXd> computeArcLengthCoordinates(
      const Foil& foil, double stagnationSst,
      const SidePair<BoundaryLayerLattice>& lattice);
  static Eigen::VectorXd computeWakeGap(
      const Foil& foil, const BoundaryLayerLattice& bottom,
      const Eigen::VectorXd& bottomArcLengths);
  SidePair<Eigen::Matrix2Xd> uicalc(const Eigen::Matrix2Xd& qinv_matrix) const;
  bool blkin(BoundaryLayerState& state);
  bool tesys(const BoundaryLayerSideProfiles& top_profiles,
             const BoundaryLayerSideProfiles& bottom_profiles,
             const Edge& edge);
  double calcHtarg(int ibl, int is, bool wake);
  bool trchek(XFoil& xfoil);
  double xifset(const Foil& foil, const StagnationResult& stagnation,
                int is) const;

private:
  static double adjustDisplacementForHkLimit(double displacementThickness,
                                             double momentumThickness,
                                             double msq, double hklim);
  double computeTransitionLocation(double weightingFactor) const;
  void copyStationState(int side, int destination, int source);

public:
  Eigen::VectorXd wgap;
  SidePair<BoundaryLayerLattice> lattice;
  BlSystemCoeffs blc;
  BoundaryLayerState state;
  blDiff xt;
  int stagnationIndex = 0;
  double stagnationSst = 0.0;
};

template <typename StationContext>
inline void BoundaryLayerWorkflow::resetStationKinematicsAfterFailure(
    int side, int stationIndex, StationContext& ctx,
    EdgeVelocityFallbackMode edgeMode) {
  if (ctx.dmax <= 0.1 || stationIndex < 2) {
    return;
  }

  if (stationIndex <= lattice.get(side).trailingEdgeIndex) {
    const double ratio =
        lattice.get(side).arcLengthCoordinates[stationIndex] /
        lattice.get(side).arcLengthCoordinates[stationIndex - 1];
    const double scale = std::sqrt(ratio);
    ctx.thi = lattice.get(side).profiles.momentumThickness[stationIndex - 1] * scale;
    ctx.dsi = lattice.get(side).profiles.displacementThickness[stationIndex - 1] * scale;
  } else {
    if (stationIndex == lattice.get(side).trailingEdgeIndex + 1) {
      ctx.cti = ctx.cte;
      ctx.thi = ctx.tte;
      ctx.dsi = ctx.dte;
    } else {
      ctx.thi = lattice.get(side).profiles.momentumThickness[stationIndex - 1];
      const double ratlen =
          (lattice.get(side).arcLengthCoordinates[stationIndex] -
           lattice.get(side).arcLengthCoordinates[stationIndex - 1]) /
          (10.0 * lattice.get(side).profiles.displacementThickness[stationIndex - 1]);
      ctx.dsi =
          (lattice.get(side).profiles.displacementThickness[stationIndex - 1] + ctx.thi * ratlen) /
          (1.0 + ratlen);
    }
  }

  ctx.uei = fallbackEdgeVelocity(side, stationIndex, edgeMode);

  if (stationIndex == lattice.get(side).profiles.transitionIndex) {
    ctx.cti = 0.05;
  }
  if (stationIndex > lattice.get(side).profiles.transitionIndex) {
    ctx.cti = lattice.get(side).profiles.skinFrictionCoeff[stationIndex - 1];
  }
}
