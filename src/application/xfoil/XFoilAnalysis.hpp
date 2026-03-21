#pragma once

#include <stdexcept>

#include "Eigen/Core"
#include "Eigen/Dense"

#include "application/xfoil/XFoilResult.hpp"
#include "application/xfoil/XFoilSharedTypes.hpp"
#include "application/xfoil/XFoilWorkspace.hpp"
#include "model/flow_state.hpp"
#include "model/foil/FoilAerodynamicCache.hpp"
#include "numerics/side_pair.hpp"
#include "solver/xfoil/viscous_update.hpp"

struct SetblOutputView;
class InviscidSolver;

class XFoilAnalysis {
 public:
  using VectorXd = Eigen::VectorXd;
  using VectorXi = Eigen::VectorXi;
  using Vector2d = Eigen::Vector2d;
  using MatrixXd = Eigen::MatrixXd;
  using Matrix2Xd = Eigen::Matrix2Xd;
  using Matrix2d = Eigen::Matrix2d;
  template <typename T> using FullPivLU = Eigen::FullPivLU<T>;
  using Matrix3x2d = BoundaryLayerMatrix3x2d;
  using Matrix3x2dVector = BoundaryLayerMatrix3x2dVector;
  using EdgeVelocityDistribution = BoundaryLayerEdgeVelocityDistribution;
  using BoundaryLayerDelta = ::BoundaryLayerDelta;
  using BoundaryLayerMetrics = ::BoundaryLayerMetrics;
  using ReynoldsType = FlowState::ReynoldsType;
  using MachType = FlowState::MachType;
  using ClComputation = XFoilClComputation;
  using ViscalEndResult = XFoilViscalEndResult;
  using CompressibilityParams = XFoilCompressibilityParams;
  using PressureCoefficientResult = XFoilPressureCoefficientResult;
  using MixedModeStationContext = BoundaryLayerMixedModeStationContext;
  using UpdateResult = ViscousUpdateResult;
  using ClContributions = BoundaryLayerClContributions;
  enum class SideType { TOP = 1, BOTTOM = 2 };

  template <class T> class IterPair {
   public:
    T before;
    T after;

    T &get(int phase) {
      if (phase == 0) {
        return before;
      }
      if (phase == 1) {
        return after;
      }
      throw std::invalid_argument("invalid phase type");
    }
  };

  XFoilAnalysis(FlowState &analysis_state, XFoilResult &result, double &acrit,
                XFoilWorkspace &workspace);

  bool isValidFoilAngles(Matrix2Xd points);
  bool isValidFoilPointSize(Matrix2Xd points);
  bool initialize();
  bool initXFoilGeometry(int fn, const double *fx, const double *fy);
  bool initXFoilAnalysis(double Re, double alpha, double Mach, double NCrit,
                         double XtrTop, double XtrBot, ReynoldsType reType,
                         MachType maType, bool bViscous);
  bool specal();
  bool speccl();
  bool viscal();
  ViscalEndResult ViscalEnd();
  bool ViscousIter();
  bool abcopy(Matrix2Xd copyFrom);
  bool isBLInitialized() const;
  void setBLInitialized(bool bInitialized);
  bool setMach();
  double getActualMach(double cls, MachType mach_type);
  double getActualReynolds(double cls, ReynoldsType reynolds_type);
  void invalidateConvergedSolution();
  void invalidateWakeGeometry();
  void invalidatePanelMap();
  bool hasPanelMap() const;
  bool hasAirfoilInfluenceMatrix() const;
  bool hasWakeInfluenceMatrix() const;
  bool hasConvergedSolution() const;
  double cdcalc() const;
  double getXcp() const;
  ClComputation clcalc(Vector2d ref) const;
  void applyClComputation(const ClComputation &result);
  CompressibilityParams buildCompressibilityParams() const;
  PressureCoefficientResult
  computePressureCoefficient(double tangential_velocity,
                             double velocity_derivative,
                             const CompressibilityParams &params) const;
  Matrix2Xd gamqv() const;
  FoilAerodynamicCache ggcalc();
  bool qdcalc();
  VectorXd qvfue(const VectorXd &base_qvis) const;
  Matrix2Xd qwcalc(const Foil &foil, const Matrix2Xd &base_qinvu,
                   const Matrix2Xd &gamu, const Matrix2Xd &surface_vortex,
                   double alpha, double qinf) const;
  UpdateResult update(const Matrix3x2dVector &vdel) const;
  void ensureWakeTrajectoryAndInviscidVelocity();
  void ensurePanelMapAndBoundaryLayerGeometry();
  void assignCurrentInviscidEdgeVelocity();
  void ensureBoundaryLayerEdgeSeed();
  void restoreConvergedOperatingPoint(int point_count,
                                      int total_nodes_with_wake);
  void ensureSourceInfluenceMatrix();
  SetblOutputView initializeBoundaryLayerNewtonSystem();
  Matrix3x2dVector
  solveBoundaryLayerNewtonStep(const SetblOutputView &setbl_output) const;
  void applyBoundaryLayerIterationUpdate(const UpdateResult &update_result);
  void updateFreestreamForIteration(const UpdateResult &update_result);
  void refreshViscousFlowFields();
  void finalizeViscousIteration(const UpdateResult &update_result,
                                double eps1);
  void initializeDataStructures();
  void resetFlags();
  void resetVariables();
  void updateTrailingEdgeState();
  void publishPressureCoefficients(const ViscalEndResult &result);
  void publishTransitionLocations();

  const XFoilResult &result() const { return result_; }
  FlowState &analysisState() { return analysis_state_; }
  const FlowState &analysisState() const { return analysis_state_; }
  XFoilInternalState &internalState() { return state_; }
  const XFoilInternalState &internalState() const { return state_; }
  XFoilWorkspace &workspace() { return workspace_; }
  const XFoilWorkspace &workspace() const { return workspace_; }

  static double VAccel();
  static void setVAccel(double accel);
  inline static double vaccel_ = 0.01;

  friend class InviscidSolver;

 private:
  FlowState &analysis_state_;
  XFoilResult &result_;
  AeroCoefficients &aero_coeffs_;
  double &acrit_;
  XFoilWorkspace &workspace_;
  XFoilInternalState &state_;
  Foil &foil;
  BoundaryLayer &boundaryLayer;
  Eigen::VectorXd &cpi_;
  Eigen::VectorXd &cpv_;
  const Vector2d cmref = Vector2d{0.25, 0.0};
};
