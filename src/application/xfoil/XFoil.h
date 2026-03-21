#pragma once

#include <memory>

#include "Eigen/Core"

#include "application/xfoil/XFoilSharedTypes.hpp"
#include "application/xfoil/XFoilResult.hpp"
#include "application/xfoil/XFoilWorkspace.hpp"
#include "model/flow_state.hpp"

class XFoilAnalysis;

class XFoil {
 public:
  using VectorXd = Eigen::VectorXd;
  using Vector2d = Eigen::Vector2d;
  using ReynoldsType = FlowState::ReynoldsType;
  using MachType = FlowState::MachType;
  using ViscalEndResult = XFoilViscalEndResult;

  XFoil();
  ~XFoil();
  XFoil(const XFoil &) = delete;
  XFoil &operator=(const XFoil &) = delete;
  XFoil(XFoil &&) = delete;
  XFoil &operator=(XFoil &&) = delete;

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

  bool isBLInitialized() const;
  void setBLInitialized(bool bInitialized);
  bool hasConvergedSolution() const;
  void invalidateConvergedSolution();
  void invalidateWakeGeometry();
  void invalidatePanelMap();

  double QInf() const;
  void setQInf(double v);
  double alpha() const;
  void setAlpha(double aoa);
  void setControlByAlpha(bool enabled);
  double ClSpec() const;
  void setClSpec(double cl);
  bool setMach();

  double Cl() const;
  double Cd() const;
  double Cm() const;
  double getXcp() const;

  const VectorXd &inviscidCp() const;
  const VectorXd &viscousCp() const;
  double transitionLocationTop() const;
  double transitionLocationBottom() const;
  const XFoilResult &result() const;

  FlowState &analysisState();
  const FlowState &analysisState() const;

  static double VAccel();
  static void setVAccel(double accel);

 private:
  XFoilAnalysis &analysis();
  const XFoilAnalysis &analysis() const;

  FlowState analysis_state_;
  XFoilResult result_;
  double acrit_ = 9.0;
  std::unique_ptr<XFoilWorkspace> workspace_;
  std::unique_ptr<XFoilAnalysis> analysis_;
};
