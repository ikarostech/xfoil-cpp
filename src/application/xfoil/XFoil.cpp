#include "application/xfoil/XFoil.h"

#include "application/xfoil/XFoilAnalysis.hpp"
#include "application/xfoil/XFoilWorkspace.hpp"

XFoil::XFoil()
    : workspace_(std::make_unique<XFoilWorkspace>()),
      analysis_(std::make_unique<XFoilAnalysis>(analysis_state_, result_,
                                                acrit_, *workspace_)) {}

XFoil::~XFoil() = default;

XFoilAnalysis &XFoil::analysis() { return *analysis_; }

const XFoilAnalysis &XFoil::analysis() const { return *analysis_; }

bool XFoil::initialize() { return analysis().initialize(); }

bool XFoil::initXFoilGeometry(int fn, const double *fx, const double *fy) {
  return analysis().initXFoilGeometry(fn, fx, fy);
}

bool XFoil::initXFoilAnalysis(double Re, double alpha, double Mach,
                              double NCrit, double XtrTop, double XtrBot,
                              ReynoldsType reType, MachType maType,
                              bool bViscous) {
  return analysis().initXFoilAnalysis(Re, alpha, Mach, NCrit, XtrTop, XtrBot,
                                      reType, maType, bViscous);
}

bool XFoil::specal() { return analysis().specal(); }

bool XFoil::speccl() { return analysis().speccl(); }

bool XFoil::viscal() { return analysis().viscal(); }

XFoil::ViscalEndResult XFoil::ViscalEnd() {
  const auto result = analysis().ViscalEnd();
  analysis().publishPressureCoefficients(result);
  analysis().publishTransitionLocations();
  return result;
}

bool XFoil::ViscousIter() { return analysis().ViscousIter(); }

bool XFoil::isBLInitialized() const { return analysis().isBLInitialized(); }

void XFoil::setBLInitialized(bool bInitialized) {
  analysis().setBLInitialized(bInitialized);
}

bool XFoil::hasConvergedSolution() const {
  return analysis().hasConvergedSolution();
}

void XFoil::invalidateConvergedSolution() {
  analysis().invalidateConvergedSolution();
}

void XFoil::invalidateWakeGeometry() { analysis().invalidateWakeGeometry(); }

void XFoil::invalidatePanelMap() { analysis().invalidatePanelMap(); }

double XFoil::QInf() const { return analysis_state_.qinf; }

void XFoil::setQInf(double v) { analysis_state_.qinf = v; }

double XFoil::alpha() const { return analysis_state_.alpha; }

void XFoil::setAlpha(double aoa) { analysis_state_.alpha = aoa; }

void XFoil::setControlByAlpha(bool enabled) {
  analysis_state_.controlByAlpha = enabled;
}

double XFoil::ClSpec() const { return analysis_state_.clspec; }

void XFoil::setClSpec(double cl) { analysis_state_.clspec = cl; }

bool XFoil::setMach() { return analysis().setMach(); }

double XFoil::Cl() const { return result_.aeroCoefficients.cl; }

double XFoil::Cd() const { return result_.aeroCoefficients.cd; }

double XFoil::Cm() const { return result_.aeroCoefficients.cm; }

double XFoil::getXcp() const { return result_.aeroCoefficients.xcp; }

const XFoil::VectorXd &XFoil::inviscidCp() const { return result_.inviscidCp; }

const XFoil::VectorXd &XFoil::viscousCp() const { return result_.viscousCp; }

double XFoil::transitionLocationTop() const {
  return result_.transitionLocations.top;
}

double XFoil::transitionLocationBottom() const {
  return result_.transitionLocations.bottom;
}

const XFoilResult &XFoil::result() const { return result_; }

FlowState &XFoil::analysisState() { return analysis_state_; }

const FlowState &XFoil::analysisState() const { return analysis_state_; }

double XFoil::VAccel() { return XFoilAnalysis::VAccel(); }

void XFoil::setVAccel(double accel) { XFoilAnalysis::setVAccel(accel); }
