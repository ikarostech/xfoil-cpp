#pragma once

#include "simulation/boundary_layer_state.hpp"
#include "domain/coefficient/bl_newton.hpp"
#include "domain/boundary_layer/boundary_layer_variables_solver.hpp"

class XFoil;
enum class FlowRegimeEnum;

class BoundaryLayerWorkflow {
 public:
  BoundaryLayerVariablesSolver boundaryLayerVariablesSolver;
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

  bool isStartOfWake(const XFoil& xfoil, int side, int stationIndex);
  void updateSystemMatricesForStation(XFoil& xfoil, int side,
                                      int stationIndex,
                                      MixedModeStationContext& ctx);
  void initializeFirstIterationState(XFoil& xfoil, int side,
                                     int stationIndex,
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
  blData blvar(blData data, FlowRegimeEnum flowRegimeType);

  bool iblpan(XFoil& xfoil);
  bool iblsys(XFoil& xfoil);
  bool stfind(XFoil& xfoil);
  bool stmove(XFoil& xfoil);
  bool tesys(XFoil& xfoil, double cte, double tte, double dte);
  bool trchek(XFoil& xfoil);

public:
  BoundaryLayerLattice lattice;
  BlSystemCoeffs blc;
  BoundaryLayerState state;
};
