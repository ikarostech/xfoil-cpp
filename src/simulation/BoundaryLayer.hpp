#pragma once

#include "simulation/boundary_layer_state.hpp"

class XFoil;

class BoundaryLayerWorkflow {
 public:
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
  void configureSimilarityRow(XFoil& xfoil, double ueref);
  void configureViscousRow(XFoil& xfoil, double hkref, double ueref,
                           double senswt, bool resetSensitivity,
                           bool averageSensitivity, double& sens,
                           double& sennew);
  bool applyMixedModeNewtonStep(XFoil& xfoil, int side, int stationIndex,
                                double deps, double& ami,
                                MixedModeStationContext& ctx);

  bool iblpan(XFoil& xfoil);
  bool iblsys(XFoil& xfoil);
  bool stfind(XFoil& xfoil);
  bool stmove(XFoil& xfoil);
  bool tesys(XFoil& xfoil, double cte, double tte, double dte);
  bool trchek(XFoil& xfoil);

 public:
  BoundaryLayerLattice lattice;
};
