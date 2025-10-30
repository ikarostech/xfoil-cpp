#pragma once

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

  static bool isStartOfWake(const XFoil& xfoil, int side, int stationIndex);
  static void updateSystemMatricesForStation(XFoil& xfoil, int side,
                                             int stationIndex,
                                             MixedModeStationContext& ctx);
  static void initializeFirstIterationState(XFoil& xfoil, int side,
                                            int stationIndex,
                                            int previousTransition,
                                            MixedModeStationContext& ctx,
                                            double& ueref, double& hkref,
                                            double& ami);
  static void configureSimilarityRow(XFoil& xfoil, double ueref);
  static void configureViscousRow(XFoil& xfoil, double hkref, double ueref,
                                  double senswt, bool resetSensitivity,
                                  bool averageSensitivity, double& sens,
                                  double& sennew);
  static bool applyMixedModeNewtonStep(XFoil& xfoil, int side, int stationIndex,
                                       double deps, double& ami,
                                       MixedModeStationContext& ctx);

  static bool iblpan(XFoil& xfoil);
  static bool iblsys(XFoil& xfoil);
  static bool stfind(XFoil& xfoil);
  static bool stmove(XFoil& xfoil);
  static bool tesys(XFoil& xfoil, double cte, double tte, double dte);
  static bool trchek(XFoil& xfoil);
};
