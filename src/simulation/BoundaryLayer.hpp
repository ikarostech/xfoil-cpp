#pragma once

#include <cmath>
#include <sstream>

#include "simulation/boundary_layer_state.hpp"
#include "domain/coefficient/bl_newton.hpp"
#include "domain/boundary_layer/boundary_layer_variables_solver.hpp"
#include "simulation/skin_friction_coefficients.hpp"

class XFoil;
enum class FlowRegimeEnum;

class BoundaryLayerWorkflow {
 public:
  BoundaryLayerVariablesSolver boundaryLayerVariablesSolver;
  enum class EdgeVelocityFallbackMode {
    UsePreviousStation,
    AverageNeighbors
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
  SkinFrictionCoefficients blmid(XFoil& xfoil,
                                 FlowRegimeEnum flowRegimeType);
  blData blprv(XFoil& xfoil, blData data, double xsi, double ami, double cti,
               double thi, double dsi, double dswaki, double uei) const;
  bool blsys(XFoil& xfoil);
  bool mrchdu(XFoil& xfoil);
  bool mrchdu(BoundaryLayerState& state, BoundaryLayerLattice& lattice,
              XFoil& xfoil);
  int resetSideState(int side, XFoil& xfoil);
  bool marchBoundaryLayerSide(BoundaryLayerState& state, int side,
                              double deps, double senswt, double& sens,
                              double& sennew, double& ami, XFoil& xfoil);
  bool processBoundaryLayerStation(BoundaryLayerState& state, int side,
                                   int stationIndex, int previousTransition,
                                   double deps, double senswt, double& sens,
                                   double& sennew, double& ami, XFoil& xfoil);
  bool mrchue(XFoil& xfoil);
  bool mrchue(BoundaryLayerState& state, BoundaryLayerLattice& lattice,
              XFoil& xfoil);
  bool marchMrchueSide(BoundaryLayerState& state, int side,
                       XFoil& xfoil, std::stringstream& ss);
  void initializeMrchueSide(int side, double& thi, double& dsi,
                            double& ami, double& cti, XFoil& xfoil);
  void prepareMrchueStationContext(int side, int stationIndex,
                                   MrchueStationContext& ctx, double thi,
                                   double dsi, double ami, double cti,
                                   XFoil& xfoil);
  bool performMrchueNewtonLoop(int side, int stationIndex,
                               MrchueStationContext& ctx, XFoil& xfoil,
                               std::stringstream& ss);
  void handleMrchueStationFailure(int side, int stationIndex,
                                  MrchueStationContext& ctx, XFoil& xfoil,
                                  std::stringstream& ss);
  void storeMrchueStationState(int side, int stationIndex,
                               const MrchueStationContext& ctx, XFoil& xfoil);
  void storeStationStateCommon(int side, int stationIndex, double ami,
                               double cti, double thi, double dsi, double uei,
                               double xsi, double dswaki, XFoil& xfoil);
  template <typename StationContext>
  void resetStationKinematicsAfterFailure(int side, int stationIndex,
                                          StationContext& ctx,
                                          EdgeVelocityFallbackMode edgeMode);
  double fallbackEdgeVelocity(int side, int stationIndex,
                              EdgeVelocityFallbackMode edgeMode) const;
  void syncStationRegimeStates(int side, int stationIndex, bool wake,
                               XFoil& xfoil);

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

template <typename StationContext>
inline void BoundaryLayerWorkflow::resetStationKinematicsAfterFailure(
    int side, int stationIndex, StationContext& ctx,
    EdgeVelocityFallbackMode edgeMode) {
  if (ctx.dmax <= 0.1 || stationIndex < 2) {
    return;
  }

  if (stationIndex <= lattice.trailingEdgeIndex.get(side)) {
    const double ratio =
        lattice.xssi.get(side)[stationIndex] /
        lattice.xssi.get(side)[stationIndex - 1];
    const double scale = std::sqrt(ratio);
    ctx.thi = lattice.thet.get(side)[stationIndex - 1] * scale;
    ctx.dsi = lattice.dstr.get(side)[stationIndex - 1] * scale;
  } else {
    if (stationIndex == lattice.trailingEdgeIndex.get(side) + 1) {
      ctx.cti = ctx.cte;
      ctx.thi = ctx.tte;
      ctx.dsi = ctx.dte;
    } else {
      ctx.thi = lattice.thet.get(side)[stationIndex - 1];
      const double ratlen =
          (lattice.xssi.get(side)[stationIndex] -
           lattice.xssi.get(side)[stationIndex - 1]) /
          (10.0 * lattice.dstr.get(side)[stationIndex - 1]);
      ctx.dsi =
          (lattice.dstr.get(side)[stationIndex - 1] + ctx.thi * ratlen) /
          (1.0 + ratlen);
    }
  }

  ctx.uei = fallbackEdgeVelocity(side, stationIndex, edgeMode);

  if (stationIndex == lattice.transitionIndex.get(side)) {
    ctx.cti = 0.05;
  }
  if (stationIndex > lattice.transitionIndex.get(side)) {
    ctx.cti = lattice.ctau.get(side)[stationIndex - 1];
  }
}
