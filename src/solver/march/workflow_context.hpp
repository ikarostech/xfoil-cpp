#pragma once

#include "solver/boundary_layer/march_access.hpp"
#include "solver/march/context.hpp"

class WorkflowMarchContext final : public MarchContext {
 public:
  explicit WorkflowMarchContext(BoundaryLayerMarchAccess access)
      : access_(access) {}

  double readCurrentShapeFactor() const override {
    return access_.readCurrentShapeFactor();
  }
  void refreshCurrentStationState(double xsi, double ami, double cti,
                                  double thi, double dsi, double dswaki,
                                  double uei) override {
    access_.refreshCurrentStationState(xsi, ami, cti, thi, dsi, dswaki, uei);
  }

  int readSideStationCount(int side) const override {
    return access_.readSideStationCount(side);
  }
  StationReadModel readStationModel(int side, int stationIndex) const override {
    return access_.readStationModel(side, stationIndex);
  }
  TrailingEdgeReadModel readTrailingEdgeModel() const override {
    return access_.readTrailingEdgeModel();
  }

  void emitMarchInfoLog(std::string_view message) const override {
    access_.emitMarchInfoLog(message);
  }

  double readNewtonRhs(int row) const override { return access_.readNewtonRhs(row); }
  void solveMrchueDirectNewtonSystem() override {
    access_.solveMrchueDirectNewtonSystem();
  }
  void solveMrchueInverseNewtonSystem(double htarg) override {
    access_.solveMrchueInverseNewtonSystem(htarg);
  }
  void runTransitionCheckForMrchue(int side, int stationIndex, double &ami,
                                   double &cti) override {
    access_.runTransitionCheckForMrchue(side, stationIndex, ami, cti);
  }
  bool solveTeSystemForCurrentProfiles(const Edge &edge) override {
    return access_.solveTeSystemForCurrentProfiles(edge);
  }

  double readBlCompressibilityHstinv() const override {
    return access_.blCompressibility().hstinv();
  }
  double readBlCompressibilityGm1bl() const override {
    return access_.blCompressibility().gm1bl();
  }
  double readBlReynoldsReybl() const override {
    return access_.blReynolds().reybl;
  }

  bool isStartOfWake(int side, int stationIndex) override {
    return access_.isStartOfWake(side, stationIndex);
  }
  FlowRegimeEnum applyFlowRegimeCandidate(FlowRegimeEnum candidate) override {
    return access_.applyFlowRegimeCandidate(candidate);
  }
  FlowRegimeEnum currentFlowRegime() const override {
    return access_.currentFlowRegime();
  }
  FlowRegimeEnum determineRegimeForStation(int side,
                                           int stationIndex) const override {
    return access_.determineRegimeForStation(side, stationIndex);
  }
  int resetSideState(int side, const Foil &foil,
                     const StagnationResult &stagnation) override {
    return access_.resetSideState(side, foil, stagnation);
  }

  void updateSystemMatricesForStation(const Edge &edge, int side,
                                      int stationIndex,
                                      MixedModeStationContext &ctx) override {
    access_.updateSystemMatricesForStation(edge, side, stationIndex, ctx);
  }
  void initializeFirstIterationState(int side, int stationIndex,
                                     int previousTransition,
                                     MixedModeStationContext &ctx,
                                     double &ueref, double &hkref) override {
    access_.initializeFirstIterationState(side, stationIndex,
                                          previousTransition, ctx, ueref, hkref);
  }
  void configureSimilarityRow(double ueref) override {
    access_.configureSimilarityRow(ueref);
  }
  void configureViscousRow(double hkref, double ueref, double senswt,
                           bool resetSensitivity, bool averageSensitivity,
                           double &sens, double &sennew) override {
    access_.configureViscousRow(hkref, ueref, senswt, resetSensitivity,
                                averageSensitivity, sens, sennew);
  }
  bool applyMixedModeNewtonStep(int side, int stationIndex, double &ami,
                                MixedModeStationContext &ctx) override {
    return access_.applyMixedModeNewtonStep(side, stationIndex, ami, ctx);
  }
  void checkTransitionIfNeeded(int side, int stationIndex, bool skipCheck,
                               int laminarAdvance, double &ami) override {
    access_.checkTransitionIfNeeded(side, stationIndex, skipCheck,
                                    laminarAdvance, ami);
  }
  void storeStationStateCommon(int side, int stationIndex,
                               const MixedModeStationContext &ctx) override {
    access_.storeStationStateCommon(side, stationIndex, ctx);
  }
  void recoverStationAfterFailure(
      int side, int stationIndex, MixedModeStationContext &ctx, double &ami,
      EdgeVelocityFallbackMode edgeMode, int laminarAdvance) override {
    access_.recoverStationAfterFailure(side, stationIndex, ctx, ami, edgeMode,
                                       laminarAdvance);
  }

  bool blsys() override { return access_.blsys(); }
  double calcHtarg(int ibl, int is, bool wake) override {
    return access_.calcHtarg(ibl, is, wake);
  }

 private:
  BoundaryLayerMarchAccess access_;
};
