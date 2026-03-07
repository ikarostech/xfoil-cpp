#pragma once

#include "simulation/BoundaryLayer.hpp"
#include "simulation/march/context.hpp"

class WorkflowMarchContext final : public MarchContext
{
public:
  explicit WorkflowMarchContext(BoundaryLayerWorkflow &workflow)
      : workflow_(workflow) {}

  BoundaryLayerState &mutableState() override { return workflow_.state; }

  int readSideStationCount(int side) const override {
    return workflow_.readSideStationCount(side);
  }
  StationReadModel readStationModel(int side, int stationIndex) const override {
    return workflow_.readStationModel(side, stationIndex);
  }
  TrailingEdgeReadModel readTrailingEdgeModel() const override {
    return workflow_.readTrailingEdgeModel();
  }

  void emitMarchInfoLog(std::string_view message) const override {
    workflow_.emitMarchInfoLog(message);
  }

  double readNewtonRhs(int row) const override {
    return workflow_.readNewtonRhs(row);
  }
  void solveMrchueDirectNewtonSystem() override {
    workflow_.solveMrchueDirectNewtonSystem();
  }
  void solveMrchueInverseNewtonSystem(double htarg) override {
    workflow_.solveMrchueInverseNewtonSystem(htarg);
  }
  void runTransitionCheckForMrchue(int side, int stationIndex, double &ami,
                                   double &cti) override {
    workflow_.runTransitionCheckForMrchue(side, stationIndex, ami, cti);
  }
  bool solveTeSystemForCurrentProfiles(const Edge &edge) override {
    return workflow_.solveTeSystemForCurrentProfiles(edge);
  }

  double readBlCompressibilityHstinv() const override {
    return workflow_.blCompressibility.hstinv;
  }
  double readBlCompressibilityGm1bl() const override {
    return workflow_.blCompressibility.gm1bl;
  }
  double readBlReynoldsReybl() const override {
    return workflow_.blReynolds.reybl;
  }

  bool isStartOfWake(int side, int stationIndex) override {
    return workflow_.isStartOfWake(side, stationIndex);
  }
  FlowRegimeEnum
  applyFlowRegimeCandidate(FlowRegimeEnum candidate) override {
    return workflow_.applyFlowRegimeCandidate(candidate);
  }
  FlowRegimeEnum currentFlowRegime() const override {
    return workflow_.currentFlowRegime();
  }
  FlowRegimeEnum determineRegimeForStation(int side,
                                           int stationIndex) const override {
    return workflow_.determineRegimeForStation(side, stationIndex);
  }
  int resetSideState(int side, const Foil &foil,
                     const StagnationResult &stagnation) override {
    return workflow_.resetSideState(side, foil, stagnation);
  }

  void updateSystemMatricesForStation(const Edge &edge, int side,
                                      int stationIndex,
                                      MixedModeStationContext &ctx) override {
    workflow_.updateSystemMatricesForStation(edge, side, stationIndex, ctx);
  }
  void initializeFirstIterationState(int side, int stationIndex,
                                     int previousTransition,
                                     MixedModeStationContext &ctx,
                                     double &ueref, double &hkref) override {
    workflow_.initializeFirstIterationState(side, stationIndex,
                                            previousTransition, ctx, ueref,
                                            hkref);
  }
  void configureSimilarityRow(double ueref) override {
    workflow_.configureSimilarityRow(ueref);
  }
  void configureViscousRow(double hkref, double ueref, double senswt,
                           bool resetSensitivity, bool averageSensitivity,
                           double &sens, double &sennew) override {
    workflow_.configureViscousRow(hkref, ueref, senswt, resetSensitivity,
                                  averageSensitivity, sens, sennew);
  }
  bool applyMixedModeNewtonStep(int side, int stationIndex, double &ami,
                                MixedModeStationContext &ctx) override {
    return workflow_.applyMixedModeNewtonStep(side, stationIndex, ami, ctx);
  }
  void checkTransitionIfNeeded(int side, int stationIndex, bool skipCheck,
                               int laminarAdvance, double &ami) override {
    workflow_.checkTransitionIfNeeded(side, stationIndex, skipCheck,
                                      laminarAdvance, ami);
  }
  void storeStationStateCommon(int side, int stationIndex,
                               const MixedModeStationContext &ctx) override {
    workflow_.storeStationStateCommon(side, stationIndex, ctx);
  }
  void recoverStationAfterFailure(
      int side, int stationIndex, MixedModeStationContext &ctx, double &ami,
      EdgeVelocityFallbackMode edgeMode, int laminarAdvance) override {
    workflow_.recoverStationAfterFailure(side, stationIndex, ctx, ami, edgeMode,
                                         laminarAdvance);
  }

  blData blprv(blData data, double xsi, double ami, double cti, double thi,
               double dsi, double dswaki, double uei) const override {
    return workflow_.blprv(data, xsi, ami, cti, thi, dsi, dswaki, uei);
  }
  bool blkin(BoundaryLayerState &state) override { return workflow_.blkin(state); }
  bool blsys() override { return workflow_.blsys(); }
  double calcHtarg(int ibl, int is, bool wake) override {
    return workflow_.calcHtarg(ibl, is, wake);
  }

private:
  BoundaryLayerWorkflow &workflow_;
};
