#pragma once

#include <cmath>

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
    workflow_.transitionSolver.trchek();
    ami = workflow_.state.station2.param.amplz;
    if (workflow_.flowRegime == FlowRegimeEnum::Transition) {
      workflow_.lattice.get(side).profiles.transitionIndex = stationIndex;
      if (cti <= 0.0) {
        cti = 0.03;
        workflow_.state.station2.param.sz = cti;
      }
    } else {
      workflow_.lattice.get(side).profiles.transitionIndex = stationIndex + 2;
    }
  }
  bool solveTeSystemForCurrentProfiles(const Edge &edge) override {
    return workflow_.tesys(workflow_.lattice.top.profiles,
                           workflow_.lattice.bottom.profiles, edge);
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
    if (ibl < workflow_.lattice.get(is).profiles.transitionIndex) {
      return workflow_.state.station1.hkz.scalar +
             0.03 * (workflow_.state.station2.param.xz -
                     workflow_.state.station1.param.xz) /
                 workflow_.state.station1.param.tz;
    }

    if (ibl == workflow_.lattice.get(is).profiles.transitionIndex) {
      return workflow_.state.station1.hkz.scalar +
             (0.03 * (workflow_.xt.scalar - workflow_.state.station1.param.xz) -
              0.15 * (workflow_.state.station2.param.xz -
                      workflow_.xt.scalar)) /
                 workflow_.state.station1.param.tz;
    }

    if (wake) {
      const double cst =
          0.03 * (workflow_.state.station2.param.xz -
                  workflow_.state.station1.param.xz) /
          workflow_.state.station1.param.tz;
      auto euler = [](double hk2, double hk1, double cst_local) {
        return hk2 - (hk2 + cst_local * std::pow(hk2 - 1, 3) - hk1) /
                         (1 + 3 * cst_local * std::pow(hk2 - 1, 2));
      };
      workflow_.state.station2.hkz.scalar = workflow_.state.station1.hkz.scalar;
      for (int i = 0; i < 3; ++i) {
        workflow_.state.station2.hkz.scalar =
            euler(workflow_.state.station2.hkz.scalar,
                  workflow_.state.station1.hkz.scalar, cst);
      }
      return workflow_.state.station2.hkz.scalar;
    }

    return workflow_.state.station1.hkz.scalar -
           0.15 * (workflow_.state.station2.param.xz -
                   workflow_.state.station1.param.xz) /
               workflow_.state.station1.param.tz;
  }

private:
  BoundaryLayerWorkflow &workflow_;
};
