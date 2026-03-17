#pragma once

#include <cmath>

#include "solver/boundary_layer/march_access.hpp"
#include "solver/march/context.hpp"

class WorkflowMarchContext final : public MarchContext {
 public:
  explicit WorkflowMarchContext(BoundaryLayerMarchAccess access)
      : access_(access) {}

  BoundaryLayerState &mutableState() override { return access_.state(); }

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
    access_.transitionSolver().trchek();
    ami = access_.state().station2.param.amplz;
    if (access_.flowRegime() == FlowRegimeEnum::Transition) {
      access_.lattice().get(side).profiles.transitionIndex = stationIndex;
      if (cti <= 0.0) {
        cti = 0.03;
        access_.state().station2.param.sz = cti;
      }
    } else {
      access_.lattice().get(side).profiles.transitionIndex = stationIndex + 2;
    }
  }
  bool solveTeSystemForCurrentProfiles(const Edge &edge) override {
    return access_.tesys(access_.lattice().top.profiles,
                         access_.lattice().bottom.profiles, edge);
  }

  double readBlCompressibilityHstinv() const override {
    return access_.blCompressibility().hstinv;
  }
  double readBlCompressibilityGm1bl() const override {
    return access_.blCompressibility().gm1bl;
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

  blData blprv(blData data, double xsi, double ami, double cti, double thi,
               double dsi, double dswaki, double uei) const override {
    return access_.blprv(data, xsi, ami, cti, thi, dsi, dswaki, uei);
  }
  bool blkin(BoundaryLayerState &state) override { return access_.blkin(state); }
  bool blsys() override { return access_.blsys(); }
  double calcHtarg(int ibl, int is, bool wake) override {
    if (ibl < access_.lattice().get(is).profiles.transitionIndex) {
      return access_.state().station1.hkz.scalar +
             0.03 * (access_.state().station2.param.xz -
                     access_.state().station1.param.xz) /
                 access_.state().station1.param.tz;
    }

    if (ibl == access_.lattice().get(is).profiles.transitionIndex) {
      return access_.state().station1.hkz.scalar +
             (0.03 * (access_.xt().scalar - access_.state().station1.param.xz) -
              0.15 * (access_.state().station2.param.xz - access_.xt().scalar)) /
                 access_.state().station1.param.tz;
    }

    if (wake) {
      const double cst =
          0.03 * (access_.state().station2.param.xz -
                  access_.state().station1.param.xz) /
          access_.state().station1.param.tz;
      auto euler = [](double hk2, double hk1, double cst_local) {
        return hk2 - (hk2 + cst_local * std::pow(hk2 - 1, 3) - hk1) /
                         (1 + 3 * cst_local * std::pow(hk2 - 1, 2));
      };
      access_.state().station2.hkz.scalar = access_.state().station1.hkz.scalar;
      for (int i = 0; i < 3; ++i) {
        access_.state().station2.hkz.scalar =
            euler(access_.state().station2.hkz.scalar,
                  access_.state().station1.hkz.scalar, cst);
      }
      return access_.state().station2.hkz.scalar;
    }

    return access_.state().station1.hkz.scalar -
           0.15 * (access_.state().station2.param.xz -
                   access_.state().station1.param.xz) /
               access_.state().station1.param.tz;
  }

 private:
  BoundaryLayerMarchAccess access_;
};
