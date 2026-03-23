#pragma once

#include <string_view>

#include "infrastructure/logger.hpp"
#include "solver/boundary_layer/boundary_layer.hpp"
#include "solver/boundary_layer/initialization/setbl_access.hpp"
#include "solver/boundary_layer/runtime/state.hpp"
#include "solver/march/mrchue_linear_system.hpp"

class BoundaryLayerMarchAccess {
  public:
    explicit BoundaryLayerMarchAccess(BoundaryLayer &boundaryLayer) : boundaryLayer_(boundaryLayer) {}

    Eigen::VectorXd &wgap() const {
        return boundaryLayer_.wakeGap();
    }
    BlSystemCoeffs &blc() const {
        return boundaryLayer_.systemCoefficients();
    }
    FlowRegimeEnum &flowRegime() const {
        return boundaryLayer_.flowRegime();
    }
    BlCompressibilityParams &blCompressibility() const {
        return boundaryLayer_.compressibility();
    }
    BlReynoldsParams &blReynolds() const {
        return boundaryLayer_.reynolds();
    }
    blDiff &xt() const {
        return boundaryLayer_.transitionSensitivity();
    }

    double readCurrentShapeFactor() const {
        return boundaryLayer_.readCurrentShapeFactor();
    }
    void refreshCurrentStationState(double xsi, double ami, double cti, double thi, double dsi, double dswaki,
                                    double uei) const {
        boundaryLayer_.refreshCurrentStationState(xsi, ami, cti, thi, dsi, dswaki, uei);
    }
    int readSideStationCount(int side) const {
        return boundaryLayer_.readSideStationCount(side);
    }
    BoundaryLayerStationReadModel readStationModel(int side, int stationIndex) const {
        return boundaryLayer_.readStationModel(side, stationIndex);
    }
    BoundaryLayerTrailingEdgeReadModel readTrailingEdgeModel() const {
        return boundaryLayer_.readTrailingEdgeModel();
    }

    void emitMarchInfoLog(std::string_view message) const {
        Logger::instance().write(std::string(message));
    }

    double readNewtonRhs(int row) const {
        return boundaryLayer_.readNewtonRhs(row);
    }
    void solveMrchueDirectNewtonSystem() const {
        boundaryLayer_.solveDirectNewtonSystem();
    }
    void solveMrchueInverseNewtonSystem(double htarg) const {
        boundaryLayer_.solveInverseNewtonSystem(htarg);
    }
    void runTransitionCheckForMrchue(int side, int stationIndex, double &ami, double &cti) const {
        boundaryLayer_.runTransitionCheckForMrchue(side, stationIndex, ami, cti);
    }
    double calcHtarg(int stationIndex, int side, bool wake) const {
        return boundaryLayer_.calcHtarg(stationIndex, side, wake);
    }

    bool isStartOfWake(int side, int stationIndex) const {
        return boundaryLayer_.isStartOfWake(side, stationIndex);
    }
    FlowRegimeEnum applyFlowRegimeCandidate(FlowRegimeEnum candidate) const {
        boundaryLayer_.flowRegime() = candidate;
        return boundaryLayer_.flowRegime();
    }
    FlowRegimeEnum currentFlowRegime() const {
        return boundaryLayer_.flowRegime();
    }
    FlowRegimeEnum determineRegimeForStation(int side, int stationIndex) const {
        return boundaryLayer_.determineRegimeForStation(side, stationIndex);
    }
    int resetSideState(int side, const Foil &foil, const StagnationResult &stagnation) const {
        return boundaryLayer_.resetSideState(side, foil, stagnation);
    }

    bool tesys(const BoundaryLayerSideState &top_profiles, const BoundaryLayerSideState &bottom_profiles,
               const Edge &edge) const {
        return boundaryLayer_.tesys(top_profiles, bottom_profiles, edge);
    }
    bool blsys() const {
        return boundaryLayer_.blsys();
    }

    void updateSystemMatricesForStation(const Edge &edge, int side, int stationIndex,
                                        BoundaryLayerMixedModeStationContext &ctx) const {
        boundaryLayer_.updateSystemMatricesForStation(edge, side, stationIndex, ctx);
    }
    void initializeFirstIterationState(int side, int stationIndex, int previousTransition,
                                       BoundaryLayerMixedModeStationContext &ctx, double &ueref, double &hkref) const {
        boundaryLayer_.initializeFirstIterationState(side, stationIndex, previousTransition, ctx, ueref, hkref);
    }
    void configureSimilarityRow(double ueref) const {
        boundaryLayer_.configureSimilarityRow(ueref);
    }
    void configureViscousRow(double hkref, double ueref, double senswt, bool resetSensitivity, bool averageSensitivity,
                             double &sens, double &sennew) const {
        boundaryLayer_.configureViscousRow(hkref, ueref, senswt, resetSensitivity, averageSensitivity, sens, sennew);
    }
    bool applyMixedModeNewtonStep(int side, int stationIndex, double &ami,
                                  BoundaryLayerMixedModeStationContext &ctx) const {
        return boundaryLayer_.applyMixedModeNewtonStep(side, stationIndex, ami, ctx);
    }
    void checkTransitionIfNeeded(int side, int stationIndex, bool skipCheck, int laminarAdvance, double &ami) const {
        boundaryLayer_.checkTransitionIfNeeded(side, stationIndex, skipCheck, laminarAdvance, ami);
    }
    void storeStationStateCommon(int side, int stationIndex, const BoundaryLayerMixedModeStationContext &ctx) const {
        boundaryLayer_.storeStationStateCommon(side, stationIndex, ctx);
    }
    void recoverStationAfterFailure(int side, int stationIndex, BoundaryLayerMixedModeStationContext &ctx, double &ami,
                                    BoundaryLayerEdgeVelocityFallbackMode edgeMode, int laminarAdvance) const {
        boundaryLayer_.recoverStationAfterFailure(side, stationIndex, ctx, ami, edgeMode, laminarAdvance);
    }
    bool solveTeSystemForCurrentProfiles(const Edge &edge) const {
        return boundaryLayer_.solveTeSystemForCurrentProfiles(edge);
    }

  private:
    BoundaryLayer &boundaryLayer_;
};
