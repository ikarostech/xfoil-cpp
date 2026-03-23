#pragma once

#include <cmath>
#include <stdexcept>
#include <utility>

#include <Eigen/Core>

#include "model/boundary_layer.hpp"
#include "model/boundary_layer/runtime_types.hpp"

/**
 * @brief Model-owned two-station window used for local boundary-layer updates.
 *
 * The boundary-layer model advances by operating on a previous/current station
 * pair. This class keeps that pair cohesive and hides the station bookkeeping
 * behind a small set of operations.
 */
class BoundaryLayerStationWindow {
  public:
    BoundaryLayerStationState station1;
    BoundaryLayerStationState station2;

    BoundaryLayerStationState &previous() {
        return station1;
    }
    const BoundaryLayerStationState &previous() const {
        return station1;
    }

    BoundaryLayerStationState &current() {
        return station2;
    }
    const BoundaryLayerStationState &current() const {
        return station2;
    }

    BoundaryLayerStationState &station(int index) {
        switch (index) {
        case 1:
            return station1;
        case 2:
            return station2;
        default:
            throw std::invalid_argument("BoundaryLayerStationWindow::station expects 1 or 2");
        }
    }

    const BoundaryLayerStationState &station(int index) const {
        switch (index) {
        case 1:
            return station1;
        case 2:
            return station2;
        default:
            throw std::invalid_argument("BoundaryLayerStationWindow::station expects 1 or 2");
        }
    }

    void swapStations() {
        std::swap(station1, station2);
    }
    void advance() {
        station1 = station2;
    }
    void stepbl() {
        advance();
    }
};

/**
 * @brief Base state for one boundary-layer side.
 *
 * This corresponds to the FoilShape role in the foil model: it stores the
 * primary side-wise arrays and exposes a few derived queries used by the rest
 * of the boundary-layer model.
 */
class BoundaryLayerSideState {
  public:
    BoundaryLayerSideState() = default;
    explicit BoundaryLayerSideState(int size)
        : edgeVelocity(size), skinFrictionCoeff(size), momentumThickness(size), displacementThickness(size),
          massFlux(size), skinFrictionCoeffHistory(size) {}

    Eigen::VectorXd edgeVelocity;
    Eigen::VectorXd skinFrictionCoeff;
    Eigen::VectorXd momentumThickness;
    Eigen::VectorXd displacementThickness;
    Eigen::VectorXd massFlux;
    Eigen::VectorXd skinFrictionCoeffHistory;
    int transitionIndex = 0;

    void resize(int size) {
        skinFrictionCoeff.resize(size);
        momentumThickness.resize(size);
        displacementThickness.resize(size);
        edgeVelocity.resize(size);
        massFlux.resize(size);
        skinFrictionCoeffHistory.resize(size);
    }

    void zeroStateVectors() {
        if (edgeVelocity.size() > 0)
            edgeVelocity.setZero();
        if (skinFrictionCoeff.size() > 0)
            skinFrictionCoeff.setZero();
        if (momentumThickness.size() > 0)
            momentumThickness.setZero();
        if (displacementThickness.size() > 0)
            displacementThickness.setZero();
        if (massFlux.size() > 0)
            massFlux.setZero();
    }

    bool hasFiniteThicknessForStationCount(int stationCount) const {
        if (stationCount <= 1 || momentumThickness.size() < stationCount ||
            displacementThickness.size() < stationCount) {
            return false;
        }

        const int count  = stationCount - 1;
        const auto theta = momentumThickness.head(count);
        const auto delta = displacementThickness.head(count);
        return theta.allFinite() && delta.allFinite() && theta.maxCoeff() > 0.0 && delta.maxCoeff() > 0.0;
    }
};

/**
 * @brief Model base state for one side's boundary-layer lattice.
 *
 * The lattice combines topology and the side state arrays for one boundary
 * layer side. It stays in the model because it represents domain state, not a
 * transient solver-owned container.
 */
class BoundaryLayerLattice {
  public:
    BoundaryLayerLattice() = default;
    explicit BoundaryLayerLattice(int size)
        : stationToPanel(size), stationToSystem(size), profiles(size), arcLengthCoordinates(size),
          inviscidEdgeVelocityMatrix(2, size), panelInfluenceFactor(size) {}

    Eigen::VectorXi stationToPanel;
    Eigen::VectorXi stationToSystem;
    int trailingEdgeIndex = 0;
    int stationCount      = 0;

    double transitionLocation = 0.0;

    BoundaryLayerSideState profiles;
    Eigen::VectorXd arcLengthCoordinates;
    Eigen::Matrix2Xd inviscidEdgeVelocityMatrix;
    Eigen::VectorXd panelInfluenceFactor;

    void resize(int size) {
        stationToPanel.resize(size);
        stationToSystem.resize(size);
        profiles.resize(size);
        arcLengthCoordinates.resize(size);
        inviscidEdgeVelocityMatrix.resize(2, size);
        panelInfluenceFactor.resize(size);
    }

    int trailingEdgeSystemIndex() const {
        return stationToSystem[trailingEdgeIndex];
    }

    bool isStartOfWake(int stationIndex) const {
        return stationIndex == trailingEdgeIndex + 1;
    }

    BoundaryLayerStationReadModel readStationModel(const Eigen::VectorXd &wakeGap, int stationIndex) const {
        BoundaryLayerStationReadModel model;
        model.stationCount          = stationCount;
        model.trailingEdgeIndex     = trailingEdgeIndex;
        model.transitionIndex       = profiles.transitionIndex;
        model.arcLength             = arcLengthCoordinates[stationIndex];
        model.edgeVelocity          = profiles.edgeVelocity[stationIndex];
        model.momentumThickness     = profiles.momentumThickness[stationIndex];
        model.displacementThickness = profiles.displacementThickness[stationIndex];
        model.skinFrictionCoeff     = profiles.skinFrictionCoeff[stationIndex];

        if (stationIndex > trailingEdgeIndex) {
            const int wakeIndex = stationIndex - trailingEdgeIndex;
            model.wakeGap       = wakeGap[wakeIndex - 1];
        }

        return model;
    }

    BoundaryLayerSideReadModel readSideModel() const {
        BoundaryLayerSideReadModel model;
        model.stationCount       = stationCount;
        model.trailingEdgeIndex  = trailingEdgeIndex;
        model.transitionIndex    = profiles.transitionIndex;
        model.hasStations        = stationCount > 1;
        model.hasFiniteThickness = profiles.hasFiniteThicknessForStationCount(stationCount);

        if (!model.hasStations) {
            return model;
        }

        const int lastIndex             = stationCount - 2;
        model.lastEdgeVelocity          = profiles.edgeVelocity[lastIndex];
        model.lastMomentumThickness     = profiles.momentumThickness[lastIndex];
        model.lastDisplacementThickness = profiles.displacementThickness[lastIndex];
        return model;
    }

    bool hasValidPanelMap(int totalNodes) const {
        if (stationCount <= 1 || stationToPanel.size() < stationCount || panelInfluenceFactor.size() < stationCount) {
            return false;
        }
        for (int i = 0; i < stationCount - 1; ++i) {
            const int panel = stationToPanel[i];
            if (panel < 0 || panel >= totalNodes) {
                return false;
            }
        }
        return trailingEdgeIndex >= 0 && trailingEdgeIndex < stationCount;
    }
};
