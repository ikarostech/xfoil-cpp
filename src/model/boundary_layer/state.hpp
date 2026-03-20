#pragma once

#include <cmath>
#include <stdexcept>
#include <utility>

#include <Eigen/Core>

#include "model/boundary_layer.hpp"
#include "model/boundary_layer/runtime_types.hpp"

/**
 * @brief Bundles the boundary-layer state at the two stations used by the
 *        marching routines (historically known as "1" and "2").
 */
struct BoundaryLayerState {
  blData station1;
  blData station2;

  blData &previous() { return station1; }
  const blData &previous() const { return station1; }

  blData &current() { return station2; }
  const blData &current() const { return station2; }

  blData &station(int index) {
    switch (index) {
    case 1:
      return station1;
    case 2:
      return station2;
    default:
      throw std::invalid_argument("BoundaryLayerState::station expects 1 or 2");
    }
  }

  const blData &station(int index) const {
    switch (index) {
    case 1:
      return station1;
    case 2:
      return station2;
    default:
      throw std::invalid_argument("BoundaryLayerState::station expects 1 or 2");
    }
  }

  void swapStations() { std::swap(station1, station2); }
  void stepbl() { station1 = station2; }
};

/**
 * @brief Holds the per-station vectors that describe a single side of the
 *        boundary layer.
 */
struct BoundaryLayerSideProfiles {
  Eigen::VectorXd edgeVelocity;
  Eigen::VectorXd skinFrictionCoeff;
  Eigen::VectorXd momentumThickness;
  Eigen::VectorXd displacementThickness;
  Eigen::VectorXd massFlux;
  Eigen::VectorXd skinFrictionCoeffHistory;
  int transitionIndex = 0;

  void clear() {
    skinFrictionCoeff.resize(0);
    momentumThickness.resize(0);
    displacementThickness.resize(0);
    edgeVelocity.resize(0);
    massFlux.resize(0);
    skinFrictionCoeffHistory.resize(0);
    transitionIndex = 0;
  }

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

    const int count = stationCount - 1;
    const auto theta = momentumThickness.head(count);
    const auto delta = displacementThickness.head(count);
    return theta.allFinite() && delta.allFinite() && theta.maxCoeff() > 0.0 &&
           delta.maxCoeff() > 0.0;
  }
};

/**
 * @brief Aggregates per-side arrays and indices that describe the
 * boundary-layer lattice along the airfoil and wake.
 */
struct BoundaryLayerLattice {
  Eigen::VectorXi stationToPanel;
  Eigen::VectorXi stationToSystem;
  int trailingEdgeIndex = 0;
  int stationCount = 0;

  double transitionLocation = 0.0;

  BoundaryLayerSideProfiles profiles;
  Eigen::VectorXd arcLengthCoordinates;
  Eigen::Matrix2Xd inviscidEdgeVelocityMatrix;
  Eigen::VectorXd panelInfluenceFactor;

  void clear() {
    stationToPanel.resize(0);
    stationToSystem.resize(0);
    trailingEdgeIndex = 0;
    stationCount = 0;
    transitionLocation = 0.0;

    profiles.clear();
    arcLengthCoordinates.resize(0);
    inviscidEdgeVelocityMatrix.resize(0, 0);
    panelInfluenceFactor.resize(0);
  }

  void resize(int size) {
    stationToPanel.resize(size);
    stationToSystem.resize(size);
    profiles.resize(size);
    arcLengthCoordinates.resize(size);
    inviscidEdgeVelocityMatrix.resize(2, size);
    panelInfluenceFactor.resize(size);
  }

  BoundaryLayerLattice() = default;
  explicit BoundaryLayerLattice(int size) { resize(size); }

  int trailingEdgeSystemIndex() const {
    return stationToSystem[trailingEdgeIndex];
  }

  bool isStartOfWake(int stationIndex) const {
    return stationIndex == trailingEdgeIndex + 1;
  }

  BoundaryLayerStationReadModel readStationModel(
      const Eigen::VectorXd &wakeGap, int stationIndex) const {
    BoundaryLayerStationReadModel model;
    model.stationCount = stationCount;
    model.trailingEdgeIndex = trailingEdgeIndex;
    model.transitionIndex = profiles.transitionIndex;
    model.arcLength = arcLengthCoordinates[stationIndex];
    model.edgeVelocity = profiles.edgeVelocity[stationIndex];
    model.momentumThickness = profiles.momentumThickness[stationIndex];
    model.displacementThickness = profiles.displacementThickness[stationIndex];
    model.skinFrictionCoeff = profiles.skinFrictionCoeff[stationIndex];

    if (stationIndex > trailingEdgeIndex) {
      const int wakeIndex = stationIndex - trailingEdgeIndex;
      model.wakeGap = wakeGap[wakeIndex - 1];
    }

    return model;
  }

  BoundaryLayerSideReadModel readSideModel() const {
    BoundaryLayerSideReadModel model;
    model.stationCount = stationCount;
    model.trailingEdgeIndex = trailingEdgeIndex;
    model.transitionIndex = profiles.transitionIndex;
    model.hasStations = stationCount > 1;
    model.hasFiniteThickness = profiles.hasFiniteThicknessForStationCount(stationCount);

    if (!model.hasStations) {
      return model;
    }

    const int lastIndex = stationCount - 2;
    model.lastEdgeVelocity = profiles.edgeVelocity[lastIndex];
    model.lastMomentumThickness = profiles.momentumThickness[lastIndex];
    model.lastDisplacementThickness = profiles.displacementThickness[lastIndex];
    return model;
  }

  bool hasValidPanelMap(int totalNodes) const {
    if (stationCount <= 1 || stationToPanel.size() < stationCount ||
        panelInfluenceFactor.size() < stationCount) {
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
