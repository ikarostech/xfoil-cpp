#pragma once

#include <stdexcept>
#include <utility>

#include <Eigen/Core>

#include "core/side_pair.hpp"
#include "domain/boundary_layer.hpp"

/**
 * @brief Bundles the boundary-layer state at the two stations used by the
 *        marching routines (historically known as "1" and "2").
 */
struct BoundaryLayerState {
  blData station1;
  blData station2;

  blData& previous() { return station1; }
  const blData& previous() const { return station1; }

  blData& current() { return station2; }
  const blData& current() const { return station2; }

  blData& station(int index) {
    switch (index) {
      case 1:
        return station1;
      case 2:
        return station2;
      default:
        throw std::invalid_argument("BoundaryLayerState::station expects 1 or 2");
    }
  }

  const blData& station(int index) const {
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
};


/**
 * @brief Aggregates per-side arrays and indices that describe the boundary-layer
 *        lattice along the airfoil and wake.
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
  BoundaryLayerLattice(const int size) { resize(size); }
};
