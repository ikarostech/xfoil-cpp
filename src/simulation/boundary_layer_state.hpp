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
 * @brief Aggregates per-side arrays and indices that describe the boundary-layer
 *        lattice along the airfoil and wake.
 */
struct BoundaryLayerLattice {
  Eigen::VectorXi stationToPanel;
  Eigen::VectorXi stationToSystem;
  int trailingEdgeIndex = 0;
  int stationCount = 0;

  double transitionLocation = 0.0;

  Eigen::VectorXd ctau;
  Eigen::VectorXd thet;
  Eigen::VectorXd dstr;
  Eigen::VectorXd uedg;
  Eigen::VectorXd ctq;
  Eigen::VectorXd xssi;
  Eigen::VectorXd uinv;
  Eigen::VectorXd uinv_a;
  Eigen::VectorXd mass;
  Eigen::VectorXd vti;

  int transitionIndex = 0;

  void clear() {
    stationToPanel.resize(0);
    stationToSystem.resize(0);
    trailingEdgeIndex = 0;
    stationCount = 0;
    transitionLocation = 0.0;
    transitionIndex = 0;

    auto resetVector = [](Eigen::VectorXd& vector) { vector.resize(0); };

    resetVector(ctau);
    resetVector(thet);
    resetVector(dstr);
    resetVector(uedg);
    resetVector(ctq);
    resetVector(xssi);
    resetVector(uinv);
    resetVector(uinv_a);
    resetVector(mass);
    resetVector(vti);
  }

  void resize(int size) {
    stationToPanel.resize(size);
    stationToSystem.resize(size);
    ctau.resize(size);
    thet.resize(size);
    dstr.resize(size);
    uedg.resize(size);
    ctq.resize(size);
    xssi.resize(size);
    uinv.resize(size);
    uinv_a.resize(size);
    mass.resize(size);
    vti.resize(size);
  }
};
