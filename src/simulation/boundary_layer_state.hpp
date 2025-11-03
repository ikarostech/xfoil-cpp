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
  SidePair<Eigen::VectorXi> stationToPanel;
  SidePair<Eigen::VectorXi> stationToSystem;
  SidePair<int> trailingEdgeIndex;
  SidePair<int> stationCount;

  SidePair<double> transitionLocation;

  SidePair<Eigen::VectorXd> ctau;
  SidePair<Eigen::VectorXd> thet;
  SidePair<Eigen::VectorXd> dstr;
  SidePair<Eigen::VectorXd> uedg;
  SidePair<Eigen::VectorXd> ctq;
  SidePair<Eigen::VectorXd> xssi;
  SidePair<Eigen::VectorXd> uinv;
  SidePair<Eigen::VectorXd> uinv_a;
  SidePair<Eigen::VectorXd> mass;
  SidePair<Eigen::VectorXd> vti;

  SidePair<int> transitionIndex;

  void clear() {
    stationToPanel.top.resize(0);
    stationToPanel.bottom.resize(0);
    stationToSystem.top.resize(0);
    stationToSystem.bottom.resize(0);
    trailingEdgeIndex.top = 0;
    trailingEdgeIndex.bottom = 0;
    stationCount.top = 0;
    stationCount.bottom = 0;
    transitionLocation.top = 0.0;
    transitionLocation.bottom = 0.0;
    transitionIndex.top = 0;
    transitionIndex.bottom = 0;

    auto resetVectorPair = [](SidePair<Eigen::VectorXd>& pair) {
      pair.top.resize(0);
      pair.bottom.resize(0);
    };

    resetVectorPair(ctau);
    resetVectorPair(thet);
    resetVectorPair(dstr);
    resetVectorPair(uedg);
    resetVectorPair(ctq);
    resetVectorPair(xssi);
    resetVectorPair(uinv);
    resetVectorPair(uinv_a);
    resetVectorPair(mass);
    resetVectorPair(vti);
  }

  void resizeSide(int side, int size) {
    stationToPanel.get(side).resize(size);
    stationToSystem.get(side).resize(size);
    ctau.get(side).resize(size);
    thet.get(side).resize(size);
    dstr.get(side).resize(size);
    uedg.get(side).resize(size);
    ctq.get(side).resize(size);
    xssi.get(side).resize(size);
    uinv.get(side).resize(size);
    uinv_a.get(side).resize(size);
    mass.get(side).resize(size);
    vti.get(side).resize(size);
  }
};
