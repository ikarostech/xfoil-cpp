#include "solver/boundary_layer/runtime/state.hpp"

#include <algorithm>
#include <sstream>

#include "numerics/spline.hpp"
#include "infrastructure/logger.hpp"

using Eigen::Vector2d;
using Eigen::VectorXd;

int BoundaryLayerRuntimeStateOps::resetSideState(
    SidePair<BoundaryLayerLattice> &lattice, BlTransitionParams &bl_transition,
    FlowRegimeEnum &flow_regime, int side, const Foil &foil,
    const StagnationResult &stagnation) {
  const int previousTransition = lattice.get(side).profiles.transitionIndex;
  bl_transition.xiforc = xifset(lattice, foil, stagnation, side);
  flow_regime = FlowRegimeEnum::Laminar;
  lattice.get(side).profiles.transitionIndex = lattice.get(side).trailingEdgeIndex;
  return previousTransition;
}

int BoundaryLayerRuntimeStateOps::readSideStationCount(
    const SidePair<BoundaryLayerLattice> &lattice, int side) {
  return lattice.get(side).stationCount;
}

BoundaryLayerStationReadModel
BoundaryLayerRuntimeStateOps::readStationModel(
    const SidePair<BoundaryLayerLattice> &lattice, const Eigen::VectorXd &wgap,
    int side, int stationIndex) {
  return lattice.get(side).readStationModel(wgap, stationIndex);
}

BoundaryLayerSideReadModel BoundaryLayerRuntimeStateOps::readSideModel(
    const SidePair<BoundaryLayerLattice> &lattice, int side) {
  return lattice.get(side).readSideModel();
}

BoundaryLayerTrailingEdgeReadModel
BoundaryLayerRuntimeStateOps::readTrailingEdgeModel(
    const SidePair<BoundaryLayerLattice> &lattice) {
  BoundaryLayerTrailingEdgeReadModel model;
  model.topTrailingEdgeIndex = lattice.top.trailingEdgeIndex;
  model.bottomTrailingEdgeIndex = lattice.bottom.trailingEdgeIndex;

  model.topMomentumThickness =
      lattice.top.profiles.momentumThickness[model.topTrailingEdgeIndex];
  model.bottomMomentumThickness =
      lattice.bottom.profiles.momentumThickness[model.bottomTrailingEdgeIndex];
  model.topDisplacementThickness =
      lattice.top.profiles.displacementThickness[model.topTrailingEdgeIndex];
  model.bottomDisplacementThickness =
      lattice.bottom.profiles.displacementThickness[model.bottomTrailingEdgeIndex];
  model.topSkinFrictionCoeff =
      lattice.top.profiles.skinFrictionCoeff[model.topTrailingEdgeIndex];
  model.bottomSkinFrictionCoeff =
      lattice.bottom.profiles.skinFrictionCoeff[model.bottomTrailingEdgeIndex];
  return model;
}

double BoundaryLayerRuntimeStateOps::xifset(
    const SidePair<BoundaryLayerLattice> &lattice, const Foil &foil,
    const StagnationResult &stagnation, int side) {
  VectorXd w1 = VectorXd::Zero(foil.foil_shape.n);
  double str;

  if (lattice.get(side).transitionLocation >= 1.0) {
    return lattice.get(side).arcLengthCoordinates[lattice.get(side).trailingEdgeIndex];
  }

  Vector2d point_chord = foil.edge.point_te - foil.edge.point_le;

  for (int i = 0; i < foil.foil_shape.n; i++) {
    w1[i] = (foil.foil_shape.points.col(i) - foil.edge.point_le)
                .dot(point_chord.normalized());
  }

  VectorXd w3 = spline::splind(w1, foil.foil_shape.spline_length);
  if (side == 1) {
    str = foil.edge.sle + (foil.foil_shape.spline_length[0] - foil.edge.sle) *
                              lattice.top.transitionLocation;
  } else {
    str =
        foil.edge.sle +
        (foil.foil_shape.spline_length[foil.foil_shape.n - 1] - foil.edge.sle) *
            lattice.bottom.transitionLocation;
  }
  str = spline::sinvrt(str, lattice.get(side).transitionLocation, w1,
                       w3, foil.foil_shape.spline_length, foil.foil_shape.n);
  double xiforc = std::min(
      (str - stagnation.sst),
      lattice.get(side).arcLengthCoordinates[lattice.get(side).trailingEdgeIndex]);
  if (xiforc < 0.0) {
    std::stringstream ss;
    ss << " ***  stagnation point is past trip on side " << side << "\n";
    Logger::instance().write(ss.str());
    return lattice.get(side).arcLengthCoordinates[lattice.get(side).trailingEdgeIndex];
  }

  return xiforc;
}
