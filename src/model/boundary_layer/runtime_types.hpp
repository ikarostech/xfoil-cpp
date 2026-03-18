#pragma once

#include "model/flow_regime.hpp"

struct BoundaryLayerStationReadModel {
  int stationCount = 0;
  int trailingEdgeIndex = 0;
  int transitionIndex = 0;
  double arcLength = 0.0;
  double edgeVelocity = 0.0;
  double momentumThickness = 0.0;
  double displacementThickness = 0.0;
  double skinFrictionCoeff = 0.0;
  double wakeGap = 0.0;
};

struct BoundaryLayerTrailingEdgeReadModel {
  int topTrailingEdgeIndex = 0;
  int bottomTrailingEdgeIndex = 0;
  double topMomentumThickness = 0.0;
  double bottomMomentumThickness = 0.0;
  double topDisplacementThickness = 0.0;
  double bottomDisplacementThickness = 0.0;
  double topSkinFrictionCoeff = 0.0;
  double bottomSkinFrictionCoeff = 0.0;
};
