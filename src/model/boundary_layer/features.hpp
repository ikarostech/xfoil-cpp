#pragma once

struct BoundaryLayerStagnationFeature {
  int index = 0;
  double sst = 0.0;
  double sst_go = 0.0;
  double sst_gp = 0.0;
  bool found = true;
};

struct BoundaryLayerTrailingEdgeFeature {
  int topIndex = 0;
  int bottomIndex = 0;
  double topMomentumThickness = 0.0;
  double bottomMomentumThickness = 0.0;
  double topDisplacementThickness = 0.0;
  double bottomDisplacementThickness = 0.0;
  double topSkinFrictionCoeff = 0.0;
  double bottomSkinFrictionCoeff = 0.0;
};
