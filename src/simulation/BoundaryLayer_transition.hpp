#pragma once

#include "simulation/BoundaryLayerStore.hpp"

class BoundaryLayerWorkflow;
class XFoil;

class BoundaryLayerTransitionSolver {
 public:
  explicit BoundaryLayerTransitionSolver(BoundaryLayerWorkflow& workflow);

  bool trchek(XFoil& xfoil);
  bool trdif();

 private:
  double computeTransitionLocation(double weightingFactor) const;

  BoundaryLayerWorkflow* workflow_;
  BoundaryLayerStore boundaryLayerStore;
};
