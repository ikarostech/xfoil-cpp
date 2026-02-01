#pragma once

#include "domain/boundary_layer/boundary_layer_diff_solver.hpp"
#include "domain/boundary_layer/boundary_layer_variables_solver.hpp"
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
  BoundaryLayerVariablesSolver boundaryLayerVariablesSolver;
  BlDiffSolver blDiffSolver;
};
