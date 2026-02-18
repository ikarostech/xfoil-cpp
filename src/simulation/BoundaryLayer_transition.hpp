#pragma once

#include "domain/boundary_layer/boundary_layer_diff_solver.hpp"
#include "domain/boundary_layer/boundary_layer_variables_solver.hpp"
#include "simulation/BoundaryLayerStore.hpp"

class BoundaryLayerWorkflow;
class XFoil;

class BoundaryLayerTransitionSolver {
 public:
  explicit BoundaryLayerTransitionSolver(BoundaryLayerWorkflow& workflow);

  bool trchek();
  bool trdif();

 private:
  struct TrchekData;
  struct TrdifData;

  double computeTransitionLocation(double weightingFactor) const;
  bool iterateAmplification(TrchekData& data);
  bool resolveTransitionLocationAndSensitivities(TrchekData& data);
  void setupLaminarTransitionSystem(TrdifData& data);
  void setupTurbulentTransitionSystem(TrdifData& data);
  void mergeTransitionSystems(TrdifData& data);

  BoundaryLayerWorkflow* workflow_;
  BoundaryLayerStore boundaryLayerStore;
  BoundaryLayerVariablesSolver boundaryLayerVariablesSolver;
  BlDiffSolver blDiffSolver;
};
