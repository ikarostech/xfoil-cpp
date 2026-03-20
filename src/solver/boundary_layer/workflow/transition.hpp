#pragma once

#include "model/boundary_layer/bl_compressibility_params.hpp"
#include "model/boundary_layer/physics.hpp"
#include "model/boundary_layer/bl_reynolds_params.hpp"
#include "model/boundary_layer/bl_transition_params.hpp"
#include "model/boundary_layer/state.hpp"
#include "numerics/boundary_layer/diff_system.hpp"
#include "numerics/boundary_layer/variables.hpp"
#include "solver/boundary_layer/boundary_layer_store.hpp"

class BoundaryLayer;
class XFoil;

class BoundaryLayerTransitionSolver {
public:
  struct Context {
    BoundaryLayerState &state;
    BlTransitionParams &blTransition;
    FlowRegimeEnum &flowRegime;
    blDiff &xt;
    BlSystemCoeffs &blc;
    const BlCompressibilityParams &blCompressibility;
    const BlReynoldsParams &blReynolds;
  };

  explicit BoundaryLayerTransitionSolver(Context context);

  bool trchek();
  bool trdif();

private:
  struct TrchekData;
  struct TrdifData;

  double computeTransitionLocation(const BoundaryLayerState &state,
                                   double weightingFactor) const;
  bool iterateAmplification(TrchekData &data);
  bool resolveTransitionLocationAndSensitivities(TrchekData &data);
  void setupLaminarTransitionSystem(TrdifData &data);
  void setupTurbulentTransitionSystem(TrdifData &data);
  void mergeTransitionSystems(TrdifData &data);

  Context context_;
  BoundaryLayerStore boundaryLayerStore;
  BoundaryLayerVariablesSolver boundaryLayerVariablesSolver;
  BlDiffSolver blDiffSolver;
};
