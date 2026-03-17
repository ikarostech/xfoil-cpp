#pragma once

#include "solver/boundary_layer/boundary_layer_diff_solver.hpp"
#include "solver/boundary_layer/boundary_layer_variables_solver.hpp"
#include "model/boundary_layer/bl_compressibility_params.hpp"
#include "model/boundary_layer/bl_reynolds_params.hpp"
#include "model/boundary_layer/bl_transition_params.hpp"
#include "solver/boundary_layer/boundary_layer_store.hpp"
#include "solver/boundary_layer/skin_friction_coefficients.hpp"

class BoundaryLayerWorkflow;
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
  bool blkin(BoundaryLayerState &state) const;
  SkinFrictionCoefficients blmid(FlowRegimeEnum flowRegimeType) const;

  Context context_;
  BoundaryLayerStore boundaryLayerStore;
  BoundaryLayerVariablesSolver boundaryLayerVariablesSolver;
  BlDiffSolver blDiffSolver;
};
