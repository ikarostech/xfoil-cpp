#include "solver/boundary_layer/viscous_iteration_ops.hpp"

#include "solver/boundary_layer/blsolve.hpp"
#include "solver/boundary_layer/initialization/viscous_initializer.hpp"

namespace boundary_layer_iteration_ops {

SetblOutputView initializeNewtonSystem(BoundaryLayer &boundary_layer,
                                       const FlowState &analysis_state,
                                       const AeroCoefficients &aero_coeffs,
                                       double acrit, const Foil &foil,
                                       const StagnationResult &stagnation,
                                       const Eigen::MatrixXd &dij,
                                       bool bl_initialized) {
  auto setbl_output = BoundaryLayerInitializer::run(
      boundary_layer, analysis_state, aero_coeffs, acrit, foil, stagnation, dij,
      bl_initialized);
  BoundaryLayerInitializer::applyOutput(boundary_layer, setbl_output);
  return setbl_output;
}

BoundaryLayerMatrix3x2dVector solveNewtonStep(
    const BoundaryLayer &boundary_layer, double vaccel,
    const SetblOutputView &setbl_output) {
  Blsolve solver;
  SidePair<int> ivte{boundary_layer.trailingEdgeSystemIndex(1),
                     boundary_layer.trailingEdgeSystemIndex(2)};
  return solver.solve(boundary_layer.systemSize(), ivte, vaccel,
                      setbl_output.bl_newton_system);
}

void applyIterationUpdate(BoundaryLayer &boundary_layer,
                          const ViscousUpdateResult &update_result) {
  boundary_layer.applyProfiles(update_result.profiles);
}

} // namespace boundary_layer_iteration_ops
