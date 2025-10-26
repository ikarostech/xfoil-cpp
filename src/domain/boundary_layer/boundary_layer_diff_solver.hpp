#pragma once

#include <Eigen/Core>

#include "domain/coefficient/bl_newton.hpp"
#include "domain/flow_regime.hpp"
#include "simulation/boundary_layer_state.hpp"
#include "simulation/skin_friction_coefficients.hpp"

class BlDiffSolver {
    public:
    BlDiffSolver() = default;
    BlSystemCoeffs solve(FlowRegimeEnum flowRegimeType,
                         BoundaryLayerState boundaryLayerState,
                         const SkinFrictionCoefficients& skinFriction,
                         double amcrit);
    private:
    void bldifLaminar(BoundaryLayerState& boundaryLayerState, double amcrit,
                      BlSystemCoeffs& coeffs);

    void bldifTurbulent(BoundaryLayerState& boundaryLayerState, FlowRegimeEnum flowRegimeType,
                        double upw, const Eigen::Vector3d& upw1, const Eigen::Vector3d& upw2,
                        double upw_ms, double ulog, BlSystemCoeffs& coeffs);

    void bldifMomentum(BoundaryLayerState& boundaryLayerState, double xlog, double ulog,
                       double tlog, double ddlog, const SkinFrictionCoefficients& skinFriction,
                       BlSystemCoeffs& coeffs);

    void bldifShape(BoundaryLayerState& boundaryLayerState, double upw, double xlog, double ulog,
                    double hlog, double ddlog, const Eigen::Vector3d& upw1,
                    const Eigen::Vector3d& upw2, double upw_ms, BlSystemCoeffs& coeffs);
};
