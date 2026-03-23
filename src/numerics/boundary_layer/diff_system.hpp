#pragma once

#include <Eigen/Core>

#include "model/boundary_layer/skin_friction_coefficients.hpp"
#include "model/boundary_layer/state.hpp"
#include "model/flow_regime.hpp"
#include "numerics/coefficient/bl_newton.hpp"

class BlDiffSolver {
  public:
    BlDiffSolver() = default;
    BlSystemCoeffs solve(FlowRegimeEnum flowRegimeType, BoundaryLayerStationWindow BoundaryLayerStationWindow,
                         const SkinFrictionCoefficients &skinFriction, double amcrit);

  private:
    void bldifLaminar(BoundaryLayerStationWindow &BoundaryLayerStationWindow, double amcrit, BlSystemCoeffs &coeffs);

    void bldifTurbulent(BoundaryLayerStationWindow &BoundaryLayerStationWindow, FlowRegimeEnum flowRegimeType,
                        double upw, const Eigen::Vector3d &upw1, const Eigen::Vector3d &upw2, double upw_ms,
                        double ulog, BlSystemCoeffs &coeffs);

    void bldifMomentum(BoundaryLayerStationWindow &BoundaryLayerStationWindow, double xlog, double ulog, double tlog,
                       double ddlog, const SkinFrictionCoefficients &skinFriction, BlSystemCoeffs &coeffs);

    void bldifShape(BoundaryLayerStationWindow &BoundaryLayerStationWindow, double upw, double xlog, double ulog,
                    double hlog, double ddlog, const Eigen::Vector3d &upw1, const Eigen::Vector3d &upw2, double upw_ms,
                    BlSystemCoeffs &coeffs);
};
