#pragma once

#include "model/boundary_layer/reference/bl_compressibility_params.hpp"
#include "model/boundary_layer/reference/bl_reynolds_params.hpp"
#include "model/boundary_layer/reference/bl_transition_params.hpp"
#include "model/boundary_layer/skin_friction_coefficients.hpp"
#include "model/boundary_layer/state.hpp"
#include "model/coefficient/aero_coefficients.hpp"
#include "model/flow_regime.hpp"
#include "model/flow_state.hpp"

struct BoundaryLayerReferenceParams {
    double currentMach = 0.0;
    double currentRe   = 0.0;
    double re_clmr     = 0.0;
    double msq_clmr    = 0.0;
    BlCompressibilityParams blCompressibility;
    BlReynoldsParams blReynolds{};
    double amcrit = 0.0;
};

class BoundaryLayerPhysics {
  public:
    static constexpr double kHvrat = 0.35;

    static bool blkin(BoundaryLayerStationWindow &state, const BlCompressibilityParams &compressibility,
                      const BlReynoldsParams &reynolds);
    static void refreshCurrentStation(BoundaryLayerStationWindow &state, const BlCompressibilityParams &compressibility,
                                      const BlReynoldsParams &reynolds, double xsi, double ami, double cti, double thi,
                                      double dsi, double dswaki, double uei);
    static SkinFrictionCoefficients blmid(BoundaryLayerStationWindow &state, FlowRegimeEnum flowRegimeType);
    static BoundaryLayerStationState blprv(BoundaryLayerStationState data,
                                           const BlCompressibilityParams &compressibility, double xsi, double ami,
                                           double cti, double thi, double dsi, double dswaki, double uei);
    static double adjustDisplacementForHkLimit(double displacementThickness, double momentumThickness, double msq,
                                               double hklim);
    static BoundaryLayerReferenceParams buildReferenceParams(const FlowState &analysisState,
                                                             const AeroCoefficients &aeroCoefficients, double acrit);
};
