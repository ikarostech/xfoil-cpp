#pragma once

#include "model/coefficient/aero_coefficients.hpp"
#include "model/flow_state.hpp"
#include "solver/boundary_layer/boundary_layer.hpp"
#include "solver/boundary_layer/boundary_layer_aerodynamics.hpp"

struct ViscousUpdateResult {
    double rmsbl = 0.0;
    FlowState analysis_state;
    AeroCoefficients aero_coeffs;
    SidePair<BoundaryLayerSideState> profiles;
};

struct ViscousUpdateInput {
    const BoundaryLayer &boundaryLayer;
    const FlowState &analysisState;
    const AeroCoefficients &aeroCoeffs;
    double machPerLift = 0.0;
    BoundaryLayerAerodynamicContext aerodynamicContext;
};

class BoundaryLayerViscousUpdate {
  public:
    static ViscousUpdateResult run(const ViscousUpdateInput &input, const BoundaryLayerMatrix3x2dVector &vdel);

  private:
    static double computeAcChange(double clnew, double cl_current, double cl_target, double cl_ac, double cl_a,
                                  double cl_ms, bool controlByAlpha, double currentMach, double machPerLift);
    static double computeRelaxation(const FlowState &analysis_state, const AeroCoefficients &aero_coeffs, double dac);
};
