#include "solver/xfoil/viscous_update.hpp"

#include <cmath>
#include <numbers>

#include "numerics/math_util.hpp"

ViscousUpdateResult BoundaryLayerViscousUpdate::run(
    const ViscousUpdateInput &input,
    const BoundaryLayerMatrix3x2dVector &vdel) {
  ViscousUpdateResult result;
  result.analysis_state = input.analysisState;
  result.aero_coeffs = input.aeroCoeffs;

  const auto ue_distribution =
      input.boundaryLayer.computeNewUeDistribution(input.aerodynamicContext,
                                                   vdel);
  const auto cl_contributions =
      input.boundaryLayer.computeClFromEdgeVelocityDistribution(
          input.aerodynamicContext, ue_distribution);

  const double cl_target = input.analysisState.controlByAlpha
                               ? input.aeroCoeffs.cl
                               : input.analysisState.clspec;
  const double dac = computeAcChange(
      cl_contributions.cl, input.aeroCoeffs.cl, cl_target,
      cl_contributions.cl_ac, cl_contributions.cl_a, cl_contributions.cl_ms,
      input.analysisState.controlByAlpha, input.analysisState.currentMach,
      input.machPerLift);

  double rlx = computeRelaxation(input.analysisState, input.aeroCoeffs, dac);

  double rmsbl = 0.0;
  const double dhi = 1.5;
  const double dlo = -0.5;
  const double gamma = 1.4;

  for (int side = 1; side <= 2; ++side) {
    const auto delta = input.boundaryLayer.buildBoundaryLayerDelta(
        side, ue_distribution.unew.get(side), ue_distribution.u_ac.get(side),
        dac, vdel);
    const auto metrics = input.boundaryLayer.evaluateSegmentRelaxation(
        side, delta, dhi, dlo, rlx);
    rmsbl += metrics.rmsContribution;

    const double hstinv =
        (gamma - 1) *
        MathUtil::pow(input.analysisState.currentMach / input.analysisState.qinf,
                      2) /
        (1.0 + 0.5 * (gamma - 1) * input.analysisState.currentMach *
                   input.analysisState.currentMach);

    result.profiles.get(side) = input.boundaryLayer.applyBoundaryLayerDelta(
        side, delta, rlx, hstinv, gamma - 1);
  }

  rmsbl = sqrt(rmsbl /
               (4.0 * double(input.boundaryLayer.readSideStationCount(1) +
                             input.boundaryLayer.readSideStationCount(2))));

  if (input.analysisState.controlByAlpha) {
    result.aero_coeffs.cl = input.aeroCoeffs.cl + rlx * dac;
  } else {
    result.analysis_state.alpha = input.analysisState.alpha + rlx * dac;
  }

  for (int kbl = 1;
       kbl <= input.boundaryLayer.readSideStationCount(2) -
                  (input.boundaryLayer.trailingEdgeIndex(2) + 1);
       ++kbl) {
    const int top_index = input.boundaryLayer.trailingEdgeIndex(1) + kbl;
    const int bottom_index = input.boundaryLayer.trailingEdgeIndex(2) + kbl;
    result.profiles.top.skinFrictionCoeff[top_index] =
        result.profiles.bottom.skinFrictionCoeff[bottom_index];
    result.profiles.top.momentumThickness[top_index] =
        result.profiles.bottom.momentumThickness[bottom_index];
    result.profiles.top.displacementThickness[top_index] =
        result.profiles.bottom.displacementThickness[bottom_index];
    result.profiles.top.edgeVelocity[top_index] =
        result.profiles.bottom.edgeVelocity[bottom_index];
    result.profiles.top.skinFrictionCoeffHistory[top_index] =
        result.profiles.bottom.skinFrictionCoeffHistory[bottom_index];
  }

  result.rmsbl = rmsbl;
  return result;
}

double BoundaryLayerViscousUpdate::computeAcChange(
    double clnew, double cl_current, double cl_target, double cl_ac,
    double cl_a, double cl_ms, bool controlByAlpha, double currentMach,
    double machPerLift) {
  if (controlByAlpha) {
    return (clnew - cl_current) /
           (1.0 - cl_ac - cl_ms * 2.0 * currentMach * machPerLift);
  }
  return (clnew - cl_target) / (0.0 - cl_ac - cl_a);
}

double BoundaryLayerViscousUpdate::computeRelaxation(
    const FlowState &analysis_state, const AeroCoefficients &aero_coeffs,
    double dac) {
  const double dtor = std::numbers::pi / 180.0;
  double dalmax = 0.5 * dtor;
  double dalmin = -0.5 * dtor;
  double dclmax = 0.5;
  double dclmin = -0.5;

  if (analysis_state.machType != FlowState::MachType::CONSTANT) {
    dclmin = std::max(-0.5, -0.9 * aero_coeffs.cl);
  }

  auto clampRelaxationForGlobalChange = [&](double relaxation, double lower,
                                            double upper) {
    if (dac == 0.0) {
      return relaxation;
    }
    if (relaxation * dac > upper) {
      relaxation = upper / dac;
    }
    if (relaxation * dac < lower) {
      relaxation = lower / dac;
    }
    return relaxation;
  };

  if (analysis_state.controlByAlpha) {
    return clampRelaxationForGlobalChange(1.0, dclmin, dclmax);
  }

  return clampRelaxationForGlobalChange(1.0, dalmin, dalmax);
}
