#include "model/coefficient/xfoil_postprocess.hpp"

#include <algorithm>
#include <cmath>

#include "numerics/math_util.hpp"

namespace xfoil_postprocess {

XFoilCompressibilityParams buildCompressibilityParams(
    const FlowState &analysis_state) {
  const double current_mach = analysis_state.currentMach;
  const double beta = std::sqrt(1.0 - current_mach * current_mach);
  const double beta_msq = -0.5 / beta;
  const double prandtlGlauertFactor =
      0.5 * current_mach * current_mach / (1.0 + beta);
  const double prandtlGlauertFactor_msq =
      0.5 / (1.0 + beta) - prandtlGlauertFactor / (1.0 + beta) * beta_msq;
  const double karmanTsienFactor =
      (current_mach / (1.0 + beta)) * (current_mach / (1.0 + beta));
  const double karmanTsienFactor_msq =
      1.0 / ((1.0 + beta) * (1.0 + beta)) -
      2.0 * karmanTsienFactor / (1.0 + beta) * beta_msq;
  return {beta,
          beta_msq,
          karmanTsienFactor,
          karmanTsienFactor_msq,
          prandtlGlauertFactor,
          prandtlGlauertFactor_msq};
}

XFoilPressureCoefficientResult computePressureCoefficient(
    double tangential_velocity, double velocity_derivative, double qinf,
    const XFoilCompressibilityParams &params) {
  const double velocity_ratio = tangential_velocity / qinf;
  const double cginc = 1.0 - velocity_ratio * velocity_ratio;
  const double denom = params.beta + params.prandtlGlauertFactor * cginc;
  const double pressure_coefficient = cginc / denom;
  const double cp_msq =
      -pressure_coefficient / denom *
      (params.beta_msq + params.prandtlGlauertFactor_msq * cginc);

  double cp_velocity_derivative = 0.0;
  if (velocity_derivative != 0.0) {
    const double cpi =
        -2.0 * tangential_velocity / (qinf * qinf);
    const double cpc_cpi =
        (1.0 - params.prandtlGlauertFactor * pressure_coefficient) / denom;
    cp_velocity_derivative = cpc_cpi * cpi * velocity_derivative;
  }

  return {pressure_coefficient, cp_msq, cp_velocity_derivative};
}

XFoilClComputation computeCl(const Foil &foil, const FlowState &analysis_state,
                             const XFoilInternalState::InviscidFieldState &state,
                             const Eigen::Vector2d &ref) {
  XFoilClComputation result;
  double xcp_accumulator = 0.0;

  const auto compressibility = buildCompressibilityParams(analysis_state);
  const Eigen::Matrix2d rotateMatrix =
      MathUtil::getRotateMatrix(analysis_state.alpha);
  const int point_count = foil.foil_shape.n;

  const auto cp_first = computePressureCoefficient(
      state.surfaceVortex(0, 0), state.surfaceVortex(1, 0), analysis_state.qinf,
      compressibility);

  double cpg1 = cp_first.cp;
  double cpg1_msq = cp_first.cp_msq;
  double cpg1_alf = cp_first.cp_velocity_derivative;

  for (int i = 0; i < point_count; i++) {
    const int ip = (i + 1) % point_count;
    const auto cp_next = computePressureCoefficient(
        state.surfaceVortex(0, ip), state.surfaceVortex(1, ip),
        analysis_state.qinf, compressibility);

    const double cpg2 = cp_next.cp;
    const double cpg2_msq = cp_next.cp_msq;
    const double cpg2_alf = cp_next.cp_velocity_derivative;

    const Eigen::Vector2d delta =
        foil.foil_shape.points.col(ip) - foil.foil_shape.points.col(i);
    const Eigen::Vector2d dpoint = rotateMatrix * delta;
    const double dg = cpg2 - cpg1;

    const Eigen::Vector2d apoint =
        rotateMatrix *
        ((foil.foil_shape.points.col(ip) + foil.foil_shape.points.col(i)) / 2 -
         ref);
    const double ag = 0.5 * (cpg2 + cpg1);

    const double dx_alf =
        MathUtil::cross2(delta, rotateMatrix.row(0).transpose());
    const double ag_alf = 0.5 * (cpg2_alf + cpg1_alf);
    const double ag_msq = 0.5 * (cpg2_msq + cpg1_msq);

    result.cl += dpoint.x() * ag;
    result.cm -= dpoint.dot(ag * apoint + dg * dpoint / 12.0);

    xcp_accumulator +=
        dpoint.x() * ag *
        (foil.foil_shape.points.col(ip).x() + foil.foil_shape.points.col(i).x()) /
        2.0;

    result.cl_alf += dpoint.x() * ag_alf + ag * dx_alf;
    result.cl_msq += dpoint.x() * ag_msq;

    cpg1 = cpg2;
    cpg1_alf = cpg2_alf;
    cpg1_msq = cpg2_msq;
  }

  result.xcp = std::fabs(result.cl) > 0.0 ? xcp_accumulator / result.cl : 0.0;
  return result;
}

double computeCd(const FlowState &analysis_state,
                 const BoundaryLayer &boundary_layer, bool bl_initialized) {
  if (!(analysis_state.viscous && bl_initialized)) {
    return 0.0;
  }

  const double beta =
      std::sqrt(std::max(0.0, 1.0 - analysis_state.currentMach * analysis_state.currentMach));
  const double tklam_local =
      MathUtil::pow(analysis_state.currentMach / (1.0 + beta), 2);

  const auto bottom = boundary_layer.readSideModel(2);
  const double thwake = bottom.lastMomentumThickness;
  const double edgeVelocityBottom = bottom.lastEdgeVelocity;
  const double urat = edgeVelocityBottom / analysis_state.qinf;
  const double uewake = edgeVelocityBottom * (1.0 - tklam_local) /
                        (1.0 - tklam_local * urat * urat);
  const double shwake =
      bottom.lastDisplacementThickness / bottom.lastMomentumThickness;

  const double exponent = 0.5 * (5.0 + shwake);
  const double wake_ratio = uewake / analysis_state.qinf;
  const double wake_term = std::pow(wake_ratio, exponent);
  return 2.0 * thwake * wake_term;
}

void applyClComputation(AeroCoefficients &aero_coeffs,
                        const XFoilClComputation &result) {
  aero_coeffs.cl = result.cl;
  aero_coeffs.cm = result.cm;
  aero_coeffs.cl_alf = result.cl_alf;
  aero_coeffs.cl_msq = result.cl_msq;
  aero_coeffs.xcp = result.xcp;
}

} // namespace xfoil_postprocess
