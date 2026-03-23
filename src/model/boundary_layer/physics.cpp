#include "model/boundary_layer/physics.hpp"

#include <algorithm>
#include <cmath>

#include "numerics/coefficient/skin_friction.hpp"
#include "numerics/math_util.hpp"

bool BoundaryLayerPhysics::blkin(BoundaryLayerStationWindow &state, const BlCompressibilityParams &compressibility,
                                 const BlReynoldsParams &reynolds) {
    BoundaryLayerStationState &current = state.current();
    current.param.mz =
        current.param.uz * current.param.uz * compressibility.hstinv() /
        (compressibility.gm1bl() * (1.0 - 0.5 * current.param.uz * current.param.uz * compressibility.hstinv()));
    const double tr2    = 1.0 + 0.5 * compressibility.gm1bl() * current.param.mz;
    current.param.mz_uz = 2.0 * current.param.mz * tr2 / current.param.uz;
    current.param.mz_ms =
        current.param.uz * current.param.uz * tr2 /
        (compressibility.gm1bl() * (1.0 - 0.5 * current.param.uz * current.param.uz * compressibility.hstinv())) *
        compressibility.hstinv_ms();

    current.param.rz    = compressibility.rstbl() * std::pow(tr2, (-1.0 / compressibility.gm1bl()));
    current.param.rz_uz = -current.param.rz / tr2 * 0.5 * current.param.mz_uz;
    current.param.rz_ms = -current.param.rz / tr2 * 0.5 * current.param.mz_ms +
                          compressibility.rstbl_ms() * std::pow(tr2, (-1.0 / compressibility.gm1bl()));

    current.param.hz    = current.param.dz / current.param.tz;
    current.param.hz_dz = 1.0 / current.param.tz;
    current.param.hz_tz = -current.param.hz / current.param.tz;

    const double herat = 1.0 - 0.5 * current.param.uz * current.param.uz * compressibility.hstinv();
    const double he_u2 = -current.param.uz * compressibility.hstinv();
    const double he_ms = -0.5 * current.param.uz * current.param.uz * compressibility.hstinv_ms();
    const double v2_he = (1.5 / herat - 1.0 / (herat + kHvrat));

    const auto hkin_result = boundary_layer::hkin(current.param.hz, current.param.mz);
    current.hkz.scalar     = hkin_result.hk;
    current.hkz.u()        = hkin_result.hk_msq * current.param.mz_uz;
    current.hkz.t()        = hkin_result.hk_h * current.param.hz_tz;
    current.hkz.d()        = hkin_result.hk_h * current.param.hz_dz;
    current.hkz.ms()       = hkin_result.hk_msq * current.param.mz_ms;

    current.rtz.scalar = current.param.rz * current.param.uz * current.param.tz /
                         (std::sqrt(herat * herat * herat) * (1.0 + kHvrat) / (herat + kHvrat) / reynolds.reybl);
    current.rtz.u() =
        current.rtz.scalar * (1.0 / current.param.uz + current.param.rz_uz / current.param.rz - v2_he * he_u2);
    current.rtz.t()  = current.rtz.scalar / current.param.tz;
    current.rtz.ms() = current.rtz.scalar * (current.param.rz_ms / current.param.rz +
                                             (1.0 / reynolds.reybl * reynolds.reybl_ms - v2_he * he_ms));
    current.rtz.re() = current.rtz.scalar * (reynolds.reybl_re / reynolds.reybl);

    return true;
}

void BoundaryLayerPhysics::refreshCurrentStation(BoundaryLayerStationWindow &state,
                                                 const BlCompressibilityParams &compressibility,
                                                 const BlReynoldsParams &reynolds, double xsi, double ami, double cti,
                                                 double thi, double dsi, double dswaki, double uei) {
    state.current() = blprv(state.current(), compressibility, xsi, ami, cti, thi, dsi, dswaki, uei);
    blkin(state, compressibility, reynolds);
}

SkinFrictionCoefficients BoundaryLayerPhysics::blmid(BoundaryLayerStationWindow &state, FlowRegimeEnum flowRegimeType) {
    BoundaryLayerStationState &previous = state.previous();
    BoundaryLayerStationState &current  = state.current();

    if (flowRegimeType == FlowRegimeEnum::Similarity) {
        previous.hkz         = current.hkz;
        previous.rtz         = current.rtz;
        previous.param.mz    = current.param.mz;
        previous.param.mz_uz = current.param.mz_uz;
        previous.param.mz_ms = current.param.mz_ms;
    }

    const double hka = 0.5 * (previous.hkz.scalar + current.hkz.scalar);
    const double rta = 0.5 * (previous.rtz.scalar + current.rtz.scalar);
    const double ma  = 0.5 * (previous.param.mz + current.param.mz);

    const auto cf_res = skin_friction::getSkinFriction(hka, rta, ma, flowRegimeType);

    SkinFrictionCoefficients coeffs;
    coeffs.cfm           = cf_res.cf;
    const double cfm_hka = cf_res.hk;
    const double cfm_rta = cf_res.rt;
    const double cfm_ma  = cf_res.msq;

    coeffs.cfm_u1 = 0.5 * (cfm_hka * previous.hkz.u() + cfm_ma * previous.param.mz_uz + cfm_rta * previous.rtz.u());
    coeffs.cfm_t1 = 0.5 * (cfm_hka * previous.hkz.t() + cfm_rta * previous.rtz.t());
    coeffs.cfm_d1 = 0.5 * (cfm_hka * previous.hkz.d());

    coeffs.cfm_u2 = 0.5 * (cfm_hka * current.hkz.u() + cfm_ma * current.param.mz_uz + cfm_rta * current.rtz.u());
    coeffs.cfm_t2 = 0.5 * (cfm_hka * current.hkz.t() + cfm_rta * current.rtz.t());
    coeffs.cfm_d2 = 0.5 * (cfm_hka * current.hkz.d());

    coeffs.cfm_ms = 0.5 * (cfm_hka * previous.hkz.ms() + cfm_ma * previous.param.mz_ms + cfm_rta * previous.rtz.ms() +
                           cfm_hka * current.hkz.ms() + cfm_ma * current.param.mz_ms + cfm_rta * current.rtz.ms());
    coeffs.cfm_re = 0.5 * (cfm_rta * previous.rtz.re() + cfm_rta * current.rtz.re());

    return coeffs;
}

BoundaryLayerStationState BoundaryLayerPhysics::blprv(BoundaryLayerStationState data,
                                                      const BlCompressibilityParams &compressibility, double xsi,
                                                      double ami, double cti, double thi, double dsi, double dswaki,
                                                      double uei) {
    data.assignPrimaryStationData(xsi, ami, cti, thi, dsi, dswaki);
    data.assignCompressibleEdgeVelocity(uei, compressibility.qinfbl(), compressibility.tkbl(), compressibility.tkbl_ms());
    return data;
}

double BoundaryLayerPhysics::adjustDisplacementForHkLimit(double displacementThickness, double momentumThickness,
                                                          double msq, double hklim) {
    const double h         = displacementThickness / momentumThickness;
    const auto hkin_result = boundary_layer::hkin(h, msq);
    const double dh        = std::max(0.0, hklim - hkin_result.hk) / hkin_result.hk_h;
    return displacementThickness + dh * momentumThickness;
}

BoundaryLayerReferenceParams BoundaryLayerPhysics::buildReferenceParams(const FlowState &analysisState,
                                                                        const AeroCoefficients &aeroCoefficients,
                                                                        double acrit) {
    BoundaryLayerReferenceParams params;
    const double clmr = analysisState.controlByAlpha ? aeroCoefficients.cl : analysisState.clspec;
    const double cla  = std::max(clmr, 0.000001);

    switch (analysisState.machType) {
    case FlowState::MachType::CONSTANT:
        params.currentMach = analysisState.referenceMach;
        break;
    case FlowState::MachType::FIXED_LIFT:
        params.currentMach = analysisState.referenceMach / std::sqrt(cla);
        break;
    case FlowState::MachType::FIXED_LIFT_AND_DYNAMIC_PRESSURE:
        params.currentMach = analysisState.referenceMach;
        break;
    default:
        params.currentMach = analysisState.currentMach;
        break;
    }

    double ma_clmr = 0.0;
    if (analysisState.machType == FlowState::MachType::FIXED_LIFT) {
        ma_clmr = -0.5 * params.currentMach / cla;
    }

    switch (analysisState.reynoldsType) {
    case FlowState::ReynoldsType::CONSTANT:
        params.currentRe = analysisState.referenceRe;
        break;
    case FlowState::ReynoldsType::FIXED_LIFT:
        params.currentRe = analysisState.referenceRe / std::sqrt(cla);
        params.re_clmr   = -0.5 * params.currentRe / cla;
        break;
    case FlowState::ReynoldsType::FIXED_LIFT_AND_DYNAMIC_PRESSURE:
        params.currentRe = analysisState.referenceRe / cla;
        params.re_clmr   = -params.currentRe / cla;
        break;
    default:
        params.currentRe = analysisState.currentRe;
        break;
    }

    params.msq_clmr = 2.0 * params.currentMach * ma_clmr;

    params.blCompressibility = BlCompressibilityParams(analysisState.qinf, params.currentMach);

    const double herat =
        1.0 - 0.5 * params.blCompressibility.qinfbl() * params.blCompressibility.qinfbl() *
                    params.blCompressibility.hstinv();
    const double herat_ms =
        -0.5 * params.blCompressibility.qinfbl() * params.blCompressibility.qinfbl() *
            params.blCompressibility.hstinv_ms();

    params.blReynolds.reybl = params.currentRe * std::sqrt(herat * herat * herat) * (1.0 + kHvrat) / (herat + kHvrat);
    params.blReynolds.reybl_re = std::sqrt(herat * herat * herat) * (1.0 + kHvrat) / (herat + kHvrat);
    params.blReynolds.reybl_ms = params.blReynolds.reybl * (1.5 / herat - 1.0 / (herat + kHvrat)) * herat_ms;

    params.amcrit = acrit;
    return params;
}
