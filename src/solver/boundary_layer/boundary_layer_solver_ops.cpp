#include "solver/boundary_layer/boundary_layer_solver_ops.hpp"

#include <cmath>

#include "numerics/boundary_layer_util.hpp"
#include "numerics/coefficient/skin_friction.hpp"

bool BoundaryLayerSolverOps::blkin(BoundaryLayerState &state) const {
  const auto &compressibility = context_.blCompressibility;
  const auto &reynolds = context_.blReynolds;
  blData &current = state.current();
  current.param.mz = current.param.uz * current.param.uz *
                     compressibility.hstinv /
                     (compressibility.gm1bl *
                      (1.0 - 0.5 * current.param.uz * current.param.uz *
                                 compressibility.hstinv));
  double tr2 = 1.0 + 0.5 * compressibility.gm1bl * current.param.mz;
  current.param.mz_uz = 2.0 * current.param.mz * tr2 / current.param.uz;
  current.param.mz_ms =
      current.param.uz * current.param.uz * tr2 /
      (compressibility.gm1bl *
       (1.0 - 0.5 * current.param.uz * current.param.uz *
                  compressibility.hstinv)) *
      compressibility.hstinv_ms;

  current.param.rz = compressibility.rstbl *
                     pow(tr2, (-1.0 / compressibility.gm1bl));
  current.param.rz_uz = -current.param.rz / tr2 * 0.5 * current.param.mz_uz;
  current.param.rz_ms =
      -current.param.rz / tr2 * 0.5 * current.param.mz_ms +
      compressibility.rstbl_ms * pow(tr2, (-1.0 / compressibility.gm1bl));

  current.param.hz = current.param.dz / current.param.tz;
  current.param.hz_dz = 1.0 / current.param.tz;
  current.param.hz_tz = -current.param.hz / current.param.tz;

  double herat = 1.0 - 0.5 * current.param.uz * current.param.uz *
                           compressibility.hstinv;
  double he_u2 = -current.param.uz * compressibility.hstinv;
  double he_ms = -0.5 * current.param.uz * current.param.uz *
                 compressibility.hstinv_ms;
  double v2_he = (1.5 / herat -
                  1.0 / (herat + BoundaryLayerWorkflow::kHvrat));

  boundary_layer::KineticShapeParameterResult hkin_result =
      boundary_layer::hkin(current.param.hz, current.param.mz);
  current.hkz.scalar = hkin_result.hk;

  current.hkz.u() = hkin_result.hk_msq * current.param.mz_uz;
  current.hkz.t() = hkin_result.hk_h * current.param.hz_tz;
  current.hkz.d() = hkin_result.hk_h * current.param.hz_dz;
  current.hkz.ms() = hkin_result.hk_msq * current.param.mz_ms;

  current.rtz.scalar =
      current.param.rz * current.param.uz * current.param.tz /
      (sqrt(herat * herat * herat) * (1.0 + BoundaryLayerWorkflow::kHvrat) /
       (herat + BoundaryLayerWorkflow::kHvrat) / reynolds.reybl);
  current.rtz.u() =
      current.rtz.scalar *
      (1.0 / current.param.uz + current.param.rz_uz / current.param.rz -
       v2_he * he_u2);
  current.rtz.t() = current.rtz.scalar / current.param.tz;
  current.rtz.ms() =
      current.rtz.scalar *
      (current.param.rz_ms / current.param.rz +
       (1 / reynolds.reybl * reynolds.reybl_ms - v2_he * he_ms));
  current.rtz.re() = current.rtz.scalar * (reynolds.reybl_re / reynolds.reybl);

  return true;
}

SkinFrictionCoefficients BoundaryLayerSolverOps::blmid(
    FlowRegimeEnum flowRegimeType) const {
  blData &previous = context_.state.previous();
  blData &current = context_.state.current();

  if (flowRegimeType == FlowRegimeEnum::Similarity) {
    previous.hkz = current.hkz;
    previous.rtz = current.rtz;
    previous.param.mz = current.param.mz;
    previous.param.mz_uz = current.param.mz_uz;
    previous.param.mz_ms = current.param.mz_ms;
  }

  const double hka = 0.5 * (previous.hkz.scalar + current.hkz.scalar);
  const double rta = 0.5 * (previous.rtz.scalar + current.rtz.scalar);
  const double ma = 0.5 * (previous.param.mz + current.param.mz);

  skin_friction::C_f cf_res =
      skin_friction::getSkinFriction(hka, rta, ma, flowRegimeType);

  SkinFrictionCoefficients coeffs;
  coeffs.cfm = cf_res.cf;
  const double cfm_hka = cf_res.hk;
  const double cfm_rta = cf_res.rt;
  const double cfm_ma = cf_res.msq;

  coeffs.cfm_u1 =
      0.5 * (cfm_hka * previous.hkz.u() + cfm_ma * previous.param.mz_uz +
             cfm_rta * previous.rtz.u());
  coeffs.cfm_t1 =
      0.5 * (cfm_hka * previous.hkz.t() + cfm_rta * previous.rtz.t());
  coeffs.cfm_d1 = 0.5 * (cfm_hka * previous.hkz.d());

  coeffs.cfm_u2 =
      0.5 * (cfm_hka * current.hkz.u() + cfm_ma * current.param.mz_uz +
             cfm_rta * current.rtz.u());
  coeffs.cfm_t2 =
      0.5 * (cfm_hka * current.hkz.t() + cfm_rta * current.rtz.t());
  coeffs.cfm_d2 = 0.5 * (cfm_hka * current.hkz.d());

  coeffs.cfm_ms =
      0.5 * (cfm_hka * previous.hkz.ms() + cfm_ma * previous.param.mz_ms +
             cfm_rta * previous.rtz.ms() + cfm_hka * current.hkz.ms() +
             cfm_ma * current.param.mz_ms + cfm_rta * current.rtz.ms());
  coeffs.cfm_re =
      0.5 * (cfm_rta * previous.rtz.re() + cfm_rta * current.rtz.re());

  return coeffs;
}

blData BoundaryLayerSolverOps::blprv(blData data, double xsi, double ami,
                                     double cti, double thi, double dsi,
                                     double dswaki, double uei) const {
  const auto &compressibility = context_.blCompressibility;
  data.param.xz = xsi;
  data.param.amplz = ami;
  data.param.sz = cti;
  data.param.tz = thi;
  data.param.dz = dsi - dswaki;
  data.param.dwz = dswaki;

  data.param.uz =
      uei * (1.0 - compressibility.tkbl) /
      (1.0 - compressibility.tkbl * (uei / compressibility.qinfbl) *
                 (uei / compressibility.qinfbl));
  data.param.uz_uei =
      (1.0 + compressibility.tkbl *
                 (2.0 * data.param.uz * uei / compressibility.qinfbl /
                      compressibility.qinfbl -
                  1.0)) /
      (1.0 - compressibility.tkbl * (uei / compressibility.qinfbl) *
                 (uei / compressibility.qinfbl));
  data.param.uz_ms =
      (data.param.uz * (uei / compressibility.qinfbl) *
           (uei / compressibility.qinfbl) -
       uei) *
      compressibility.tkbl_ms /
      (1.0 - compressibility.tkbl * (uei / compressibility.qinfbl) *
                 (uei / compressibility.qinfbl));
  return data;
}

bool BoundaryLayerSolverOps::blsys() const {
  blData &previous = context_.state.previous();
  blData &current = context_.state.current();

  SkinFrictionCoefficients skinFriction = blmid(context_.flowRegime);
  current = context_.boundaryLayerVariablesSolver.solve(current,
                                                        context_.flowRegime);

  if (context_.flowRegime == FlowRegimeEnum::Similarity) {
    context_.state.stepbl();
  }

  if (context_.flowRegime == FlowRegimeEnum::Transition) {
    context_.transitionSolver.trdif();
  } else {
    context_.blc = context_.blDiffSolver.solve(context_.flowRegime, context_.state,
                                               skinFriction,
                                               context_.blTransition.amcrit);
  }

  if (context_.flowRegime == FlowRegimeEnum::Similarity) {
    context_.blc.a2 += context_.blc.a1;
    context_.blc.a1.setZero();
  }

  for (int k = 0; k < 4; ++k) {
    double res_u1 = context_.blc.a1(k, 3);
    double res_u2 = context_.blc.a2(k, 3);
    double res_ms = context_.blc.d_msq[k];

    context_.blc.a1(k, 3) *= previous.param.uz_uei;
    context_.blc.a2(k, 3) *= current.param.uz_uei;
    context_.blc.d_msq[k] =
        res_u1 * previous.param.uz_ms + res_u2 * current.param.uz_ms + res_ms;
  }

  return true;
}

bool BoundaryLayerSolverOps::tesys(
    const BoundaryLayerSideProfiles &top_profiles,
    const BoundaryLayerSideProfiles &bottom_profiles, const Edge &edge) const {
  context_.blc.clear();

  context_.state.station2 = context_.boundaryLayerVariablesSolver.solve(
      context_.state.station2, FlowRegimeEnum::Wake);

  const int top_te = context_.lattice.top.trailingEdgeIndex;
  const int bottom_te = context_.lattice.bottom.trailingEdgeIndex;
  const double tte = top_profiles.momentumThickness[top_te] +
                     bottom_profiles.momentumThickness[bottom_te];
  const double dte = top_profiles.displacementThickness[top_te] +
                     bottom_profiles.displacementThickness[bottom_te] +
                     edge.ante;
  const double cte =
      (top_profiles.skinFrictionCoeff[top_te] *
           top_profiles.momentumThickness[top_te] +
       bottom_profiles.skinFrictionCoeff[bottom_te] *
           bottom_profiles.momentumThickness[bottom_te]) /
      tte;

  context_.blc.a1(0, 0) = -1.0;
  context_.blc.a2(0, 0) = 1.0;
  context_.blc.rhs[0] = cte - context_.state.station2.param.sz;

  context_.blc.a1(1, 1) = -1.0;
  context_.blc.a2(1, 1) = 1.0;
  context_.blc.rhs[1] = tte - context_.state.station2.param.tz;

  context_.blc.a1(2, 2) = -1.0;
  context_.blc.a2(2, 2) = 1.0;
  context_.blc.rhs[2] =
      dte - context_.state.station2.param.dz - context_.state.station2.param.dwz;

  return true;
}
