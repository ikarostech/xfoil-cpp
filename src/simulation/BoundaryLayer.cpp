#include "BoundaryLayer.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

#include "XFoil.h"
#include "domain/coefficient/skin_friction.hpp"

using BoundaryContext = BoundaryLayerWorkflow::MixedModeStationContext;
using Eigen::Matrix;
using Eigen::Vector;

bool BoundaryLayerWorkflow::isStartOfWake(const XFoil& xfoil, int side,
                                          int stationIndex) {
  return stationIndex == lattice.get(side).trailingEdgeIndex + 1;
}

void BoundaryLayerWorkflow::updateSystemMatricesForStation(
    XFoil& xfoil, int side, int stationIndex, BoundaryContext& ctx) {
  if (isStartOfWake(xfoil, side, stationIndex)) {
    ctx.tte = lattice.get(1).profiles.momentumThickness[lattice.top.trailingEdgeIndex] +
              lattice.get(2).profiles.momentumThickness[lattice.bottom.trailingEdgeIndex];
    ctx.dte = lattice.get(1).profiles.displacementThickness[lattice.top.trailingEdgeIndex] +
              lattice.get(2).profiles.displacementThickness[lattice.bottom.trailingEdgeIndex] + xfoil.foil.edge.ante;
    ctx.cte =
        (lattice.get(1).profiles.skinFrictionCoeff[lattice.top.trailingEdgeIndex] *
             lattice.get(1).profiles.momentumThickness[lattice.top.trailingEdgeIndex] +
         lattice.get(2).profiles.skinFrictionCoeff[lattice.bottom.trailingEdgeIndex] *
             lattice.get(2).profiles.momentumThickness[lattice.bottom.trailingEdgeIndex]) /
        ctx.tte;
    tesys(xfoil, ctx.cte, ctx.tte, ctx.dte);
  } else {
    blsys(xfoil);
  }
}

void BoundaryLayerWorkflow::initializeFirstIterationState(
    XFoil& xfoil, int side, int stationIndex, int previousTransition,
    BoundaryContext& ctx, double& ueref, double& hkref, double& ami) {
  ueref = state.station2.param.uz;
  hkref = state.station2.hkz.scalar;

  const bool inLaminarWindow =
      stationIndex < lattice.get(side).transitionIndex && stationIndex >= previousTransition;
  if (inLaminarWindow) {
    double uem;
    double dsm;
    double thm;
    if (stationIndex > 0) {
      uem = lattice.get(side).profiles.edgeVelocity[stationIndex - 1];
      dsm = lattice.get(side).profiles.displacementThickness[stationIndex - 1];
      thm = lattice.get(side).profiles.momentumThickness[stationIndex - 1];
    } else {
      uem = lattice.get(side).profiles.edgeVelocity[stationIndex];
      dsm = lattice.get(side).profiles.displacementThickness[stationIndex];
      thm = lattice.get(side).profiles.momentumThickness[stationIndex];
    }
    const double uem_sq = uem * uem;
    const double msq =
        uem_sq * xfoil.hstinv /
        (xfoil.gm1bl * (1.0 - 0.5 * uem_sq * xfoil.hstinv));
    const auto hkin_result =
        boundary_layer::hkin(dsm / thm, msq);
    hkref = hkin_result.hk;
  }

  if (stationIndex < previousTransition) {
    if (xfoil.flowRegime == FlowRegimeEnum::Transition) {
      lattice.get(side).profiles.skinFrictionCoeff[stationIndex] = 0.03;
    }
    if (xfoil.flowRegime == FlowRegimeEnum::Turbulent || xfoil.flowRegime == FlowRegimeEnum::Wake) {
      const double prev =
          (stationIndex >= 1) ? lattice.get(side).profiles.skinFrictionCoeff[stationIndex - 1]
                              : lattice.get(side).profiles.skinFrictionCoeff[stationIndex];
      lattice.get(side).profiles.skinFrictionCoeff[stationIndex] = prev;
    }
    if (xfoil.flowRegime == FlowRegimeEnum::Transition || xfoil.flowRegime == FlowRegimeEnum::Turbulent || xfoil.flowRegime == FlowRegimeEnum::Wake) {
      ctx.cti = lattice.get(side).profiles.skinFrictionCoeff[stationIndex - 1];
      state.station2.param.sz = ctx.cti;
    }
  }
}

void BoundaryLayerWorkflow::configureSimilarityRow(double ueref) {
  blc.a2(3, 0) = 0.0;
  blc.a2(3, 1) = 0.0;
  blc.a2(3, 2) = 0.0;
  blc.a2(3, 3) = state.station2.param.uz_uei;
  blc.rhs[3] = ueref - state.station2.param.uz;
}

void BoundaryLayerWorkflow::configureViscousRow(double hkref,
                                                double ueref, double senswt,
                                                bool resetSensitivity,
                                                bool averageSensitivity,
                                                double& sens, double& sennew) {
  blc.a2(3, 0) = 0.0;
  blc.a2(3, 1) = state.station2.hkz.t();
  blc.a2(3, 2) = state.station2.hkz.d();
  blc.a2(3, 3) = state.station2.hkz.u() * state.station2.param.uz_uei;
  blc.rhs[3] = 1.0;

  const double delta_sen =
      blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs)[3];

  sennew = senswt * delta_sen * hkref / ueref;
  if (resetSensitivity) {
    sens = sennew;
  } else if (averageSensitivity) {
    sens = 0.5 * (sens + sennew);
  }

  blc.a2(3, 1) = state.station2.hkz.t() * hkref;
  blc.a2(3, 2) = state.station2.hkz.d() * hkref;
  blc.a2(3, 3) =
      (state.station2.hkz.u() * hkref + sens / ueref) * state.station2.param.uz_uei;
  blc.rhs[3] =
      -(hkref * hkref) * (state.station2.hkz.scalar / hkref - 1.0) -
      sens * (state.station2.param.uz / ueref - 1.0);
}

bool BoundaryLayerWorkflow::applyMixedModeNewtonStep(
    XFoil& xfoil, int side, int stationIndex, double deps, double& ami,
    BoundaryContext& ctx) {
  blc.rhs =
      blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);

  ctx.dmax = std::max(std::fabs(blc.rhs[1] / ctx.thi),
                      std::fabs(blc.rhs[2] / ctx.dsi));
  if (stationIndex >= lattice.get(side).transitionIndex) {
    ctx.dmax = std::max(ctx.dmax,
                        std::fabs(blc.rhs[0] / (10.0 * ctx.cti)));
  }

  xfoil.rlx = 1.0;
  if (ctx.dmax > 0.3) {
    xfoil.rlx = 0.3 / ctx.dmax;
  }

  if (stationIndex < lattice.get(side).transitionIndex) {
    ami += xfoil.rlx * blc.rhs[0];
    ctx.ami = ami;
  }
  if (stationIndex >= lattice.get(side).transitionIndex) {
    ctx.cti += xfoil.rlx * blc.rhs[0];
  }
  ctx.thi += xfoil.rlx * blc.rhs[1];
  ctx.dsi += xfoil.rlx * blc.rhs[2];
  ctx.uei += xfoil.rlx * blc.rhs[3];

  if (stationIndex >= lattice.get(side).transitionIndex) {
    ctx.cti = std::clamp(ctx.cti, 0.0000001, 0.30);
  }

  const double hklim =
      (stationIndex <= lattice.get(side).trailingEdgeIndex) ? 1.02 : 1.00005;
  const double uei_sq = ctx.uei * ctx.uei;
  const double msq = uei_sq * xfoil.hstinv /
                     (xfoil.gm1bl * (1.0 - 0.5 * uei_sq * xfoil.hstinv));
  double dsw = ctx.dsi - ctx.dswaki;
  xfoil.dslim(dsw, ctx.thi, msq, hklim);
  ctx.dsi = dsw + ctx.dswaki;

  return ctx.dmax <= deps;
}

blData BoundaryLayerWorkflow::blvar(blData data, FlowRegimeEnum flowRegimeType) {
  return this->boundaryLayerVariablesSolver.solve(data, flowRegimeType);
}

SkinFrictionCoefficients BoundaryLayerWorkflow::blmid(
    XFoil& xfoil, FlowRegimeEnum flowRegimeType) {
  BoundaryLayerState& state = this->state;
  blData& previous = state.previous();
  blData& current = state.current();

  if (xfoil.flowRegime == FlowRegimeEnum::Similarity) {
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

  coeffs.cfm_u1 = 0.5 * (cfm_hka * previous.hkz.u() +
                         cfm_ma * previous.param.mz_uz +
                         cfm_rta * previous.rtz.u());
  coeffs.cfm_t1 = 0.5 * (cfm_hka * previous.hkz.t() +
                         cfm_rta * previous.rtz.t());
  coeffs.cfm_d1 = 0.5 * (cfm_hka * previous.hkz.d());

  coeffs.cfm_u2 = 0.5 * (cfm_hka * current.hkz.u() +
                         cfm_ma * current.param.mz_uz +
                         cfm_rta * current.rtz.u());
  coeffs.cfm_t2 = 0.5 * (cfm_hka * current.hkz.t() +
                         cfm_rta * current.rtz.t());
  coeffs.cfm_d2 = 0.5 * (cfm_hka * current.hkz.d());

  coeffs.cfm_ms =
      0.5 * (cfm_hka * previous.hkz.ms() + cfm_ma * previous.param.mz_ms +
             cfm_rta * previous.rtz.ms() + cfm_hka * current.hkz.ms() +
             cfm_ma * current.param.mz_ms + cfm_rta * current.rtz.ms());
  coeffs.cfm_re = 0.5 * (cfm_rta * previous.rtz.re() +
                         cfm_rta * current.rtz.re());

  return coeffs;
}

blData BoundaryLayerWorkflow::blprv(XFoil& xfoil, blData data, double xsi,
                                    double ami, double cti, double thi,
                                    double dsi, double dswaki,
                                    double uei) const {
  data.param.xz = xsi;
  data.param.amplz = ami;
  data.param.sz = cti;
  data.param.tz = thi;
  data.param.dz = dsi - dswaki;
  data.param.dwz = dswaki;

  data.param.uz =
      uei * (1.0 - xfoil.tkbl) /
      (1.0 - xfoil.tkbl * (uei / xfoil.qinfbl) * (uei / xfoil.qinfbl));
  data.param.uz_uei =
      (1.0 + xfoil.tkbl *
                (2.0 * data.param.uz * uei / xfoil.qinfbl / xfoil.qinfbl -
                 1.0)) /
      (1.0 - xfoil.tkbl * (uei / xfoil.qinfbl) * (uei / xfoil.qinfbl));
  data.param.uz_ms =
      (data.param.uz * (uei / xfoil.qinfbl) * (uei / xfoil.qinfbl) - uei) *
      xfoil.tkbl_ms /
      (1.0 - xfoil.tkbl * (uei / xfoil.qinfbl) * (uei / xfoil.qinfbl));
  return data;
}

bool BoundaryLayerWorkflow::blsys(XFoil& xfoil) {
  blData& previous = state.previous();
  blData& current = state.current();

  SkinFrictionCoefficients skinFriction;

  if (xfoil.flowRegime == FlowRegimeEnum::Wake) {
    current = blvar(current, FlowRegimeEnum::Wake);
    skinFriction = blmid(xfoil, FlowRegimeEnum::Wake);
  } else if (xfoil.flowRegime == FlowRegimeEnum::Transition || xfoil.flowRegime == FlowRegimeEnum::Turbulent || xfoil.flowRegime == FlowRegimeEnum::Wake) {
    current = blvar(current, FlowRegimeEnum::Turbulent);
    skinFriction = blmid(xfoil, FlowRegimeEnum::Turbulent);
  } else {
    current = blvar(current, FlowRegimeEnum::Laminar);
    skinFriction = blmid(xfoil, FlowRegimeEnum::Laminar);
  }

  if (xfoil.flowRegime == FlowRegimeEnum::Similarity) {
    state.stepbl();
  }

  if (xfoil.flowRegime == FlowRegimeEnum::Transition) {
    trdif(xfoil);
  } else if (xfoil.flowRegime == FlowRegimeEnum::Similarity) {
    blc = xfoil.blDiffSolver.solve(FlowRegimeEnum::Similarity, state,
                                   skinFriction, xfoil.amcrit);
  } else if (!(xfoil.flowRegime == FlowRegimeEnum::Turbulent || xfoil.flowRegime == FlowRegimeEnum::Wake)) {
    blc = xfoil.blDiffSolver.solve(FlowRegimeEnum::Laminar, state,
                                   skinFriction, xfoil.amcrit);
  } else if (xfoil.flowRegime == FlowRegimeEnum::Wake) {
    blc = xfoil.blDiffSolver.solve(FlowRegimeEnum::Wake, state,
                                   skinFriction, xfoil.amcrit);
  } else {
    blc = xfoil.blDiffSolver.solve(FlowRegimeEnum::Turbulent, state,
                                   skinFriction, xfoil.amcrit);
  }

  if (xfoil.flowRegime == FlowRegimeEnum::Similarity) {
    blc.a2 += blc.a1;
    blc.a1.setZero();
  }

  for (int k = 0; k < 4; ++k) {
    double res_u1 = blc.a1(k, 3);
    double res_u2 = blc.a2(k, 3);
    double res_ms = blc.d_msq[k];

    blc.a1(k, 3) *= previous.param.uz_uei;
    blc.a2(k, 3) *= current.param.uz_uei;
    blc.d_msq[k] =
        res_u1 * previous.param.uz_ms + res_u2 * current.param.uz_ms + res_ms;
  }

  return true;
}

bool BoundaryLayerWorkflow::trdif(XFoil& xfoil) {
  //-----------------------------------------------
  //     sets up the newton system governing the
  //     transition interval.  equations governing
  //     the  laminar  part  x1 < xi < xt  and
  //     the turbulent part  xt < xi < x2
  //     are simply summed.
  //-----------------------------------------------
  Matrix<double, 4, 5> bl1, bl2, bt1, bt2;
  Vector<double, 4> blrez, blm, blr, blx, btrez, btm, btr, btx;

  double tt, tt_a1, tt_x1, tt_x2, tt_t1, tt_t2, tt_d1, tt_d2, tt_u1, tt_u2;
  double tt_ms, tt_re, tt_xf, dt, dt_a1, dt_x1, dt_x2, dt_t1, dt_t2;
  double dt_d1, dt_d2, dt_u1, dt_u2, dt_ms, dt_re, dt_xf;
  double ut, ut_a1, ut_x1, ut_x2, ut_t1, ut_t2, ut_d1, ut_d2, ut_u1, ut_u2;
  double ut_ms, ut_re, ut_xf;
  double st, st_tt, st_dt, st_ut, st_ms, st_re, st_a1, st_x1, st_x2, st_t1,
      st_t2;
  double st_d1, st_d2, st_u1, st_u2, st_xf;
  double ctr, ctr_hk2;

  xfoil.saveblData(1);
  xfoil.saveblData(2);

  //---- weighting factors for linear interpolation to transition point
  double wf2 = (xfoil.xt - state.station1.param.xz) / (state.station2.param.xz - state.station1.param.xz);
  double wf2_xt = 1.0 / (state.station2.param.xz - state.station1.param.xz);

  double wf2_a1 = wf2_xt * xfoil.xt_a1;
  double wf2_x1 =
      wf2_xt * xfoil.xt_x1 + (wf2 - 1.0) / (state.station2.param.xz - state.station1.param.xz);
  double wf2_x2 = wf2_xt * xfoil.xt_x2 - wf2 / (state.station2.param.xz - state.station1.param.xz);
  double wf2_t1 = wf2_xt * xfoil.xt_t1;
  double wf2_t2 = wf2_xt * xfoil.xt_t2;
  double wf2_d1 = wf2_xt * xfoil.xt_d1;
  double wf2_d2 = wf2_xt * xfoil.xt_d2;
  double wf2_u1 = wf2_xt * xfoil.xt_u1;
  double wf2_u2 = wf2_xt * xfoil.xt_u2;
  double wf2_ms = wf2_xt * xfoil.xt_ms;
  double wf2_re = wf2_xt * xfoil.xt_re;
  double wf2_xf = wf2_xt * xfoil.xt_xf;

  double wf1 = 1.0 - wf2;
  double wf1_a1 = -wf2_a1;
  double wf1_x1 = -wf2_x1;
  double wf1_x2 = -wf2_x2;
  double wf1_t1 = -wf2_t1;
  double wf1_t2 = -wf2_t2;
  double wf1_d1 = -wf2_d1;
  double wf1_d2 = -wf2_d2;
  double wf1_u1 = -wf2_u1;
  double wf1_u2 = -wf2_u2;
  double wf1_ms = -wf2_ms;
  double wf1_re = -wf2_re;
  double wf1_xf = -wf2_xf;

  //-----interpolate primary variables to transition point
  tt = state.station1.param.tz * wf1 + state.station2.param.tz * wf2;
  tt_a1 = state.station1.param.tz * wf1_a1 + state.station2.param.tz * wf2_a1;
  tt_x1 = state.station1.param.tz * wf1_x1 + state.station2.param.tz * wf2_x1;
  tt_x2 = state.station1.param.tz * wf1_x2 + state.station2.param.tz * wf2_x2;
  tt_t1 = state.station1.param.tz * wf1_t1 + state.station2.param.tz * wf2_t1 + wf1;
  tt_t2 = state.station1.param.tz * wf1_t2 + state.station2.param.tz * wf2_t2 + wf2;
  tt_d1 = state.station1.param.tz * wf1_d1 + state.station2.param.tz * wf2_d1;
  tt_d2 = state.station1.param.tz * wf1_d2 + state.station2.param.tz * wf2_d2;
  tt_u1 = state.station1.param.tz * wf1_u1 + state.station2.param.tz * wf2_u1;
  tt_u2 = state.station1.param.tz * wf1_u2 + state.station2.param.tz * wf2_u2;
  tt_ms = state.station1.param.tz * wf1_ms + state.station2.param.tz * wf2_ms;
  tt_re = state.station1.param.tz * wf1_re + state.station2.param.tz * wf2_re;
  tt_xf = state.station1.param.tz * wf1_xf + state.station2.param.tz * wf2_xf;

  dt = state.station1.param.dz * wf1 + state.station2.param.dz * wf2;
  dt_a1 = state.station1.param.dz * wf1_a1 + state.station2.param.dz * wf2_a1;
  dt_x1 = state.station1.param.dz * wf1_x1 + state.station2.param.dz * wf2_x1;
  dt_x2 = state.station1.param.dz * wf1_x2 + state.station2.param.dz * wf2_x2;
  dt_t1 = state.station1.param.dz * wf1_t1 + state.station2.param.dz * wf2_t1;
  dt_t2 = state.station1.param.dz * wf1_t2 + state.station2.param.dz * wf2_t2;
  dt_d1 = state.station1.param.dz * wf1_d1 + state.station2.param.dz * wf2_d1 + wf1;
  dt_d2 = state.station1.param.dz * wf1_d2 + state.station2.param.dz * wf2_d2 + wf2;
  dt_u1 = state.station1.param.dz * wf1_u1 + state.station2.param.dz * wf2_u1;
  dt_u2 = state.station1.param.dz * wf1_u2 + state.station2.param.dz * wf2_u2;
  dt_ms = state.station1.param.dz * wf1_ms + state.station2.param.dz * wf2_ms;
  dt_re = state.station1.param.dz * wf1_re + state.station2.param.dz * wf2_re;
  dt_xf = state.station1.param.dz * wf1_xf + state.station2.param.dz * wf2_xf;

  ut = state.station1.param.uz * wf1 + state.station2.param.uz * wf2;
  ut_a1 = state.station1.param.uz * wf1_a1 + state.station2.param.uz * wf2_a1;
  ut_x1 = state.station1.param.uz * wf1_x1 + state.station2.param.uz * wf2_x1;
  ut_x2 = state.station1.param.uz * wf1_x2 + state.station2.param.uz * wf2_x2;
  ut_t1 = state.station1.param.uz * wf1_t1 + state.station2.param.uz * wf2_t1;
  ut_t2 = state.station1.param.uz * wf1_t2 + state.station2.param.uz * wf2_t2;
  ut_d1 = state.station1.param.uz * wf1_d1 + state.station2.param.uz * wf2_d1;
  ut_d2 = state.station1.param.uz * wf1_d2 + state.station2.param.uz * wf2_d2;
  ut_u1 = state.station1.param.uz * wf1_u1 + state.station2.param.uz * wf2_u1 + wf1;
  ut_u2 = state.station1.param.uz * wf1_u2 + state.station2.param.uz * wf2_u2 + wf2;
  ut_ms = state.station1.param.uz * wf1_ms + state.station2.param.uz * wf2_ms;
  ut_re = state.station1.param.uz * wf1_re + state.station2.param.uz * wf2_re;
  ut_xf = state.station1.param.uz * wf1_xf + state.station2.param.uz * wf2_xf;

  //---- set primary "t" variables at xt  (really placed into "2" variables)
  state.station2.param.xz = xfoil.xt;
  state.station2.param.tz = tt;
  state.station2.param.dz = dt;
  state.station2.param.uz = ut;

  state.station2.param.amplz = xfoil.amcrit;
  state.station2.param.sz = 0.0;

  //---- calculate laminar secondary "t" variables
  xfoil.blkin(state);
  state.station2 = blvar(state.station2, FlowRegimeEnum::Laminar);

  //---- calculate x1-xt midpoint cfm value
  SkinFrictionCoefficients laminarSkinFriction =
      blmid(xfoil, FlowRegimeEnum::Laminar);

  //=    at this point, all "2" variables are really "t" variables at xt

  //---- set up newton system for dam, dth, dds, due, dxi  at  x1 and xt
  blc = xfoil.blDiffSolver.solve(FlowRegimeEnum::Laminar, state, laminarSkinFriction, xfoil.amcrit);

  //---- the current newton system is in terms of "1" and "t" variables,
  //-    so calculate its equivalent in terms of "1" and "2" variables.
  //-    in other words, convert residual sensitivities wrt "t" variables
  //-    into sensitivities wrt "1" and "2" variables.  the amplification
  //-    equation is unnecessary here, so the k=1 row is left empty.
  for (int k = 1; k < 3; k++) {
    blrez[k] = blc.rhs[k];
    blm[k] = blc.d_msq[k] + blc.a2(k, 1) * tt_ms + blc.a2(k, 2) * dt_ms +
             blc.a2(k, 3) * ut_ms + blc.a2(k, 4) * xfoil.xt_ms;
    blr[k] = blc.d_re[k] + blc.a2(k, 1) * tt_re + blc.a2(k, 2) * dt_re +
             blc.a2(k, 3) * ut_re + blc.a2(k, 4) * xfoil.xt_re;
    blx[k] = blc.d_xi[k] + blc.a2(k, 1) * tt_xf + blc.a2(k, 2) * dt_xf +
             blc.a2(k, 3) * ut_xf + blc.a2(k, 4) * xfoil.xt_xf;

    bl1(k, 0) = blc.a1(k, 0) + blc.a2(k, 1) * tt_a1 + blc.a2(k, 2) * dt_a1 +
                blc.a2(k, 3) * ut_a1 + blc.a2(k, 4) * xfoil.xt_a1;
    bl1(k, 1) = blc.a1(k, 1) + blc.a2(k, 1) * tt_t1 + blc.a2(k, 2) * dt_t1 +
                blc.a2(k, 3) * ut_t1 + blc.a2(k, 4) * xfoil.xt_t1;
    bl1(k, 2) = blc.a1(k, 2) + blc.a2(k, 1) * tt_d1 + blc.a2(k, 2) * dt_d1 +
                blc.a2(k, 3) * ut_d1 + blc.a2(k, 4) * xfoil.xt_d1;
    bl1(k, 3) = blc.a1(k, 3) + blc.a2(k, 1) * tt_u1 + blc.a2(k, 2) * dt_u1 +
                blc.a2(k, 3) * ut_u1 + blc.a2(k, 4) * xfoil.xt_u1;
    bl1(k, 4) = blc.a1(k, 4) + blc.a2(k, 1) * tt_x1 + blc.a2(k, 2) * dt_x1 +
                blc.a2(k, 3) * ut_x1 + blc.a2(k, 4) * xfoil.xt_x1;

    bl2(k, 0) = 0.0;
    bl2(k, 1) = blc.a2(k, 1) * tt_t2 + blc.a2(k, 2) * dt_t2 + blc.a2(k, 3) * ut_t2 +
                blc.a2(k, 4) * xfoil.xt_t2;
    bl2(k, 2) = blc.a2(k, 1) * tt_d2 + blc.a2(k, 2) * dt_d2 + blc.a2(k, 3) * ut_d2 +
                blc.a2(k, 4) * xfoil.xt_d2;
    bl2(k, 3) = blc.a2(k, 1) * tt_u2 + blc.a2(k, 2) * dt_u2 + blc.a2(k, 3) * ut_u2 +
                blc.a2(k, 4) * xfoil.xt_u2;
    bl2(k, 4) = blc.a2(k, 1) * tt_x2 + blc.a2(k, 2) * dt_x2 + blc.a2(k, 3) * ut_x2 +
                blc.a2(k, 4) * xfoil.xt_x2;
  }

  //**** second, set up turbulent part between xt and x2  ****

  //---- calculate equilibrium shear coefficient cqt at transition point
  state.station2 = blvar(state.station2, FlowRegimeEnum::Turbulent);

  //---- set initial shear coefficient value st at transition point
  //-    ( note that cq2, cq2_t2, etc. are really "cqt", "cqt_tt", etc.)

  ctr = 1.8 * exp(-3.3 / (state.station2.hkz.scalar - 1.0));
  ctr_hk2 = ctr * 3.3 / (state.station2.hkz.scalar - 1.0) / (state.station2.hkz.scalar - 1.0);

  st = ctr * state.station2.cqz.scalar;
  st_tt =
      ctr * state.station2.cqz.t() + state.station2.cqz.scalar * ctr_hk2 * state.station2.hkz.t();
  st_dt =
      ctr * state.station2.cqz.d() + state.station2.cqz.scalar * ctr_hk2 * state.station2.hkz.d();
  st_ut =
      ctr * state.station2.cqz.u() + state.station2.cqz.scalar * ctr_hk2 * state.station2.hkz.u();
  st_ms =
      ctr * state.station2.cqz.ms() + state.station2.cqz.scalar * ctr_hk2 * state.station2.hkz.ms();
  st_re = ctr * state.station2.cqz.re();

  //---- calculate st sensitivities wrt the actual "1" and "2" variables
  st_a1 = st_tt * tt_a1 + st_dt * dt_a1 + st_ut * ut_a1;
  st_x1 = st_tt * tt_x1 + st_dt * dt_x1 + st_ut * ut_x1;
  st_x2 = st_tt * tt_x2 + st_dt * dt_x2 + st_ut * ut_x2;
  st_t1 = st_tt * tt_t1 + st_dt * dt_t1 + st_ut * ut_t1;
  st_t2 = st_tt * tt_t2 + st_dt * dt_t2 + st_ut * ut_t2;
  st_d1 = st_tt * tt_d1 + st_dt * dt_d1 + st_ut * ut_d1;
  st_d2 = st_tt * tt_d2 + st_dt * dt_d2 + st_ut * ut_d2;
  st_u1 = st_tt * tt_u1 + st_dt * dt_u1 + st_ut * ut_u1;
  st_u2 = st_tt * tt_u2 + st_dt * dt_u2 + st_ut * ut_u2;
  st_ms = st_tt * tt_ms + st_dt * dt_ms + st_ut * ut_ms + st_ms;
  st_re = st_tt * tt_re + st_dt * dt_re + st_ut * ut_re + st_re;
  st_xf = st_tt * tt_xf + st_dt * dt_xf + st_ut * ut_xf;

  state.station2.param.amplz = 0.0;
  state.station2.param.sz = st;

  //---- recalculate turbulent secondary "t" variables using proper cti
  state.station2 = blvar(state.station2, FlowRegimeEnum::Turbulent);

  state.stepbl();
  xfoil.restoreblData(2);

  //---- calculate xt-x2 midpoint cfm value
  SkinFrictionCoefficients turbulentSkinFriction =
      blmid(xfoil, FlowRegimeEnum::Turbulent);

  //---- set up newton system for dct, dth, dds, due, dxi  at  xt and x2
  blc = xfoil.blDiffSolver.solve(FlowRegimeEnum::Turbulent, state, turbulentSkinFriction, xfoil.amcrit);

  //---- convert sensitivities wrt "t" variables into sensitivities
  //-    wrt "1" and "2" variables as done before for the laminar part
  Matrix<double, 5, 5> bt1_right =
      Matrix<double, 5, 5>{{st_a1, st_t1, st_d1, st_u1, st_x1},
                           {tt_a1, tt_t1, tt_d1, tt_u1, tt_x1},
                           {dt_a1, dt_t1, dt_d1, dt_u1, dt_x1},
                           {ut_a1, ut_t1, ut_d1, ut_u1, ut_x1},
                           {xfoil.xt_a1, xfoil.xt_t1, xfoil.xt_d1, xfoil.xt_u1, xfoil.xt_x1}};

  Matrix<double, 5, 5> bt2_right =
      Matrix<double, 5, 5>{{0, st_t2, st_d2, st_u2, st_x2},
                           {0, tt_t2, tt_d2, tt_u2, tt_x2},
                           {0, dt_t2, dt_d2, dt_u2, dt_x2},
                           {0, ut_t2, ut_d2, ut_u2, ut_x2},
                           {0, xfoil.xt_t2, xfoil.xt_d2, xfoil.xt_u2, xfoil.xt_x2}};
  bt1.block(0, 0, 3, 5) = blc.a1.block(0, 0, 3, 5) * bt1_right;
  bt2.block(0, 0, 3, 5) = blc.a1.block(0, 0, 3, 5) * bt2_right;
  bt2 += blc.a2;
  for (int k = 0; k < 3; k++) {
    btrez[k] = blc.rhs[k];
    btm[k] = blc.d_msq[k] + blc.a1(k, 0) * st_ms + blc.a1(k, 1) * tt_ms +
             blc.a1(k, 2) * dt_ms + blc.a1(k, 3) * ut_ms + blc.a1(k, 4) * xfoil.xt_ms;
    btr[k] = blc.d_re[k] + blc.a1(k, 0) * st_re + blc.a1(k, 1) * tt_re +
             blc.a1(k, 2) * dt_re + blc.a1(k, 3) * ut_re + blc.a1(k, 4) * xfoil.xt_re;
    btx[k] = blc.d_xi[k] + blc.a1(k, 0) * st_xf + blc.a1(k, 1) * tt_xf +
             blc.a1(k, 2) * dt_xf + blc.a1(k, 3) * ut_xf + blc.a1(k, 4) * xfoil.xt_xf;
  }

  //---- add up laminar and turbulent parts to get final system
  //-    in terms of honest-to-god "1" and "2" variables.
  blc.rhs[0] = btrez[0];
  blc.rhs[1] = blrez[1] + btrez[1];
  blc.rhs[2] = blrez[2] + btrez[2];
  blc.d_msq[0] = btm[0];
  blc.d_msq[1] = blm[1] + btm[1];
  blc.d_msq[2] = blm[2] + btm[2];
  blc.d_re[0] = btr[0];
  blc.d_re[1] = blr[1] + btr[1];
  blc.d_re[2] = blr[2] + btr[2];
  blc.d_xi[0] = btx[0];
  blc.d_xi[1] = blx[1] + btx[1];
  blc.d_xi[2] = blx[2] + btx[2];
  blc.a1.row(0) = bt1.row(0);
  blc.a2.row(0) = bt2.row(0);
  blc.a1.middleRows(1, 2) = bl1.middleRows(1, 2) + bt1.middleRows(1, 2);
  blc.a2.middleRows(1, 2) = bl2.middleRows(1, 2) + bt2.middleRows(1, 2);

  //---- to be sanitary, restore "1" quantities which got clobbered
  //-    in all of the numerical gymnastics above.  the "2" variables
  //-    were already restored for the xt-x2 differencing part.
  //	for (icom=1; icom<=ncom;icom++){
  //		com1[icom] = c1sav[icom];
  //	}
  xfoil.restoreblData(1);

  return true;
}

bool BoundaryLayerWorkflow::iblpan(XFoil& xfoil) {
  std::stringstream ss;
  const int point_count = xfoil.foil.foil_shape.n;

  for (int i = 0; i <= stagnationIndex; i++) {
    lattice.top.stationToPanel[i] = stagnationIndex - i;
    lattice.top.panelInfluenceFactor[i] = 1.0;
  }

  lattice.top.trailingEdgeIndex = stagnationIndex;
  lattice.top.stationCount = lattice.top.trailingEdgeIndex + 2;

  for (int index = 0; index <= point_count - stagnationIndex; ++index) {
    lattice.bottom.stationToPanel[index] = stagnationIndex + 1 + index;
    lattice.bottom.panelInfluenceFactor[index] = -1.0;
  }

  lattice.bottom.trailingEdgeIndex = point_count - stagnationIndex - 2;

  for (int iw = 0; iw < xfoil.foil.wake_shape.n; iw++) {
    const int panel = point_count + iw;
    const int index = lattice.bottom.trailingEdgeIndex + iw + 2;
    lattice.bottom.stationToPanel[index - 1] = panel;
    lattice.bottom.panelInfluenceFactor[index - 1] = -1.0;
  }

  lattice.bottom.stationCount = lattice.bottom.trailingEdgeIndex + xfoil.foil.wake_shape.n + 2;

  for (int iw = 0; iw < xfoil.foil.wake_shape.n; iw++) {
    lattice.top.stationToPanel[lattice.top.trailingEdgeIndex + iw + 1] =
        lattice.bottom.stationToPanel[lattice.bottom.trailingEdgeIndex + iw + 1];
    lattice.top.panelInfluenceFactor[lattice.top.trailingEdgeIndex + iw + 1] = 1.0;
  }

  const int iblmax = std::max(lattice.top.trailingEdgeIndex, lattice.bottom.trailingEdgeIndex) +
                     xfoil.foil.wake_shape.n + 2;
  if (iblmax > IVX) {
    ss << "iblpan :  ***  bl array overflow\n";
    ss << "Increase IVX to at least " << iblmax << "\n";
    xfoil.writeString(ss.str());
    return false;
  }

  xfoil.lipan = true;
  return true;
}

bool BoundaryLayerWorkflow::iblsys(XFoil& xfoil) {
  int iv = 0;
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 0; ibl < lattice.get(is).stationCount - 1; ++ibl) {
      ++iv;
      lattice.get(is).stationToSystem[ibl] = iv;
    }
  }

  xfoil.nsys = iv;
  if (xfoil.nsys > 2 * IVX) {
    xfoil.writeString("*** iblsys: bl system array overflow. ***");
    return false;
  }

  return true;
}

bool BoundaryLayerWorkflow::stfind(XFoil& xfoil) {
  int stagnation_index = 0;
  bool found = false;
  const int point_count = xfoil.foil.foil_shape.n;
  for (int i = 0; i < point_count - 1; ++i) {
    if (xfoil.surface_vortex(0, i) >= 0.0 &&
        xfoil.surface_vortex(0, i + 1) < 0.0) {
      stagnation_index = i;
      found = true;
      break;
    }
  }

  if (!found) {
    xfoil.writeString("stfind: Stagnation point not found. Continuing ...\n");
    stagnation_index = point_count / 2;
  }

  stagnationIndex = stagnation_index;
  const double dgam = xfoil.surface_vortex(0, stagnation_index + 1) -
                      xfoil.surface_vortex(0, stagnation_index);
  const double ds = xfoil.foil.foil_shape.spline_length[stagnation_index + 1] -
                    xfoil.foil.foil_shape.spline_length[stagnation_index];

  if (xfoil.surface_vortex(0, stagnation_index) <
      -xfoil.surface_vortex(0, stagnation_index + 1)) {
    xfoil.sst = xfoil.foil.foil_shape.spline_length[stagnation_index] -
                ds * (xfoil.surface_vortex(0, stagnation_index) / dgam);
  } else {
    xfoil.sst =
        xfoil.foil.foil_shape.spline_length[stagnation_index + 1] -
        ds * (xfoil.surface_vortex(0, stagnation_index + 1) / dgam);
  }

  if (xfoil.sst <= xfoil.foil.foil_shape.spline_length[stagnation_index])
    xfoil.sst =
        xfoil.foil.foil_shape.spline_length[stagnation_index] + 0.0000001;
  if (xfoil.sst >= xfoil.foil.foil_shape.spline_length[stagnation_index + 1])
    xfoil.sst = xfoil.foil.foil_shape.spline_length[stagnation_index + 1] -
                0.0000001;

  xfoil.sst_go =
      (xfoil.sst - xfoil.foil.foil_shape.spline_length[stagnation_index + 1]) /
      dgam;
  xfoil.sst_gp =
      (xfoil.foil.foil_shape.spline_length[stagnation_index] - xfoil.sst) /
      dgam;

  return true;
}

bool BoundaryLayerWorkflow::stmove(XFoil& xfoil) {
  const int previous = stagnationIndex;
  stfind(xfoil);

  if (previous == stagnationIndex) {
    xfoil.xicalc();
  } else {
    iblpan(xfoil);
    uicalc(xfoil);
    xfoil.xicalc();
    iblsys(xfoil);

    if (stagnationIndex > previous) {
      const int delta = stagnationIndex - previous;

      lattice.top.transitionIndex += delta;
      lattice.bottom.transitionIndex -= delta;

      for (int ibl = lattice.top.stationCount - 2; ibl >= delta; --ibl) {
        copyStationState(1, ibl, ibl - delta);
      }

      const double dudx =
          lattice.top.profiles.edgeVelocity[delta] / lattice.top.arcLengthCoordinates[delta];
      for (int ibl = delta; ibl >= 1; --ibl) {
        copyStationState(1, ibl - 1, delta);
        lattice.top.profiles.edgeVelocity[ibl - 1] = dudx * lattice.top.arcLengthCoordinates[ibl - 1];
      }

      for (int ibl = 0; ibl < lattice.bottom.stationCount - 1; ++ibl) {
        copyStationState(2, ibl, ibl + delta);
      }
    } else {
      const int delta = previous - stagnationIndex;

      lattice.top.transitionIndex -= delta;
      lattice.bottom.transitionIndex += delta;

      for (int ibl = lattice.bottom.stationCount - 1; ibl >= delta + 1; --ibl) {
        copyStationState(2, ibl - 1, (ibl - delta) - 1);
      }

      const double dudx =
          lattice.bottom.profiles.edgeVelocity[delta] / lattice.bottom.arcLengthCoordinates[delta];
      for (int ibl = delta; ibl >= 1; --ibl) {
        copyStationState(2, ibl - 1, delta);
        lattice.bottom.profiles.edgeVelocity[ibl - 1] = dudx * lattice.bottom.arcLengthCoordinates[ibl - 1];
      }

      for (int ibl = 0; ibl < lattice.top.stationCount - 1; ++ibl) {
        copyStationState(1, ibl, ibl + delta);
      }
    }
  }

  for (int is = 1; is <= 2; ++is) {
    for (int ibl = 0; ibl < lattice.get(is).stationCount - 1; ++ibl) {
      lattice.get(is).profiles.massFlux[ibl] =
          lattice.get(is).profiles.displacementThickness[ibl] * lattice.get(is).profiles.edgeVelocity[ibl];
    }
  }

  return true;
}

bool BoundaryLayerWorkflow::uicalc(XFoil& xfoil) {
  //--------------------------------------------------------------
  //     sets inviscid ue from panel inviscid tangential velocity
  //--------------------------------------------------------------
  for (int side = 1; side <= 2; ++side) {
    lattice.get(side).inviscidEdgeVelocity[0] = 0.0;
    lattice.get(side).inviscidEdgeVelocityDerivative[0] = 0.0;
    for (int stationIndex = 0; stationIndex < lattice.get(side).stationCount - 1; ++stationIndex) {
      const int panelIndex = lattice.get(side).stationToPanel[stationIndex];
      lattice.get(side).inviscidEdgeVelocity[stationIndex] = lattice.get(side).panelInfluenceFactor[stationIndex] * xfoil.qinv[panelIndex];
      lattice.get(side).inviscidEdgeVelocityDerivative[stationIndex] = lattice.get(side).panelInfluenceFactor[stationIndex] * xfoil.qinv_a[panelIndex];
    }
  }

  return true;
}

bool BoundaryLayerWorkflow::tesys(XFoil& xfoil, double cte, double tte, double dte) {
  blc.clear();

  state.station2 = this->blvar(state.station2, FlowRegimeEnum::Wake);

  blc.a1(0, 0) = -1.0;
  blc.a2(0, 0) = 1.0;
  blc.rhs[0] = cte - state.station2.param.sz;

  blc.a1(1, 1) = -1.0;
  blc.a2(1, 1) = 1.0;
  blc.rhs[1] = tte - state.station2.param.tz;

  blc.a1(2, 2) = -1.0;
  blc.a2(2, 2) = 1.0;
  blc.rhs[2] = dte - state.station2.param.dz - state.station2.param.dwz;

  return true;
}

void BoundaryLayerWorkflow::copyStationState(int side, int destination, int source) {
  lattice.get(side).profiles.skinFrictionCoeff[destination] =
      lattice.get(side).profiles.skinFrictionCoeff[source];
  lattice.get(side).profiles.momentumThickness[destination] =
      lattice.get(side).profiles.momentumThickness[source];
  lattice.get(side).profiles.displacementThickness[destination] =
      lattice.get(side).profiles.displacementThickness[source];
  lattice.get(side).profiles.edgeVelocity[destination] =
      lattice.get(side).profiles.edgeVelocity[source];
}
