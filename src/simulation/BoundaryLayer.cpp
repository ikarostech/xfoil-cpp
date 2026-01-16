#include "BoundaryLayer.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

#include "XFoil.h"
#include "domain/boundary_layer/boundary_layer_builder.hpp"
#include "domain/coefficient/skin_friction.hpp"
#include "core/boundary_layer_util.hpp"

using BoundaryContext = BoundaryLayerWorkflow::MixedModeStationContext;
using Eigen::Matrix;
using Eigen::Vector;
using Eigen::Vector2d;
using Eigen::VectorXd;

double BoundaryLayerWorkflow::adjustDisplacementForHkLimit(
    double displacementThickness, double momentumThickness, double msq,
    double hklim) {
  const double h = displacementThickness / momentumThickness;

  boundary_layer::KineticShapeParameterResult hkin_result =
      boundary_layer::hkin(h, msq);

  const double dh = std::max(0.0, hklim - hkin_result.hk) / hkin_result.hk_h;
  return displacementThickness + dh * momentumThickness;
}

bool BoundaryLayerWorkflow::blkin(XFoil& xfoil, BoundaryLayerState& state) {
  //----------------------------------------------------------
  //     calculates turbulence-independent secondary "2"
  //     variables from the primary "2" variables.
  //----------------------------------------------------------
  blData& current = state.current();
  //---- set edge mach number ** 2
  current.param.mz =
      current.param.uz * current.param.uz * xfoil.blCompressibility.hstinv /
      (xfoil.blCompressibility.gm1bl *
       (1.0 - 0.5 * current.param.uz * current.param.uz * xfoil.blCompressibility.hstinv));
  double tr2 = 1.0 + 0.5 * xfoil.blCompressibility.gm1bl * current.param.mz;
  current.param.mz_uz = 2.0 * current.param.mz * tr2 / current.param.uz;
  current.param.mz_ms =
      current.param.uz * current.param.uz * tr2 /
      (xfoil.blCompressibility.gm1bl *
       (1.0 - 0.5 * current.param.uz * current.param.uz * xfoil.blCompressibility.hstinv)) *
      xfoil.blCompressibility.hstinv_ms;

  //---- set edge density (isentropic relation)
  current.param.rz =
      xfoil.blCompressibility.rstbl *
      pow(tr2, (-1.0 / xfoil.blCompressibility.gm1bl));
  current.param.rz_uz = -current.param.rz / tr2 * 0.5 * current.param.mz_uz;
  current.param.rz_ms = -current.param.rz / tr2 * 0.5 * current.param.mz_ms +
                        xfoil.blCompressibility.rstbl_ms *
                            pow(tr2, (-1.0 / xfoil.blCompressibility.gm1bl));

  //---- set shape parameter
  current.param.hz = current.param.dz / current.param.tz;
  current.param.hz_dz = 1.0 / current.param.tz;
  current.param.hz_tz = -current.param.hz / current.param.tz;

  //---- set edge static/stagnation enthalpy
  double herat =
      1.0 - 0.5 * current.param.uz * current.param.uz * xfoil.blCompressibility.hstinv;
  double he_u2 = -current.param.uz * xfoil.blCompressibility.hstinv;
  double he_ms =
      -0.5 * current.param.uz * current.param.uz * xfoil.blCompressibility.hstinv_ms;
  //---- set molecular viscosity
  double v2_he = (1.5 / herat - 1.0 / (herat + kHvrat));

  //---- set kinematic shape parameter
  boundary_layer::KineticShapeParameterResult hkin_result =
      boundary_layer::hkin(current.param.hz, current.param.mz);
  current.hkz.scalar = hkin_result.hk;

  current.hkz.u() = hkin_result.hk_msq * current.param.mz_uz;
  current.hkz.t() = hkin_result.hk_h * current.param.hz_tz;
  current.hkz.d() = hkin_result.hk_h * current.param.hz_dz;
  current.hkz.ms() = hkin_result.hk_msq * current.param.mz_ms;

  //---- set momentum thickness reynolds number
  current.rtz.scalar =
      current.param.rz * current.param.uz * current.param.tz /
      (sqrt(herat * herat * herat) * (1.0 + kHvrat) / (herat + kHvrat) /
       xfoil.blReynolds.reybl);
  current.rtz.u() = current.rtz.scalar *
                    (1.0 / current.param.uz +
                     current.param.rz_uz / current.param.rz - v2_he * he_u2);
  current.rtz.t() = current.rtz.scalar / current.param.tz;
  current.rtz.ms() =
      current.rtz.scalar * (current.param.rz_ms / current.param.rz +
                            (1 / xfoil.blReynolds.reybl * xfoil.blReynolds.reybl_ms -
                             v2_he * he_ms));
  current.rtz.re() =
      current.rtz.scalar * (xfoil.blReynolds.reybl_re / xfoil.blReynolds.reybl);

  return true;
}

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
    tesys(lattice.top.profiles, lattice.bottom.profiles, xfoil.foil.edge);
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
        uem_sq * xfoil.blCompressibility.hstinv /
        (xfoil.blCompressibility.gm1bl *
         (1.0 - 0.5 * uem_sq * xfoil.blCompressibility.hstinv));
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
  const double msq =
      uei_sq * xfoil.blCompressibility.hstinv /
      (xfoil.blCompressibility.gm1bl *
       (1.0 - 0.5 * uei_sq * xfoil.blCompressibility.hstinv));
  double dsw = ctx.dsi - ctx.dswaki;
  dsw = adjustDisplacementForHkLimit(dsw, ctx.thi, msq, hklim);
  ctx.dsi = dsw + ctx.dswaki;

  return ctx.dmax <= deps;
}

blData BoundaryLayerWorkflow::blvar(blData data, FlowRegimeEnum flowRegimeType) {
  return this->boundaryLayerVariablesSolver.solve(data, flowRegimeType);
}

SkinFrictionCoefficients BoundaryLayerWorkflow::blmid(
    FlowRegimeEnum flowRegimeType) {
  BoundaryLayerState& state = this->state;
  blData& previous = state.previous();
  blData& current = state.current();

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
      uei * (1.0 - xfoil.blCompressibility.tkbl) /
      (1.0 - xfoil.blCompressibility.tkbl * (uei / xfoil.blCompressibility.qinfbl) *
                 (uei / xfoil.blCompressibility.qinfbl));
  data.param.uz_uei =
      (1.0 + xfoil.blCompressibility.tkbl *
                (2.0 * data.param.uz * uei / xfoil.blCompressibility.qinfbl /
                     xfoil.blCompressibility.qinfbl -
                 1.0)) /
      (1.0 - xfoil.blCompressibility.tkbl *
                 (uei / xfoil.blCompressibility.qinfbl) *
                 (uei / xfoil.blCompressibility.qinfbl));
  data.param.uz_ms =
      (data.param.uz * (uei / xfoil.blCompressibility.qinfbl) *
           (uei / xfoil.blCompressibility.qinfbl) -
       uei) *
      xfoil.blCompressibility.tkbl_ms /
      (1.0 - xfoil.blCompressibility.tkbl *
                 (uei / xfoil.blCompressibility.qinfbl) *
                 (uei / xfoil.blCompressibility.qinfbl));
  return data;
}

bool BoundaryLayerWorkflow::blsys(XFoil& xfoil) {
  blData& previous = state.previous();
  blData& current = state.current();

  SkinFrictionCoefficients skinFriction = blmid(xfoil.flowRegime);
  current = blvar(current, xfoil.flowRegime);

  if (xfoil.flowRegime == FlowRegimeEnum::Similarity) {
    state.stepbl();
  }

  if (xfoil.flowRegime == FlowRegimeEnum::Transition) {
    trdif(xfoil);
  } else {
    blc = blDiffSolver.solve(xfoil.flowRegime, state, skinFriction,
                             xfoil.blTransition.amcrit);
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

  double ctr, ctr_hk2;

  boundaryLayerStore.saveblData(state.station1, 1);
  boundaryLayerStore.saveblData(state.station2, 2);

  //---- weighting factors for linear interpolation to transition point
  blDiff wf;
  wf.scalar = (xt.scalar - state.station1.param.xz) /
                  (state.station2.param.xz - state.station1.param.xz);
  
  double wf2 = (xt.scalar - state.station1.param.xz) / (state.station2.param.xz - state.station1.param.xz);
  double wf_xt = 1.0 / (state.station2.param.xz - state.station1.param.xz);
  wf.vector = xt.vector * wf_xt;
  wf.x1() += (wf2 - 1.0) / (state.station2.param.xz - state.station1.param.xz);
  wf.x2() -= wf2 / (state.station2.param.xz - state.station1.param.xz);

  double wf1 = 1.0 - wf2;

  //-----interpolate primary variables to transition point
  blDiff tt;
  tt.scalar = state.station1.param.tz * wf1 + state.station2.param.tz * wf2;
  tt.vector = (state.station2.param.tz - state.station1.param.tz) * wf.vector;
  tt.t1() += wf1;
  tt.t2() += wf2;

  blDiff dt;
  dt.scalar = state.station1.param.dz * wf1 + state.station2.param.dz * wf2;
  dt.vector = (state.station2.param.dz - state.station1.param.dz) * wf.vector;
  dt.d1() += wf1;
  dt.d2() += wf2;

  blDiff ut;
  ut.scalar = state.station1.param.uz * wf1 + state.station2.param.uz * wf2;
  ut.vector = (state.station2.param.uz - state.station1.param.uz) * wf.vector;
  ut.u1() += wf1;
  ut.u2() += wf2;

  //---- set primary "t" variables at xt  (really placed into "2" variables)
  state.station2.param.xz = xt.scalar;
  state.station2.param.tz = tt.scalar;
  state.station2.param.dz = dt.scalar;
  state.station2.param.uz = ut.scalar;

  state.station2.param.amplz = xfoil.blTransition.amcrit;
  state.station2.param.sz = 0.0;

  //---- calculate laminar secondary "t" variables
  blkin(xfoil, state);
  state.station2 = blvar(state.station2, FlowRegimeEnum::Laminar);

  //---- calculate x1-xt midpoint cfm value
  SkinFrictionCoefficients laminarSkinFriction =
      blmid(FlowRegimeEnum::Laminar);

  //=    at this point, all "2" variables are really "t" variables at xt

  //---- set up newton system for dam, dth, dds, due, dxi  at  x1 and xt
  blc = blDiffSolver.solve(FlowRegimeEnum::Laminar, state, laminarSkinFriction,
                           xfoil.blTransition.amcrit);

  //---- the current newton system is in terms of "1" and "t" variables,
  //-    so calculate its equivalent in terms of "1" and "2" variables.
  //-    in other words, convert residual sensitivities wrt "t" variables
  //-    into sensitivities wrt "1" and "2" variables.  the amplification
  //-    equation is unnecessary here, so the k=1 row is left empty.
  blrez = blc.rhs;
  for (int k = 1; k < 3; k++) {
    blm[k] = blc.d_msq[k] + blc.a2(k, 1) * tt.ms() + blc.a2(k, 2) * dt.ms() +
             blc.a2(k, 3) * ut.ms() + blc.a2(k, 4) * xt.ms();
    blr[k] = blc.d_re[k] + blc.a2(k, 1) * tt.re() + blc.a2(k, 2) * dt.re() +
             blc.a2(k, 3) * ut.re() + blc.a2(k, 4) * xt.re();
    blx[k] = blc.d_xi[k] + blc.a2(k, 1) * tt.xf() + blc.a2(k, 2) * dt.xf() +
             blc.a2(k, 3) * ut.xf() + blc.a2(k, 4) * xt.xf();
  }
  const Eigen::Matrix<double, 4, 5> bl1_transform{
    {tt.a(), tt.t1(), tt.d1(), tt.u1(), tt.x1()},
    {dt.a(), dt.t1(), dt.d1(), dt.u1(), dt.x1()},
    {ut.a(), ut.t1(), ut.d1(), ut.u1(), ut.x1()},
    {xt.a(), xt.t1(), xt.d1(), xt.u1(), xt.x1()}
  };
  bl1.block<2, 5>(1, 0) =
      blc.a1.middleRows<2>(1) + blc.a2.block<2, 4>(1, 1) * bl1_transform;

  const Eigen::Matrix<double, 4, 4> bl2_transform{
  {tt.t2(), tt.d2(), tt.u2(), tt.x2()},
  {dt.t2(), dt.d2(), dt.u2(), dt.x2()},
  {ut.t2(), ut.d2(), ut.u2(), ut.x2()},
  {xt.t2(), xt.d2(), xt.u2(), xt.x2()}
};

  bl2.block<2, 1>(1, 0).setZero();
  bl2.block<2, 4>(1, 1) = blc.a2.block<2, 4>(1, 1) * bl2_transform;

  //**** second, set up turbulent part between xt and x2  ****

  //---- calculate equilibrium shear coefficient cqt at transition point
  state.station2 = blvar(state.station2, FlowRegimeEnum::Turbulent);

  //---- set initial shear coefficient value st at transition point
  //-    ( note that cq2, cq2_t2, etc. are really "cqt", "cqt_tt", etc.)

  ctr = 1.8 * exp(-3.3 / (state.station2.hkz.scalar - 1.0));
  ctr_hk2 = ctr * 3.3 / (state.station2.hkz.scalar - 1.0) / (state.station2.hkz.scalar - 1.0);

  double st = ctr * state.station2.cqz.scalar;
  double st_tt =
      ctr * state.station2.cqz.t() + state.station2.cqz.scalar * ctr_hk2 * state.station2.hkz.t();
  double st_dt =
      ctr * state.station2.cqz.d() + state.station2.cqz.scalar * ctr_hk2 * state.station2.hkz.d();
  double st_ut =
      ctr * state.station2.cqz.u() + state.station2.cqz.scalar * ctr_hk2 * state.station2.hkz.u();
  double st_ms =
      ctr * state.station2.cqz.ms() + state.station2.cqz.scalar * ctr_hk2 * state.station2.hkz.ms();
  double st_re = ctr * state.station2.cqz.re();

  state.station2.param.amplz = 0.0;
  state.station2.param.sz = st;

  //---- recalculate turbulent secondary "t" variables using proper cti
  state.station2 = blvar(state.station2, FlowRegimeEnum::Turbulent);

  state.stepbl();
  state.station2 = boundaryLayerStore.restoreblData(2);

  //---- calculate xt-x2 midpoint cfm value
  SkinFrictionCoefficients turbulentSkinFriction =
      blmid(FlowRegimeEnum::Turbulent);

  //---- set up newton system for dct, dth, dds, due, dxi  at  xt and x2
  blc = blDiffSolver.solve(FlowRegimeEnum::Turbulent, state, turbulentSkinFriction,
                           xfoil.blTransition.amcrit);

  //---- convert sensitivities wrt "t" variables into sensitivities
  //-    wrt "1" and "2" variables as done before for the laminar part
  Eigen::VectorXd common_st = Eigen::Vector3d{st_tt, st_dt, st_ut};
  Eigen::VectorXd st1 = 
    Eigen::Matrix<double, 5, 3>{
      {tt.a(), dt.a(), ut.a()},
      {tt.t1(), dt.t1(), ut.t1()},
      {tt.d1(), dt.d1(), ut.d1()},
      {tt.u1(), dt.u1(), ut.u1()},
      {tt.x1(), dt.x1(), ut.x1()}
    } * common_st;
  Eigen::VectorXd st2 = 
    Eigen::Matrix<double, 5, 3>{
      {0, 0, 0},
      {tt.t2(), dt.t2(), ut.t2()},
      {tt.d2(), dt.d2(), ut.d2()},
      {tt.u2(), dt.u2(), ut.u2()},
      {tt.x2(), dt.x2(), ut.x2()}
    } * common_st;

  st_ms = st_tt * tt.ms() + st_dt * dt.ms() + st_ut * ut.ms() + st_ms;
  st_re = st_tt * tt.re() + st_dt * dt.re() + st_ut * ut.re() + st_re;
  double st_xf = st_tt * tt.xf() + st_dt * dt.xf() + st_ut * ut.xf();

  Matrix<double, 5, 5> bt1_right = Matrix<double, 5, 5>::Zero();
  bt1_right.block<4, 5>(1, 0) = bl1_transform;
  bt1_right.row(0) = st1.transpose();

  Matrix<double, 5, 5> bt2_right = Matrix<double, 5, 5>::Zero();
  bt2_right.block<4, 4>(1, 1) = bl2_transform;
  bt2_right.row(0) = st2.transpose();
  bt1.block(0, 0, 3, 5) = blc.a1.block(0, 0, 3, 5) * bt1_right;
  bt2.block(0, 0, 3, 5) = blc.a1.block(0, 0, 3, 5) * bt2_right;
  bt2 += blc.a2;
  for (int k = 0; k < 3; k++) {
    btrez[k] = blc.rhs[k];
    btm[k] = blc.d_msq[k] + blc.a1(k, 0) * st_ms + blc.a1(k, 1) * tt.ms() +
             blc.a1(k, 2) * dt.ms() + blc.a1(k, 3) * ut.ms() + blc.a1(k, 4) * xt.ms();
    btr[k] = blc.d_re[k] + blc.a1(k, 0) * st_re + blc.a1(k, 1) * tt.re() +
             blc.a1(k, 2) * dt.re() + blc.a1(k, 3) * ut.re() + blc.a1(k, 4) * xt.re();
    btx[k] = blc.d_xi[k] + blc.a1(k, 0) * st_xf + blc.a1(k, 1) * tt.xf() +
             blc.a1(k, 2) * dt.xf() + blc.a1(k, 3) * ut.xf() + blc.a1(k, 4) * xt.xf();
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

  state.station1 = boundaryLayerStore.restoreblData(1);

  return true;
}

double BoundaryLayerWorkflow::computeTransitionLocation(double weightingFactor) const {
  const double upstreamLocation = state.station1.param.xz;
  const double downstreamLocation = state.station2.param.xz;
  return upstreamLocation * (1.0 - weightingFactor) + downstreamLocation * weightingFactor;
}

double BoundaryLayerWorkflow::calcHtarg(int ibl, int is, bool wake) {
  if (ibl < lattice.get(is).transitionIndex) {
    return state.station1.hkz.scalar +
           0.03 * (state.station2.param.xz - state.station1.param.xz) / state.station1.param.tz;
  }

  if (ibl == lattice.get(is).transitionIndex) {
    return state.station1.hkz.scalar +
           (0.03 * (xt.scalar - state.station1.param.xz) -
            0.15 * (state.station2.param.xz - xt.scalar)) /
               state.station1.param.tz;
  }

  if (wake) {
    const double cst =
        0.03 * (state.station2.param.xz - state.station1.param.xz) / state.station1.param.tz;
    auto euler = [](double hk2, double hk1, double cst_local) {
      return hk2 - (hk2 + cst_local * pow(hk2 - 1, 3) - hk1) /
                       (1 + 3 * cst_local * pow(hk2 - 1, 2));
    };
    state.station2.hkz.scalar = state.station1.hkz.scalar;
    for (int i = 0; i < 3; i++) {
      state.station2.hkz.scalar =
          euler(state.station2.hkz.scalar, state.station1.hkz.scalar, cst);
    }
    return state.station2.hkz.scalar;
  }

  return state.station1.hkz.scalar -
         0.15 * (state.station2.param.xz - state.station1.param.xz) / state.station1.param.tz;
}

bool BoundaryLayerWorkflow::trchek(XFoil& xfoil) {
  //----------------------------------------------------------------
  //     new second-order version:  december 1994.
  //
  //     checks if transition occurs in the current interval x1..x2.
  //     if transition occurs, then set transition location xt, and
  //     its sensitivities to "1" and "2" variables.  if no transition,
  //     set amplification ampl2.
  //
  //     solves the implicit amplification equation for n2:
  //
  //       n2 - n1     n'(xt,nt) + n'(x1,n1)
  //       -------  =  ---------------------
  //       x2 - x1               2
  //
  //     in effect, a 2-point central difference is used between
  //     x1..x2 (no transition), or x1..xt (transition).  the switch
  //     is done by defining xt,nt in the equation above depending
  //     on whether n2 exceeds ncrit.
  //
  //  if n2<ncrit:  nt=n2    , xt=x2                  (no transition)
  //
  //  if n2>ncrit:  nt=ncrit , xt=(ncrit-n1)/(n2-n1)  (transition)
  //
  //----------------------------------------------------------------
  double amplt, sfa, sfa_a1, sfa_a2, sfx;
  double sfx_x1, sfx_x2, sfx_xf;
  double tt, dt, ut, amsave;
  double res = 0.0, res_a2 = 0.0;
  double da2 = 0.0, dxt = 0.0, tt_t1 = 0.0, dt_d1 = 0.0, ut_u1 = 0.0;
  double tt_t2 = 0.0, dt_d2 = 0.0, ut_u2 = 0.0, tt_a1 = 0.0, dt_a1 = 0.0;
  double ut_a1 = 0.0, tt_x1 = 0.0, dt_x1 = 0.0, ut_x1 = 0.0, tt_x2 = 0.0,
         dt_x2 = 0.0, ut_x2 = 0.0;
  double amplt_a2, wf, wf_a1, wf_a2, wf_xf, wf_x1, wf_x2;
  double xt_a2, dt_a2, tt_a2;
  double ut_a2;
  double daeps = 0.00005;

  amplt_a2 = 0.0;
  xt_a2 = dt_a2 = tt_a2 = 0.0;
  ut_a2 = 0.0;

  //---- save variables and sensitivities at ibl ("2") for future restoration
  boundaryLayerStore.saveblData(state.station2, 2);

  //---- calculate average amplification rate ax over x1..x2 interval
  BoundaryLayerUtil::AxResult ax_result =
      BoundaryLayerUtil::axset(state.station1.hkz.scalar, state.station1.param.tz, state.station1.rtz.scalar,
            state.station1.param.amplz, state.station2.hkz.scalar, state.station2.param.tz,
            state.station2.rtz.scalar, state.station2.param.amplz,
            xfoil.blTransition.amcrit);

  //---- set initial guess for iterate n2 (ampl2) at x2
  state.station2.param.amplz = state.station1.param.amplz +
                        ax_result.ax * (state.station2.param.xz - state.station1.param.xz);
  //---- solve implicit system for amplification ampl2
  auto iterateAmplification = [&]() -> bool {
    for (int itam = 0; itam < 30; itam++) {
      //---- define weighting factors wf1,wf2 for defining "t" quantities
      if (state.station2.param.amplz <= xfoil.blTransition.amcrit) {
        //------ there is no transition yet,  "t" is the same as "2"
        amplt = state.station2.param.amplz;
        amplt_a2 = 1.0;
        sfa = 1.0;
        sfa_a1 = 0.0;
        sfa_a2 = 0.0;
      } else {
        //------ there is transition in x1..x2, "t" is set from n1, n2
        amplt = xfoil.blTransition.amcrit;
        amplt_a2 = 0.0;
        sfa = (amplt - state.station1.param.amplz) /
              (state.station2.param.amplz - state.station1.param.amplz);
        sfa_a1 = (sfa - 1.0) / (state.station2.param.amplz - state.station1.param.amplz);
        sfa_a2 = (-sfa) / (state.station2.param.amplz - state.station1.param.amplz);
      }

      if (xfoil.blTransition.xiforc < state.station2.param.xz) {
        sfx = (xfoil.blTransition.xiforc - state.station1.param.xz) /
              (state.station2.param.xz - state.station1.param.xz);
        sfx_x1 = (sfx - 1.0) / (state.station2.param.xz - state.station1.param.xz);
        sfx_x2 = (-sfx) / (state.station2.param.xz - state.station1.param.xz);
        sfx_xf = 1.0 / (state.station2.param.xz - state.station1.param.xz);
      } else {
        sfx = 1.0;
        sfx_x1 = 0.0;
        sfx_x2 = 0.0;
        sfx_xf = 0.0;
      }

      //---- set weighting factor from free or forced transition
      if (sfa < sfx) {
        wf = sfa;
        wf_a1 = sfa_a1;
        wf_a2 = sfa_a2;
        wf_x1 = 0.0;
        wf_x2 = 0.0;
        wf_xf = 0.0;
      } else {
        wf = sfx;
        wf_a1 = 0.0;
        wf_a2 = 0.0;
        wf_x1 = sfx_x1;
        wf_x2 = sfx_x2;
        wf_xf = sfx_xf;
      }

      //---- interpolate bl variables to xt
      xt.scalar = computeTransitionLocation(wf);
      tt = state.station1.param.tz * (1 - wf) + state.station2.param.tz * wf;
      dt = state.station1.param.dz * (1 - wf) + state.station2.param.dz * wf;
      ut = state.station1.param.uz * (1 - wf) + state.station2.param.uz * wf;

      xt_a2 = (state.station2.param.xz - state.station1.param.xz) * wf_a2;
      tt_a2 = (state.station2.param.tz - state.station1.param.tz) * wf_a2;
      dt_a2 = (state.station2.param.dz - state.station1.param.dz) * wf_a2;
      ut_a2 = (state.station2.param.uz - state.station1.param.uz) * wf_a2;

      //---- temporarily set "2" variables from "t" for blkin
      state.station2.param.xz = xt.scalar;
      state.station2.param.tz = tt;
      state.station2.param.dz = dt;
      state.station2.param.uz = ut;

      //---- calculate laminar secondary "t" variables hkt, rtt
      blkin(xfoil, state);

      blData::blVector hkt = state.station2.hkz;
      blData::blVector rtt = state.station2.rtz;

      //---- restore clobbered "2" variables, except for ampl2
      amsave = state.station2.param.amplz;

      state.station2 = boundaryLayerStore.restoreblData(2);

      state.station2.param.amplz = amsave;

      //---- calculate amplification rate ax over current x1-xt interval
      ax_result = BoundaryLayerUtil::axset(state.station1.hkz.scalar, state.station1.param.tz,
                        state.station1.rtz.scalar, state.station1.param.amplz, hkt.scalar, tt, rtt.scalar,
                        amplt, xfoil.blTransition.amcrit);

      //---- punch out early if there is no amplification here
      if (ax_result.ax <= 0.0) {
        return true;
      }

      //---- set sensitivity of ax(a2)
      ax_result.ax_a2 =
          (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
           ax_result.ax_rt2 * rtt.t()) *
              tt_a2 +
          (ax_result.ax_hk2 * hkt.d()) * dt_a2 +
          (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_a2 +
          ax_result.ax_a2 * amplt_a2;

      //---- residual for implicit ampl2 definition (amplification equation)
      res = state.station2.param.amplz - state.station1.param.amplz -
            ax_result.ax * (state.station2.param.xz - state.station1.param.xz);
      res_a2 = 1.0 - ax_result.ax_a2 * (state.station2.param.xz - state.station1.param.xz);

      da2 = -res / res_a2;

      xfoil.rlx = 1.0;
      dxt = xt_a2 * da2;

      if (xfoil.rlx * fabs(dxt / (state.station2.param.xz - state.station1.param.xz)) > 0.05) {
        xfoil.rlx = 0.05 * fabs((state.station2.param.xz - state.station1.param.xz) / dxt);
      }

      if (xfoil.rlx * fabs(da2) > 1.0) {
        xfoil.rlx = 1.0 * fabs(1.0 / da2);
      }

      //---- check if converged
      if (fabs(da2) < daeps) {
        return true;
      }

      if ((state.station2.param.amplz > xfoil.blTransition.amcrit &&
           state.station2.param.amplz + xfoil.rlx * da2 < xfoil.blTransition.amcrit) ||
          (state.station2.param.amplz < xfoil.blTransition.amcrit &&
           state.station2.param.amplz + xfoil.rlx * da2 > xfoil.blTransition.amcrit)) {
        //------ limited newton step so ampl2 doesn't step across amcrit either
        // way
        state.station2.param.amplz = xfoil.blTransition.amcrit;
      } else {
        //------ regular newton step
        state.station2.param.amplz = state.station2.param.amplz + xfoil.rlx * da2;
      }
    }
    return false;
  };

  if (!iterateAmplification()) {
    // TRACE("trchek2 - n2 convergence failed\n");
    xfoil.writeString("trchek2 - n2 convergence failed\n");
    if (XFoil::isCancelled())
      return false;
  }

  //---- test for free or forced transition
  xfoil.trfree = (state.station2.param.amplz >= xfoil.blTransition.amcrit);
  xfoil.trforc =
      (xfoil.blTransition.xiforc > state.station1.param.xz) &&
      (xfoil.blTransition.xiforc <= state.station2.param.xz);

  //---- set transition interval flag
  const bool transitionDetected = (xfoil.trforc || xfoil.trfree);
  xfoil.flowRegime = transitionDetected ? FlowRegimeEnum::Transition
                                   : FlowRegimeEnum::Laminar;

  if (!transitionDetected)
    return false;

  //---- resolve if both forced and free transition
  if (xfoil.trfree && xfoil.trforc) {
    xfoil.trforc = xfoil.blTransition.xiforc < xt.scalar;
    xfoil.trfree = xfoil.blTransition.xiforc >= xt.scalar;
  }

  if (xfoil.trforc) {
    //----- if forced transition, then xt is prescribed,
    //-     no sense calculating the sensitivities, since we know them...
    xt.scalar = xfoil.blTransition.xiforc;
    xt.a() = 0.0;
    xt.x1() = 0.0;
    xt.t1() = 0.0;
    xt.d1() = 0.0;
    xt.u1() = 0.0;
    xt.x2() = 0.0;
    xt.t2() = 0.0;
    xt.d2() = 0.0;
    xt.u2() = 0.0;
    xt.ms() = 0.0;
    xt.re() = 0.0;
    xt.xf() = 1.0;
    return true;
  }

  //---- free transition ... set sensitivities of xt

  xt.x1() = (1 - wf);
  tt_t1 = (1 - wf);
  dt_d1 = (1 - wf);
  ut_u1 = (1 - wf);

  xt.x2() = wf;
  tt_t2 = wf;
  dt_d2 = wf;
  ut_u2 = wf;

  xt.a() = (state.station2.param.xz - state.station1.param.xz) * wf_a1;
  tt_a1 = (state.station2.param.tz - state.station1.param.tz) * wf_a1;
  dt_a1 = (state.station2.param.dz - state.station1.param.dz) * wf_a1;
  ut_a1 = (state.station2.param.uz - state.station1.param.uz) * wf_a1;

  xt.x1() += (state.station2.param.xz - state.station1.param.xz) * wf_x1;
  tt_x1 = (state.station2.param.tz - state.station1.param.tz) * wf_x1;
  dt_x1 = (state.station2.param.dz - state.station1.param.dz) * wf_x1;
  ut_x1 = (state.station2.param.uz - state.station1.param.uz) * wf_x1;

  xt.x2() += (state.station2.param.xz - state.station1.param.xz) * wf_x2;
  tt_x2 = (state.station2.param.tz - state.station1.param.tz) * wf_x2;
  dt_x2 = (state.station2.param.dz - state.station1.param.dz) * wf_x2;
  ut_x2 = (state.station2.param.uz - state.station1.param.uz) * wf_x2;

  xt.xf() = (state.station2.param.xz - state.station1.param.xz) * wf_xf;

  //---- at this point, ax = ax( hk1, t1, rt1, a1, hkt, tt, rtt, at )
  blData::blVector hkt = state.station2.hkz;
  blData::blVector rtt = state.station2.rtz;

  //---- set sensitivities of ax( t1 d1 u1 a1 t2 d2 u2 a2 ms re )
  blDiff ax;
  ax.scalar = 0.0;  // store a2 sensitivity in scalar slot to keep ax terms together

  ax.t1() = ax_result.ax_hk1 * state.station1.hkz.t() + ax_result.ax_t1 +
            ax_result.ax_rt1 * state.station1.rtz.t() +
            (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
             ax_result.ax_rt2 * rtt.t()) *
                tt_t1;
  ax.d1() =
      ax_result.ax_hk1 * state.station1.hkz.d() + (ax_result.ax_hk2 * hkt.d()) * dt_d1;
  ax.u1() = ax_result.ax_hk1 * state.station1.hkz.u() +
            ax_result.ax_rt1 * state.station1.rtz.u() +
            (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_u1;
  ax.a() = ax_result.ax_a1 +
           (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
            ax_result.ax_rt2 * rtt.t()) *
               tt_a1 +
           (ax_result.ax_hk2 * hkt.d()) * dt_a1 +
           (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_a1;
  ax.x1() =
      (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
       ax_result.ax_rt2 * rtt.t()) *
          tt_x1 +
      (ax_result.ax_hk2 * hkt.d()) * dt_x1 +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_x1;

  ax.t2() = (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
             ax_result.ax_rt2 * rtt.t()) *
            tt_t2;
  ax.d2() = (ax_result.ax_hk2 * hkt.d()) * dt_d2;
  ax.u2() = (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_u2;
  ax.scalar =
      ax_result.ax_a2 * amplt_a2 +
      (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
       ax_result.ax_rt2 * rtt.t()) *
          tt_a2 +
      (ax_result.ax_hk2 * hkt.d()) * dt_a2 +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_a2;
  ax.x2() =
      (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
       ax_result.ax_rt2 * rtt.t()) *
          tt_x2 +
      (ax_result.ax_hk2 * hkt.d()) * dt_x2 +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_x2;

  ax.ms() = ax_result.ax_hk2 * hkt.ms() + ax_result.ax_rt2 * rtt.ms() +
            ax_result.ax_hk1 * state.station1.hkz.ms() +
            ax_result.ax_rt1 * state.station1.rtz.ms();
  ax.re() = ax_result.ax_rt2 * rtt.re() + ax_result.ax_rt1 * state.station1.rtz.re();

  //---- set sensitivities of residual res
  blDiff z;
  z.scalar = -(state.station2.param.xz - state.station1.param.xz);
  z.vector = z.scalar * ax.vector;
  z.a() -= 1.0;
  z.x1() += ax_result.ax;
  z.x2() -= ax_result.ax;

  //---- set sensitivities of xt, with res being stationary for a2 constraint
  xt.vector -= (xt_a2 / (z.scalar * ax.scalar + 1.0)) * z.vector;

  return true;
}

bool BoundaryLayerWorkflow::iblpan(int point_count, int wake_point_count,
                                   std::string* error_message) {
  std::stringstream ss;
  if (error_message) {
    error_message->clear();
  }

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

  for (int iw = 0; iw < wake_point_count; iw++) {
    const int panel = point_count + iw;
    const int index = lattice.bottom.trailingEdgeIndex + iw + 2;
    lattice.bottom.stationToPanel[index - 1] = panel;
    lattice.bottom.panelInfluenceFactor[index - 1] = -1.0;
  }

  lattice.bottom.stationCount = lattice.bottom.trailingEdgeIndex + wake_point_count + 2;

  for (int iw = 0; iw < wake_point_count; iw++) {
    lattice.top.stationToPanel[lattice.top.trailingEdgeIndex + iw + 1] =
        lattice.bottom.stationToPanel[lattice.bottom.trailingEdgeIndex + iw + 1];
    lattice.top.panelInfluenceFactor[lattice.top.trailingEdgeIndex + iw + 1] = 1.0;
  }

  const int iblmax = std::max(lattice.top.trailingEdgeIndex, lattice.bottom.trailingEdgeIndex) +
                     wake_point_count + 2;
  if (iblmax > IVX) {
    ss << "iblpan :  ***  bl array overflow\n";
    ss << "Increase IVX to at least " << iblmax << "\n";
    if (error_message) {
      *error_message = ss.str();
    }
    return false;
  }

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

BoundaryLayerWorkflow::StagnationResult BoundaryLayerWorkflow::stfind(
    const Eigen::Matrix2Xd& surface_vortex,
    const Eigen::VectorXd& spline_length) const {
  int stagnation_index = 0;
  bool found = false;
  const int point_count = static_cast<int>(surface_vortex.cols());
  for (int i = 0; i < point_count - 1; ++i) {
    if (surface_vortex(0, i) >= 0.0 &&
        surface_vortex(0, i + 1) < 0.0) {
      stagnation_index = i;
      found = true;
      break;
    }
  }

  if (!found) {
    stagnation_index = point_count / 2;
  }

  StagnationResult result;
  result.stagnationIndex = stagnation_index;
  result.found = found;
  const double dgam = surface_vortex(0, stagnation_index + 1) -
                      surface_vortex(0, stagnation_index);
  const double ds = spline_length[stagnation_index + 1] -
                    spline_length[stagnation_index];

  if (surface_vortex(0, stagnation_index) <
      -surface_vortex(0, stagnation_index + 1)) {
    result.sst = spline_length[stagnation_index] -
                 ds * (surface_vortex(0, stagnation_index) / dgam);
  } else {
    result.sst =
        spline_length[stagnation_index + 1] -
        ds * (surface_vortex(0, stagnation_index + 1) / dgam);
  }

  if (result.sst <= spline_length[stagnation_index])
    result.sst = spline_length[stagnation_index] + 0.0000001;
  if (result.sst >= spline_length[stagnation_index + 1])
    result.sst = spline_length[stagnation_index + 1] - 0.0000001;

  result.sst_go =
      (result.sst - spline_length[stagnation_index + 1]) / dgam;
  result.sst_gp =
      (spline_length[stagnation_index] - result.sst) / dgam;

  return result;
}

bool BoundaryLayerWorkflow::stmove(XFoil& xfoil) {
  const int previous = stagnationIndex;
  const auto stagnation = stfind(xfoil.surface_vortex,
                                 xfoil.foil.foil_shape.spline_length);
  if (!stagnation.found) {
    std::cout << "stfind: Stagnation point not found. Continuing ..." << std::endl;
  }
  stagnationIndex = stagnation.stagnationIndex;
  xfoil.stagnation = stagnation;
  stagnationSst = stagnation.sst;

  if (previous == stagnationIndex) {
    xicalc(xfoil.foil);
  } else {
    std::string iblpan_error;
    if (iblpan(xfoil.foil.foil_shape.n, xfoil.foil.wake_shape.n, &iblpan_error)) {
      xfoil.lipan = true;
    } else if (!iblpan_error.empty()) {
      xfoil.writeString(iblpan_error);
    }
    const auto inviscid_edge_velocity = uicalc(xfoil.qinv_matrix);
    lattice.top.inviscidEdgeVelocityMatrix = inviscid_edge_velocity.top;
    lattice.bottom.inviscidEdgeVelocityMatrix = inviscid_edge_velocity.bottom;
    xicalc(xfoil.foil);
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

SidePair<Eigen::Matrix2Xd> BoundaryLayerWorkflow::uicalc(
    const Eigen::Matrix2Xd& qinv_matrix) const {
  //--------------------------------------------------------------
  //     sets inviscid ue from panel inviscid tangential velocity
  //--------------------------------------------------------------
  SidePair<Matrix2Xd> inviscid_matrix;
  for (int side = 1; side <= 2; ++side) {
    inviscid_matrix.get(side) =
        Matrix2Xd::Zero(2, lattice.get(side).stationCount);
    for (int stationIndex = 0; stationIndex < lattice.get(side).stationCount - 1; ++stationIndex) {
      const int panelIndex = lattice.get(side).stationToPanel[stationIndex];
      inviscid_matrix.get(side)(0, stationIndex) =
          lattice.get(side).panelInfluenceFactor[stationIndex] * qinv_matrix(0, panelIndex);
      inviscid_matrix.get(side)(1, stationIndex) =
          lattice.get(side).panelInfluenceFactor[stationIndex] * qinv_matrix(1, panelIndex);
    }
  }

  return inviscid_matrix;
}

bool BoundaryLayerWorkflow::tesys(const BoundaryLayerSideProfiles& top_profiles,
                                  const BoundaryLayerSideProfiles& bottom_profiles,
                                  const Edge& edge) {
  blc.clear();

  state.station2 = this->blvar(state.station2, FlowRegimeEnum::Wake);

  const int top_te = lattice.top.trailingEdgeIndex;
  const int bottom_te = lattice.bottom.trailingEdgeIndex;
  const double tte =
      top_profiles.momentumThickness[top_te] +
      bottom_profiles.momentumThickness[bottom_te];
  const double dte =
      top_profiles.displacementThickness[top_te] +
      bottom_profiles.displacementThickness[bottom_te] + edge.ante;
  const double cte =
      (top_profiles.skinFrictionCoeff[top_te] *
           top_profiles.momentumThickness[top_te] +
       bottom_profiles.skinFrictionCoeff[bottom_te] *
           bottom_profiles.momentumThickness[bottom_te]) /
      tte;

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
SetblInputView SetblInputView::fromXFoil(const XFoil& xfoil) {
  const auto& lattice = xfoil.boundaryLayerWorkflow.lattice;
  return SetblInputView{
      xfoil.lblini,
      {lattice.top.profiles.edgeVelocity, lattice.bottom.profiles.edgeVelocity},
      {lattice.top.profiles.skinFrictionCoeff, lattice.bottom.profiles.skinFrictionCoeff},
      {lattice.top.profiles.momentumThickness, lattice.bottom.profiles.momentumThickness},
      {lattice.top.profiles.displacementThickness, lattice.bottom.profiles.displacementThickness},
      {lattice.top.profiles.massFlux, lattice.bottom.profiles.massFlux},
      {lattice.top.skinFrictionCoeffHistory, lattice.bottom.skinFrictionCoeffHistory},
      {lattice.top.transitionIndex, lattice.bottom.transitionIndex}};
}

SetblOutputView SetblOutputView::fromXFoil(XFoil& xfoil) {
  auto& lattice = xfoil.boundaryLayerWorkflow.lattice;
  return SetblOutputView{
      xfoil.lblini,
      xfoil.blCompressibility,
      xfoil.blReynolds,
      xfoil.blTransition,
      {lattice.top.profiles.edgeVelocity, lattice.bottom.profiles.edgeVelocity},
      {lattice.top.profiles.skinFrictionCoeff, lattice.bottom.profiles.skinFrictionCoeff},
      {lattice.top.profiles.momentumThickness, lattice.bottom.profiles.momentumThickness},
      {lattice.top.profiles.displacementThickness, lattice.bottom.profiles.displacementThickness},
      {lattice.top.profiles.massFlux, lattice.bottom.profiles.massFlux},
      {lattice.top.skinFrictionCoeffHistory, lattice.bottom.skinFrictionCoeffHistory},
      {lattice.top.transitionIndex, lattice.bottom.transitionIndex},
      xfoil.va,
      xfoil.vb,
      xfoil.vdel,
      xfoil.vm,
      xfoil.vz,
      xfoil.flowRegime,
      };
}

XFoil::MixedModeStationContext XFoil::prepareMixedModeStation(int side, int ibl,
                                                              int itrold,
                                                              double& ami) {
  MixedModeStationContext ctx;

  ctx.simi = (ibl == 0);
  ctx.wake = ibl > boundaryLayerWorkflow.lattice.get(side).trailingEdgeIndex;
  ctx.xsi = boundaryLayerWorkflow.lattice.get(side).arcLengthCoordinates[ibl];
  ctx.uei = boundaryLayerWorkflow.lattice.get(side).profiles.edgeVelocity[ibl];
  ctx.thi = boundaryLayerWorkflow.lattice.get(side).profiles.momentumThickness[ibl];
  ctx.dsi = boundaryLayerWorkflow.lattice.get(side).profiles.displacementThickness[ibl];

  if (ibl < itrold) {
    ami = boundaryLayerWorkflow.lattice.get(side).profiles.skinFrictionCoeff[ibl];
    ctx.cti = 0.03;
  } else {
    ctx.cti = boundaryLayerWorkflow.lattice.get(side).profiles.skinFrictionCoeff[ibl];
    if (ctx.cti <= 0.0) {
      ctx.cti = 0.03;
    }
  }
  ctx.ami = ami;

  if (ctx.wake) {
    int iw = ibl - boundaryLayerWorkflow.lattice.get(side).trailingEdgeIndex;
    ctx.dswaki = boundaryLayerWorkflow.wgap[iw - 1];
  } else {
    ctx.dswaki = 0.0;
  }

  double thickness_limit = (ibl <= boundaryLayerWorkflow.lattice.get(side).trailingEdgeIndex) ? 1.02 : 1.00005;
  ctx.dsi = std::max(ctx.dsi - ctx.dswaki, thickness_limit * ctx.thi) + ctx.dswaki;

  flowRegime = boundaryLayerWorkflow.determineRegimeForStation(side, ibl, ctx.simi, ctx.wake);

  return ctx;
}

void XFoil::checkTransitionIfNeeded(int side, int ibl, bool skipCheck,
                                    int laminarAdvance, double& ami) {
  if (skipCheck || flowRegime == FlowRegimeEnum::Turbulent || flowRegime == FlowRegimeEnum::Wake) {
    return;
  }

  boundaryLayerWorkflow.trchek(*this);
  ami = boundaryLayerWorkflow.state.station2.param.amplz;
  if (flowRegime == FlowRegimeEnum::Transition) {
    boundaryLayerWorkflow.lattice.get(side).transitionIndex = ibl;
  } else {
    boundaryLayerWorkflow.lattice.get(side).transitionIndex = ibl + laminarAdvance;
  }
}

bool XFoil::performMixedModeNewtonIteration(int side, int ibl, int itrold,
                                            MixedModeStationContext& ctx,
                                            double deps, double senswt,
                                            double& sens, double& sennew,
                                            double& ami) {
  bool converged = false;
  double ueref = 0.0;
  double hkref = 0.0;

  for (int itbl = 1; itbl <= 25; ++itbl) {
    {
    blData updatedCurrent =
      boundaryLayerWorkflow.blprv(*this, boundaryLayerWorkflow.state.current(),
                                  ctx.xsi, ami, ctx.cti, ctx.thi, ctx.dsi,
                                  ctx.dswaki, ctx.uei);
      boundaryLayerWorkflow.state.current() = updatedCurrent;
    }
    boundaryLayerWorkflow.blkin(*this, boundaryLayerWorkflow.state);

    checkTransitionIfNeeded(side, ibl, ctx.simi, 1, ami);

    const bool startOfWake =
        boundaryLayerWorkflow.isStartOfWake(*this, side, ibl);
    boundaryLayerWorkflow.updateSystemMatricesForStation(*this, side, ibl,
                                                          ctx);

    if (itbl == 1) {
      boundaryLayerWorkflow.initializeFirstIterationState(
          *this, side, ibl, itrold, ctx, ueref, hkref, ami);
    }

    if (ctx.simi || startOfWake) {
      boundaryLayerWorkflow.configureSimilarityRow(ueref);
    } else {
      const bool resetSensitivity = (itbl <= 5);
      const bool averageSensitivity = (itbl > 5 && itbl <= 15);
      boundaryLayerWorkflow.configureViscousRow(
          hkref, ueref, senswt, resetSensitivity, averageSensitivity,
          sens, sennew);
    }

    if (boundaryLayerWorkflow.applyMixedModeNewtonStep(*this, side, ibl, deps,
                                                       ami, ctx)) {
      converged = true;
      break;
    }
  }

  return converged;
}

void XFoil::handleMixedModeNonConvergence(int side, int ibl,
                                          MixedModeStationContext& ctx,
                                          double& ami) {
  std::stringstream ss;
  ss << "     mrchdu: convergence failed at " << ibl << " ,  side " << side
     << ", res=" << std::setw(4) << std::fixed << std::setprecision(3)
     << ctx.dmax << "\n";
  writeString(ss.str());

  boundaryLayerWorkflow.resetStationKinematicsAfterFailure(
      side, ibl, ctx,
      BoundaryLayerWorkflow::EdgeVelocityFallbackMode::UsePreviousStation);

  {
    blData updatedCurrent =
        boundaryLayerWorkflow.blprv(
            *this, boundaryLayerWorkflow.state.current(), ctx.xsi, ami,
            ctx.cti, ctx.thi, ctx.dsi, ctx.dswaki, ctx.uei);
    boundaryLayerWorkflow.state.current() = updatedCurrent;
  }
  boundaryLayerWorkflow.blkin(*this, boundaryLayerWorkflow.state);

  checkTransitionIfNeeded(side, ibl, ctx.simi, 2, ami);

  boundaryLayerWorkflow.syncStationRegimeStates(side, ibl, ctx.wake, *this);

  ctx.ami = ami;
}

void XFoil::setupBlReferenceParams(SetblOutputView& output, double& re_clmr,
                                   double& msq_clmr) {
  double clmr = 0.0;
  double ma_clmr = 0.0;
  double herat = 0.0;
  double herat_ms = 0.0;
  //---- set the cl used to define mach, reynolds numbers
  if (analysis_state_.controlByAlpha)
    clmr = cl;
  else
    clmr = analysis_state_.clspec;

  //---- set current minf(cl)
  ma_clmr = getActualMach(clmr, analysis_state_.machType);
  re_clmr = getActualReynolds(clmr, analysis_state_.reynoldsType);
  msq_clmr = 2.0 * analysis_state_.currentMach * ma_clmr;

  //---- set compressibility parameter tklam and derivative tk_msq
  const auto compressibility = buildCompressibilityParams();
  tklam = compressibility.karmanTsienFactor;
  tkl_msq = compressibility.karmanTsienFactor_msq;

  //---- set gas constant (= cp/cv)
  output.blCompressibility.gm1bl = gamm1;

  //---- set parameters for compressibility correction
  output.blCompressibility.qinfbl = analysis_state_.qinf;
  output.blCompressibility.tkbl = tklam;
  output.blCompressibility.tkbl_ms = tkl_msq;

  //---- stagnation density and 1/enthalpy
  output.blCompressibility.rstbl =
      pow((1.0 + 0.5 * output.blCompressibility.gm1bl *
                       analysis_state_.currentMach *
                       analysis_state_.currentMach),
          (1.0 / output.blCompressibility.gm1bl));
  output.blCompressibility.rstbl_ms =
      0.5 * output.blCompressibility.rstbl /
      (1.0 + 0.5 * output.blCompressibility.gm1bl *
                 analysis_state_.currentMach *
                 analysis_state_.currentMach);
  output.blCompressibility.hstinv =
      output.blCompressibility.gm1bl *
      MathUtil::pow(analysis_state_.currentMach /
                        output.blCompressibility.qinfbl,
                    2) /
      (1.0 + 0.5 * output.blCompressibility.gm1bl *
                 analysis_state_.currentMach *
                 analysis_state_.currentMach);
  output.blCompressibility.hstinv_ms =
      output.blCompressibility.gm1bl *
          MathUtil::pow(1.0 / output.blCompressibility.qinfbl, 2) /
          (1.0 + 0.5 * output.blCompressibility.gm1bl *
                     analysis_state_.currentMach *
                     analysis_state_.currentMach) -
      0.5 * output.blCompressibility.gm1bl *
          output.blCompressibility.hstinv /
          (1.0 + 0.5 * output.blCompressibility.gm1bl *
                     analysis_state_.currentMach *
                     analysis_state_.currentMach);

  //---- set reynolds number based on freestream density, velocity, viscosity
  herat = 1.0 - 0.5 * output.blCompressibility.qinfbl *
                     output.blCompressibility.qinfbl *
                     output.blCompressibility.hstinv;
  herat_ms = -0.5 * output.blCompressibility.qinfbl *
             output.blCompressibility.qinfbl *
             output.blCompressibility.hstinv_ms;

  output.blReynolds.reybl =
      analysis_state_.currentRe * sqrt(herat * herat * herat) *
      (1.0 + BoundaryLayerWorkflow::kHvrat) /
      (herat + BoundaryLayerWorkflow::kHvrat);
  output.blReynolds.reybl_re = sqrt(herat * herat * herat) *
                               (1.0 + BoundaryLayerWorkflow::kHvrat) /
                               (herat + BoundaryLayerWorkflow::kHvrat);
  output.blReynolds.reybl_ms =
      output.blReynolds.reybl *
      (1.5 / herat - 1.0 / (herat + BoundaryLayerWorkflow::kHvrat)) *
      herat_ms;

  output.blTransition.amcrit = acrit;
}

void XFoil::initializeAndMarchBl(const SetblInputView& input,
                                 SetblOutputView& output) {
  if (!input.lblini) {
    //----- initialize bl by marching with ue (fudge at separation)
    // TRACE(" initializing bl ...\n");
    writeString("   Initializing bl ...\n");

    boundaryLayerWorkflow.mrchue(*this);
    output.lblini = true;
  }

  //---- march bl with current ue and ds to establish transition
  boundaryLayerWorkflow.mrchdu(*this);
}

void XFoil::prepareEdgeVelocityAndSensitivities(
    const SetblInputView& input, SetblOutputView& output,
    SidePair<VectorXd>& usav, int& jvte1, int& jvte2, double& dule1,
    double& dule2, VectorXd& ule1_m, VectorXd& ule2_m, VectorXd& ute1_m,
    VectorXd& ute2_m, double& ule1_a, double& ule2_a) {
  usav.top = input.edgeVelocity.top;
  usav.bottom = input.edgeVelocity.bottom;

  ueset();
  const auto swapped_edge_velocities = swapEdgeVelocities(usav);
  usav = swapped_edge_velocities.swappedUsav;
  output.edgeVelocity.top = swapped_edge_velocities.restoredUedg.top;
  output.edgeVelocity.bottom = swapped_edge_velocities.restoredUedg.bottom;
  jvte1 = boundaryLayerWorkflow.lattice.top
              .stationToSystem
                  [boundaryLayerWorkflow.lattice.top.trailingEdgeIndex];
  jvte2 = boundaryLayerWorkflow.lattice.bottom
              .stationToSystem
                  [boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex];

  dule1 = output.edgeVelocity.top[0] - usav.top[0];
  dule2 = output.edgeVelocity.bottom[0] - usav.bottom[0];

  //---- set le and te ue sensitivities wrt all m values
  const auto le_te_sensitivities = computeLeTeSensitivities(
      boundaryLayerWorkflow.lattice.get(1).stationToPanel[0],
      boundaryLayerWorkflow.lattice.get(2).stationToPanel[0],
      boundaryLayerWorkflow.lattice.get(1).stationToPanel
          [boundaryLayerWorkflow.lattice.top.trailingEdgeIndex],
      boundaryLayerWorkflow.lattice.get(2).stationToPanel
          [boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex]);
  ule1_m = le_te_sensitivities.ule1_m;
  ule2_m = le_te_sensitivities.ule2_m;
  ute1_m = le_te_sensitivities.ute1_m;
  ute2_m = le_te_sensitivities.ute2_m;

  ule1_a = boundaryLayerWorkflow.lattice.get(1).inviscidEdgeVelocityMatrix(1, 0);
  ule2_a = boundaryLayerWorkflow.lattice.get(2).inviscidEdgeVelocityMatrix(1, 0);
}

void XFoil::assembleBlJacobianForStation(
    int is, int iv, int nsys, const VectorXd& d1_m, const VectorXd& u1_m,
    const VectorXd& d2_m, const VectorXd& u2_m, double xi_ule1,
    double xi_ule2, const VectorXd& ule1_m, const VectorXd& ule2_m,
    double ule1_a, double ule2_a, double u1_a, double d1_a, double u2_a,
    double d2_a, double due1, double dds1, double due2, double dds2,
    double dule1, double dule2, double re_clmr, double msq_clmr,
    SetblOutputView& output) {
  for (int jv = 1; jv <= nsys; jv++) {
    output.vm[0][jv][iv] =
        boundaryLayerWorkflow.blc.a1(0, 2) * d1_m[jv] +
        boundaryLayerWorkflow.blc.a1(0, 3) * u1_m[jv] +
        boundaryLayerWorkflow.blc.a2(0, 2) * d2_m[jv] +
        boundaryLayerWorkflow.blc.a2(0, 3) * u2_m[jv] +
        (boundaryLayerWorkflow.blc.a1(0, 4) +
         boundaryLayerWorkflow.blc.a2(0, 4) +
         boundaryLayerWorkflow.blc.d_xi[0]) *
            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
  }

  output.vb[iv](0, 0) = boundaryLayerWorkflow.blc.a1(0, 0);
  output.vb[iv](0, 1) = boundaryLayerWorkflow.blc.a1(0, 1);

  output.va[iv](0, 0) = boundaryLayerWorkflow.blc.a2(0, 0);
  output.va[iv](0, 1) = boundaryLayerWorkflow.blc.a2(0, 1);

  if (analysis_state_.controlByAlpha)
    output.vdel[iv](0, 1) = boundaryLayerWorkflow.blc.d_re[0] * re_clmr +
                            boundaryLayerWorkflow.blc.d_msq[0] * msq_clmr;
  else
    output.vdel[iv](0, 1) =
        (boundaryLayerWorkflow.blc.a1(0, 3) * u1_a +
         boundaryLayerWorkflow.blc.a1(0, 2) * d1_a) +
        (boundaryLayerWorkflow.blc.a2(0, 3) * u2_a +
         boundaryLayerWorkflow.blc.a2(0, 2) * d2_a) +
        (boundaryLayerWorkflow.blc.a1(0, 4) +
         boundaryLayerWorkflow.blc.a2(0, 4) +
         boundaryLayerWorkflow.blc.d_xi[0]) *
            (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

  output.vdel[iv](0, 0) =
      boundaryLayerWorkflow.blc.rhs[0] +
      (boundaryLayerWorkflow.blc.a1(0, 3) * due1 +
       boundaryLayerWorkflow.blc.a1(0, 2) * dds1) +
      (boundaryLayerWorkflow.blc.a2(0, 3) * due2 +
       boundaryLayerWorkflow.blc.a2(0, 2) * dds2) +
      (boundaryLayerWorkflow.blc.a1(0, 4) +
       boundaryLayerWorkflow.blc.a2(0, 4) +
       boundaryLayerWorkflow.blc.d_xi[0]) *
          (xi_ule1 * dule1 + xi_ule2 * dule2);

  for (int jv = 1; jv <= nsys; jv++) {
    output.vm[1][jv][iv] =
        boundaryLayerWorkflow.blc.a1(1, 2) * d1_m[jv] +
        boundaryLayerWorkflow.blc.a1(1, 3) * u1_m[jv] +
        boundaryLayerWorkflow.blc.a2(1, 2) * d2_m[jv] +
        boundaryLayerWorkflow.blc.a2(1, 3) * u2_m[jv] +
        (boundaryLayerWorkflow.blc.a1(1, 4) +
         boundaryLayerWorkflow.blc.a2(1, 4) +
         boundaryLayerWorkflow.blc.d_xi[1]) *
            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
  }
  output.vb[iv](1, 0) = boundaryLayerWorkflow.blc.a1(1, 0);
  output.vb[iv](1, 1) = boundaryLayerWorkflow.blc.a1(1, 1);

  output.va[iv](1, 0) = boundaryLayerWorkflow.blc.a2(1, 0);
  output.va[iv](1, 1) = boundaryLayerWorkflow.blc.a2(1, 1);

  if (analysis_state_.controlByAlpha)
    output.vdel[iv](1, 1) = boundaryLayerWorkflow.blc.d_re[1] * re_clmr +
                            boundaryLayerWorkflow.blc.d_msq[1] * msq_clmr;
  else
    output.vdel[iv](1, 1) =
        (boundaryLayerWorkflow.blc.a1(1, 3) * u1_a +
         boundaryLayerWorkflow.blc.a1(1, 2) * d1_a) +
        (boundaryLayerWorkflow.blc.a2(1, 3) * u2_a +
         boundaryLayerWorkflow.blc.a2(1, 2) * d2_a) +
        (boundaryLayerWorkflow.blc.a1(1, 4) +
         boundaryLayerWorkflow.blc.a2(1, 4) +
         boundaryLayerWorkflow.blc.d_xi[1]) *
            (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

  output.vdel[iv](1, 0) =
      boundaryLayerWorkflow.blc.rhs[1] +
      (boundaryLayerWorkflow.blc.a1(1, 3) * due1 +
       boundaryLayerWorkflow.blc.a1(1, 2) * dds1) +
      (boundaryLayerWorkflow.blc.a2(1, 3) * due2 +
       boundaryLayerWorkflow.blc.a2(1, 2) * dds2) +
      (boundaryLayerWorkflow.blc.a1(1, 4) +
       boundaryLayerWorkflow.blc.a2(1, 4) +
       boundaryLayerWorkflow.blc.d_xi[1]) *
          (xi_ule1 * dule1 + xi_ule2 * dule2);

  // memory overlap problem
  for (int jv = 1; jv <= nsys; jv++) {
    output.vm[2][jv][iv] =
        boundaryLayerWorkflow.blc.a1(2, 2) * d1_m[jv] +
        boundaryLayerWorkflow.blc.a1(2, 3) * u1_m[jv] +
        boundaryLayerWorkflow.blc.a2(2, 2) * d2_m[jv] +
        boundaryLayerWorkflow.blc.a2(2, 3) * u2_m[jv] +
        (boundaryLayerWorkflow.blc.a1(2, 4) +
         boundaryLayerWorkflow.blc.a2(2, 4) +
         boundaryLayerWorkflow.blc.d_xi[2]) *
            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
  }

  output.vb[iv](2, 0) = boundaryLayerWorkflow.blc.a1(2, 0);
  output.vb[iv](2, 1) = boundaryLayerWorkflow.blc.a1(2, 1);

  output.va[iv](2, 0) = boundaryLayerWorkflow.blc.a2(2, 0);
  output.va[iv](2, 1) = boundaryLayerWorkflow.blc.a2(2, 1);

  if (analysis_state_.controlByAlpha)
    output.vdel[iv](2, 1) = boundaryLayerWorkflow.blc.d_re[2] * re_clmr +
                            boundaryLayerWorkflow.blc.d_msq[2] * msq_clmr;
  else
    output.vdel[iv](2, 1) =
        (boundaryLayerWorkflow.blc.a1(2, 3) * u1_a +
         boundaryLayerWorkflow.blc.a1(2, 2) * d1_a) +
        (boundaryLayerWorkflow.blc.a2(2, 3) * u2_a +
         boundaryLayerWorkflow.blc.a2(2, 2) * d2_a) +
        (boundaryLayerWorkflow.blc.a1(2, 4) +
         boundaryLayerWorkflow.blc.a2(2, 4) +
         boundaryLayerWorkflow.blc.d_xi[2]) *
            (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

  output.vdel[iv](2, 0) =
      boundaryLayerWorkflow.blc.rhs[2] +
      (boundaryLayerWorkflow.blc.a1(2, 3) * due1 +
       boundaryLayerWorkflow.blc.a1(2, 2) * dds1) +
      (boundaryLayerWorkflow.blc.a2(2, 3) * due2 +
       boundaryLayerWorkflow.blc.a2(2, 2) * dds2) +
      (boundaryLayerWorkflow.blc.a1(2, 4) +
       boundaryLayerWorkflow.blc.a2(2, 4) +
       boundaryLayerWorkflow.blc.d_xi[2]) *
          (xi_ule1 * dule1 + xi_ule2 * dule2);
}

SetblOutputView XFoil::setbl(const SetblInputView& input,
                                    SetblOutputView output) {
  //-------------------------------------------------
  //	   sets up the bl newton system coefficients for the current bl
  // variables
  //     and the edge velocities received from setup. the local bl system
  //     coefficients are then incorporated into the global newton system.
  //-------------------------------------------------

  std::stringstream ss;
  int jvte1 = 0, jvte2 = 0;
  VectorXd u1_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd u2_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd d1_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd d2_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd ule1_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd ule2_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd ute1_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd ute2_m = VectorXd::Zero(2 * IVX + 1);

  double msq_clmr = 0.0, mdi;
  double re_clmr = 0.0;
  double ule1_a = 0.0, ule2_a = 0.0, u2_a, due2, dds2;
  double xsi, cti = 0.0, uei, thi, dsi, dswaki;
  double d2_a, d2_m2, d2_u2, dte_mte1, dte_ute1, dte_mte2, dte_ute2;
  double tte, cte, dte, dule1 = 0.0, dule2 = 0.0;
  double xi_ule1, xi_ule2;
  double ami = 0.0, tte_tte1 = 0.0, tte_tte2 = 0.0, cte_tte1 = 0.0,
         cte_tte2 = 0.0, cte_cte1 = 0.0, cte_cte2 = 0.0;

  cti = 0.0; // techwinder added, otherwise variable is not initialized

  setupBlReferenceParams(output, re_clmr, msq_clmr);

  initializeAndMarchBl(input, output);

  SidePair<VectorXd> usav;
  prepareEdgeVelocityAndSensitivities(
      input, output, usav, jvte1, jvte2, dule1, dule2, ule1_m, ule2_m, ute1_m,
      ute2_m, ule1_a, ule2_a);

  //*** process each boundary layer side
  for (int is = 1; is <= 2; is++) {
    //---- there is no station "1" at similarity, so zero everything out
    for (int js = 1; js <= 2; ++js) {
      for (int jbl = 0;
           jbl < boundaryLayerWorkflow.lattice.get(js).stationCount - 1; ++jbl) {
        const int jv = boundaryLayerWorkflow.lattice.get(js).stationToSystem[jbl];
        u1_m[jv] = 0.0;
        d1_m[jv] = 0.0;
      }
    }
    double u1_a = 0.0;
    double d1_a = 0.0;

    double due1 = 0.0;
    double dds1 = 0.0;

    //---- set forced transition arc length position
    output.blTransition.xiforc = boundaryLayerWorkflow.xifset(*this, is);

    //**** sweep downstream setting up bl equation linearizations
    for (int ibl = 0; ibl < boundaryLayerWorkflow.lattice.get(is).stationCount - 1; ++ibl) {
      
      int iv = boundaryLayerWorkflow.lattice.get(is).stationToSystem[ibl];

      const bool stationIsSimilarity = (ibl == 0);
      const bool stationIsWake = (ibl > boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex);
      const bool stationIsTransitionCandidate = (ibl == output.itran.get(is));
      const bool stationIsDownstreamOfTransition = (ibl > output.itran.get(is));
      output.flowRegime =
          boundaryLayerWorkflow.determineRegimeForStation(is, ibl,
                                                          stationIsSimilarity,
                                                          stationIsWake);

      //---- set primary variables for current station
      xsi = boundaryLayerWorkflow.lattice.get(is).arcLengthCoordinates[ibl];
      if (ibl < output.itran.get(is))
        ami = output.skinFrictionCoeff.get(is)[ibl];
      else
        cti = output.skinFrictionCoeff.get(is)[ibl];
      uei = output.edgeVelocity.get(is)[ibl];
      thi = output.momentumThickness.get(is)[ibl];
      mdi = output.massFlux.get(is)[ibl];

      dsi = mdi / uei;

      if (stationIsWake) {
        int iw = ibl - boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex;
        dswaki = boundaryLayerWorkflow.wgap[iw - 1];
      } else
        dswaki = 0.0;

      //---- set derivatives of dsi (= d2)
      d2_m2 = 1.0 / uei;
      d2_u2 = -dsi / uei;

      for (int js = 1; js <= 2; js++) {
        for (int jbl = 0; jbl < boundaryLayerWorkflow.lattice.get(js).stationCount - 1; ++jbl) {
          int jv = boundaryLayerWorkflow.lattice.get(js).stationToSystem[jbl];
          u2_m[jv] = -boundaryLayerWorkflow.lattice.get(is).panelInfluenceFactor[ibl] * boundaryLayerWorkflow.lattice.get(js).panelInfluenceFactor[jbl] *
                     aerodynamicCache.dij(boundaryLayerWorkflow.lattice.get(is).stationToPanel[ibl], boundaryLayerWorkflow.lattice.get(js).stationToPanel[jbl]);
          d2_m[jv] = d2_u2 * u2_m[jv];
        }
      }
      d2_m[iv] = d2_m[iv] + d2_m2;

      u2_a = boundaryLayerWorkflow.lattice.get(is).inviscidEdgeVelocityMatrix(1, ibl);
      d2_a = d2_u2 * u2_a;

  //---- "forced" changes due to mismatch between edge velocities and
  // usav=inviscidEdgeVelocityMatrix+dij*output.massFlux
      due2 = output.edgeVelocity.get(is)[ibl] - usav.get(is)[ibl];
      dds2 = d2_u2 * due2;

  {
    blData updatedCurrent =
        boundaryLayerWorkflow.blprv(*this,
                                    boundaryLayerWorkflow.state.current(), xsi,
                                    ami, cti, thi, dsi, dswaki, uei);
    boundaryLayerWorkflow.state.current() = updatedCurrent;
  } // cti
      boundaryLayerWorkflow.blkin(*this, boundaryLayerWorkflow.state);

      //---- check for transition and set xt, etc. if found
      if (stationIsTransitionCandidate) {
        boundaryLayerWorkflow.trchek(*this);
        ami = boundaryLayerWorkflow.state.station2.param.amplz;
      }

      if (stationIsTransitionCandidate && flowRegime != FlowRegimeEnum::Transition) {
        // TRACE("setbl: xtr???  n1=%d n2=%d: \n", ampl1, ampl2);

        ss << "setbl: xtr???  n1=" << boundaryLayerWorkflow.state.station1.param.amplz
           << " n2=" << boundaryLayerWorkflow.state.station2.param.amplz << ":\n";
        writeString(ss.str());
        ss.str("");
      }

      //---- assemble 10x4 linearized system for dskinFrictionCoeff, dth, dds, due, dxi
      //	   at the previous "1" station and the current "2" station

      if (ibl == boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex + 1) {
        //----- define quantities at start of output.wake, adding te base thickness to
        // dstar
        tte = output.momentumThickness.get(1)[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] +
              output.momentumThickness.get(2)[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex];
        dte = output.displacementThickness.get(1)[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] +
              output.displacementThickness.get(2)[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex] + foil.edge.ante;
        cte = (output.skinFrictionCoeff.get(1)[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] *
                   output.momentumThickness.get(1)[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] +
               output.skinFrictionCoeff.get(2)[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex] *
                   output.momentumThickness.get(2)[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex]) /
               tte;
        boundaryLayerWorkflow.tesys(boundaryLayerWorkflow.lattice.top.profiles,
                                    boundaryLayerWorkflow.lattice.bottom.profiles,
                                    foil.edge);

        tte_tte1 = 1.0;
        tte_tte2 = 1.0;
        dte_mte1 = 1.0 / output.edgeVelocity.top[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex];
        dte_ute1 = -output.displacementThickness.get(1)[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] /
                    output.edgeVelocity.top[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex];
        dte_mte2 = 1.0 / output.edgeVelocity.bottom[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex];
        dte_ute2 = -output.displacementThickness.get(2)[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex] /
                    output.edgeVelocity.bottom[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex];
        cte_cte1 = output.momentumThickness.get(1)[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] / tte;
        cte_cte2 = output.momentumThickness.get(2)[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex] / tte;
        cte_tte1 = (output.skinFrictionCoeff.get(1)[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] - cte) / tte;
        cte_tte2 = (output.skinFrictionCoeff.get(2)[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex] - cte) / tte;

        //----- re-define d1 sensitivities wrt m since d1 depends on both te ds
        // values
      for (int js = 1; js <= 2; js++) {
        for (int jbl = 0; jbl < boundaryLayerWorkflow.lattice.get(js).stationCount - 1; ++jbl) {
            int jv = boundaryLayerWorkflow.lattice.get(js).stationToSystem[jbl];
            d1_m[jv] = dte_ute1 * ute1_m[jv] + dte_ute2 * ute2_m[jv];
          }
        }
        d1_m[jvte1] = d1_m[jvte1] + dte_mte1;
        d1_m[jvte2] = d1_m[jvte2] + dte_mte2;

        //----- "forced" changes from  output.edgeVelocity --- usav=inviscidEdgeVelocityMatrix+dij*output.massFlux mismatch
        due1 = 0.0;
        dds1 =
            dte_ute1 *
                (output.edgeVelocity.top[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] -
                 usav.top[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex]) +
            dte_ute2 *
                (output.edgeVelocity.bottom[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex] -
                 usav.bottom[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex]);
      } else {
        boundaryLayerWorkflow.blsys(*this);
      }

      //---- save wall shear and equil. max shear coefficient for plotting
      // output
      output.skinFrictionCoeffHistory.get(is)[ibl] = boundaryLayerWorkflow.state.station2.cqz.scalar;

      //---- set xi sensitivities wrt le ue changes
      if (is == 1) {
        xi_ule1 = stagnation.sst_go;
        xi_ule2 = -stagnation.sst_gp;
      } else {
        xi_ule1 = -stagnation.sst_go;
        xi_ule2 = stagnation.sst_gp;
      }

      //---- stuff bl system coefficients into main jacobian matrix
      assembleBlJacobianForStation(
          is, iv, nsys, d1_m, u1_m, d2_m, u2_m, xi_ule1, xi_ule2, ule1_m,
          ule2_m, ule1_a, ule2_a, u1_a, d1_a, u2_a, d2_a, due1, dds1, due2,
          dds2, dule1, dule2, re_clmr, msq_clmr, output);

      if (ibl == boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex + 1) {
        //----- redefine coefficients for tte, dte, etc
        output.vz[0][0] = boundaryLayerWorkflow.blc.a1(0, 0) * cte_cte1;
        output.vz[0][1] = boundaryLayerWorkflow.blc.a1(0, 0) * cte_tte1 + boundaryLayerWorkflow.blc.a1(0, 1) * tte_tte1;
        output.vb[iv](0, 0) = boundaryLayerWorkflow.blc.a1(0, 0) * cte_cte2;
        output.vb[iv](0, 1) = boundaryLayerWorkflow.blc.a1(0, 0) * cte_tte2 + boundaryLayerWorkflow.blc.a1(0, 1) * tte_tte2;

        output.vz[1][0] = boundaryLayerWorkflow.blc.a1(1, 0) * cte_cte1;
        output.vz[1][1] = boundaryLayerWorkflow.blc.a1(1, 0) * cte_tte1 + boundaryLayerWorkflow.blc.a1(1, 1) * tte_tte1;
        output.vb[iv](1, 0) = boundaryLayerWorkflow.blc.a1(1, 0) * cte_cte2;
        output.vb[iv](1, 1) = boundaryLayerWorkflow.blc.a1(1, 0) * cte_tte2 + boundaryLayerWorkflow.blc.a1(1, 1) * tte_tte2;

        output.vz[2][0] = boundaryLayerWorkflow.blc.a1(2, 0) * cte_cte1;
        output.vz[2][1] = boundaryLayerWorkflow.blc.a1(2, 0) * cte_tte1 + boundaryLayerWorkflow.blc.a1(2, 1) * tte_tte1;
        output.vb[iv](2, 0) = boundaryLayerWorkflow.blc.a1(2, 0) * cte_cte2;
        output.vb[iv](2, 1) = boundaryLayerWorkflow.blc.a1(2, 0) * cte_tte2 + boundaryLayerWorkflow.blc.a1(2, 1) * tte_tte2;
      }

      //---- turbulent intervals will follow if currently at transition interval
      if (flowRegime == FlowRegimeEnum::Transition) {
        //------ save transition location
        output.itran.get(is) = ibl;
        output.flowRegime = FlowRegimeEnum::Turbulent;
      }

      if (ibl == boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex) {
        //----- set "2" variables at te to output.wake correlations for next station

        output.flowRegime = FlowRegimeEnum::Wake;
        boundaryLayerWorkflow.state.station2 =
            boundaryLayerWorkflow.blvar(
                boundaryLayerWorkflow.state.station2,
                FlowRegimeEnum::Wake);
        boundaryLayerWorkflow.blmid(FlowRegimeEnum::Wake);
      }
      u1_m = u2_m;
      d1_m = d2_m;

      u1_a = u2_a;
      d1_a = d2_a;

      due1 = due2;
      dds1 = dds2;

      //---- set bl variables for next station
      boundaryLayerWorkflow.state.stepbl();
    }
  }

  return output;
}

bool XFoil::ueset() {
  //---------------------------------------------------------
  //     sets ue from inviscid ue plus all source influence
  //---------------------------------------------------------
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 0; ibl < boundaryLayerWorkflow.lattice.get(is).stationCount - 1; ++ibl) {
      double dui = 0.0;
      for (int js = 1; js <= 2; js++) {
        for (int jbl = 0; jbl < boundaryLayerWorkflow.lattice.get(js).stationCount - 1; ++jbl) {
          double ue_m = -boundaryLayerWorkflow.lattice.get(is).panelInfluenceFactor[ibl] * boundaryLayerWorkflow.lattice.get(js).panelInfluenceFactor[jbl] *
                        aerodynamicCache.dij(boundaryLayerWorkflow.lattice.get(is).stationToPanel[ibl],
                            boundaryLayerWorkflow.lattice.get(js).stationToPanel[jbl]);
          dui += ue_m * boundaryLayerWorkflow.lattice.get(js).profiles.massFlux[jbl];
        }
      }
      boundaryLayerWorkflow.lattice.get(is).profiles.edgeVelocity[ibl] =
          boundaryLayerWorkflow.lattice.get(is).inviscidEdgeVelocityMatrix(0, ibl) + dui;
    }
  }
  return true;
}


bool BoundaryLayerWorkflow::xicalc(const Foil& foil) {
  //-------------------------------------------------------------
  //     sets bl arc length array on each airfoil side and wake
  //-------------------------------------------------------------

  const auto arc_lengths =
      computeArcLengthCoordinates(foil, stagnationSst, lattice);
  lattice.top.arcLengthCoordinates = arc_lengths.top;
  lattice.bottom.arcLengthCoordinates = arc_lengths.bottom;
  wgap = computeWakeGap(foil, lattice.bottom, arc_lengths.bottom);
  return true;
}

SidePair<VectorXd> BoundaryLayerWorkflow::computeArcLengthCoordinates(
    const Foil& foil, double stagnationSst,
    const SidePair<BoundaryLayerLattice>& lattice) {
  SidePair<VectorXd> arc_lengths;
  arc_lengths.top = lattice.top.arcLengthCoordinates;
  arc_lengths.bottom = lattice.bottom.arcLengthCoordinates;

  for (int ibl = 0; ibl <= lattice.top.trailingEdgeIndex; ++ibl) {
    arc_lengths.top[ibl] =
        stagnationSst - foil.foil_shape.spline_length[lattice.get(1).stationToPanel[ibl]];
  }

  for (int ibl = 0; ibl <= lattice.bottom.trailingEdgeIndex; ++ibl) {
    arc_lengths.bottom[ibl] =
        foil.foil_shape.spline_length[lattice.get(2).stationToPanel[ibl]] -
        stagnationSst;
  }

  // Wake: start from TE, duplicate TE value at first wake station
  arc_lengths.bottom[lattice.bottom.trailingEdgeIndex + 1] =
      arc_lengths.bottom[lattice.bottom.trailingEdgeIndex];
  for (int ibl = lattice.bottom.trailingEdgeIndex + 2; ibl < lattice.bottom.stationCount; ++ibl) {
    arc_lengths.bottom[ibl] = arc_lengths.bottom[ibl - 1] +
                        (foil.wake_shape.points.col(lattice.get(2).stationToPanel[ibl]) -
                         foil.wake_shape.points.col(lattice.get(2).stationToPanel[ibl - 1]))
                            .norm();
  }
  return arc_lengths;
}

VectorXd BoundaryLayerWorkflow::computeWakeGap(
    const Foil& foil, const BoundaryLayerLattice& bottom,
    const VectorXd& bottomArcLengths) {
  //---- trailing edge flap length to te gap ratio
  const double telrat = 2.50;

  //---- set up parameters for te flap cubics

  const int point_count = foil.foil_shape.n;
  const double crosp = MathUtil::cross2(foil.foil_shape.dpoints_ds.col(point_count - 1).normalized(),
                                        foil.foil_shape.dpoints_ds.col(0).normalized());
  double dwdxte = crosp / sqrt(1.0 - crosp * crosp);

  //---- limit cubic to avoid absurd te gap widths
  dwdxte = std::max(dwdxte, -3.0 / telrat);
  dwdxte = std::min(dwdxte, 3.0 / telrat);

  const double aa = 3.0 + telrat * dwdxte;
  const double bb = -2.0 - telrat * dwdxte;

  VectorXd wgap_result = VectorXd::Zero(foil.wake_shape.n);
  if (foil.edge.sharp) {
    return wgap_result;
  }

  else {
    //----- set te flap (wake gap) array (0-based: iw0=0..foil.wake_shape.n-1)
    for (int iw0 = 0; iw0 < foil.wake_shape.n; iw0++) {
      const int te_bot_0b = bottom.trailingEdgeIndex; // 0-based TE for array indexing
      const double zn =
          1.0 - (bottomArcLengths[te_bot_0b + (iw0 + 1)] -
                 bottomArcLengths[te_bot_0b]) /
                    (telrat * foil.edge.ante);
      if (zn >= 0.0)
        wgap_result[iw0] =
            foil.edge.ante * (aa + bb * zn) * zn * zn;
    }
  }
  return wgap_result;
}


/** -----------------------------------------------------
 * 	   sets forced-transition bl coordinate locations.
 * ----------------------------------------------------- */
double BoundaryLayerWorkflow::xifset(const XFoil& xfoil, int is) const {
  std::stringstream ss;
  VectorXd w1 = VectorXd::Zero(xfoil.foil.foil_shape.n);
  double str;

  if (lattice.get(is).transitionLocation >= 1.0) {
    return lattice.get(is).arcLengthCoordinates[lattice.get(is).trailingEdgeIndex];
  }

  Vector2d point_chord = xfoil.foil.edge.point_te - xfoil.foil.edge.point_le;

  //---- calculate chord-based x/c, y/c
  for (int i = 0; i < xfoil.foil.foil_shape.n; i++) {
    w1[i] = (xfoil.foil.foil_shape.points.col(i) - xfoil.foil.edge.point_le).dot(point_chord.normalized());
  }

  VectorXd w3 = spline::splind(w1, xfoil.foil.foil_shape.spline_length);
  if (is == 1) {
    str = xfoil.foil.edge.sle +
          (xfoil.foil.foil_shape.spline_length[0] - xfoil.foil.edge.sle) *
              lattice.top.transitionLocation;
  } else {
    str = xfoil.foil.edge.sle +
          (xfoil.foil.foil_shape.spline_length[xfoil.foil.foil_shape.n - 1] -
           xfoil.foil.edge.sle) *
              lattice.bottom.transitionLocation;
  }
  str = spline::sinvrt(str, lattice.get(is).transitionLocation, w1, w3,
                       xfoil.foil.foil_shape.spline_length,
                       xfoil.foil.foil_shape.n);
  double xiforc = std::min((str - xfoil.stagnation.sst),
                           lattice.get(is).arcLengthCoordinates[lattice.get(is).trailingEdgeIndex]);
  if (xiforc < 0.0) {
    std::cout << " ***  stagnation point is past trip on side " << is << std::endl;
    return lattice.get(is).arcLengthCoordinates[lattice.get(is).trailingEdgeIndex];
  }

  return xiforc;
}
