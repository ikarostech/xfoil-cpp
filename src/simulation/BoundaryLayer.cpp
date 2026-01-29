#include "BoundaryLayer.hpp"

#include <algorithm>
#include <array>
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

bool BoundaryLayerWorkflow::blkin(BoundaryLayerState& state) {
  //----------------------------------------------------------
  //     calculates turbulence-independent secondary "2"
  //     variables from the primary "2" variables.
  //----------------------------------------------------------
  //BlCompressibilityParams& blCompressibility = this->blCompressibility;
  //BlReynoldsParams& blReynolds = this->blReynolds;
  blData& current = state.current();
  //---- set edge mach number ** 2
  current.param.mz =
      current.param.uz * current.param.uz * blCompressibility.hstinv /
      (blCompressibility.gm1bl *
       (1.0 - 0.5 * current.param.uz * current.param.uz * blCompressibility.hstinv));
  double tr2 = 1.0 + 0.5 * blCompressibility.gm1bl * current.param.mz;
  current.param.mz_uz = 2.0 * current.param.mz * tr2 / current.param.uz;
  current.param.mz_ms =
      current.param.uz * current.param.uz * tr2 /
      (blCompressibility.gm1bl *
       (1.0 - 0.5 * current.param.uz * current.param.uz * blCompressibility.hstinv)) *
      blCompressibility.hstinv_ms;

  //---- set edge density (isentropic relation)
  current.param.rz =
      blCompressibility.rstbl *
      pow(tr2, (-1.0 / blCompressibility.gm1bl));
  current.param.rz_uz = -current.param.rz / tr2 * 0.5 * current.param.mz_uz;
  current.param.rz_ms = -current.param.rz / tr2 * 0.5 * current.param.mz_ms +
                        blCompressibility.rstbl_ms *
                            pow(tr2, (-1.0 / blCompressibility.gm1bl));

  //---- set shape parameter
  current.param.hz = current.param.dz / current.param.tz;
  current.param.hz_dz = 1.0 / current.param.tz;
  current.param.hz_tz = -current.param.hz / current.param.tz;

  //---- set edge static/stagnation enthalpy
  double herat =
      1.0 - 0.5 * current.param.uz * current.param.uz * blCompressibility.hstinv;
  double he_u2 = -current.param.uz * blCompressibility.hstinv;
  double he_ms =
      -0.5 * current.param.uz * current.param.uz * blCompressibility.hstinv_ms;
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
       blReynolds.reybl);
  current.rtz.u() = current.rtz.scalar *
                    (1.0 / current.param.uz +
                     current.param.rz_uz / current.param.rz - v2_he * he_u2);
  current.rtz.t() = current.rtz.scalar / current.param.tz;
  current.rtz.ms() =
      current.rtz.scalar * (current.param.rz_ms / current.param.rz +
                            (1 / blReynolds.reybl * blReynolds.reybl_ms -
                             v2_he * he_ms));
  current.rtz.re() =
      current.rtz.scalar * (blReynolds.reybl_re / blReynolds.reybl);

  return true;
}

bool BoundaryLayerWorkflow::isStartOfWake(int side, int stationIndex) {
  return stationIndex == lattice.get(side).trailingEdgeIndex + 1;
}

void BoundaryLayerWorkflow::updateSystemMatricesForStation(
    XFoil& xfoil, int side, int stationIndex, BoundaryContext& ctx) {
  if (isStartOfWake(side, stationIndex)) {
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
      stationIndex < lattice.get(side).profiles.transitionIndex && stationIndex >= previousTransition;
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
        uem_sq * blCompressibility.hstinv /
        (blCompressibility.gm1bl *
         (1.0 - 0.5 * uem_sq * blCompressibility.hstinv));
    const auto hkin_result =
        boundary_layer::hkin(dsm / thm, msq);
    hkref = hkin_result.hk;
  }

  if (stationIndex < previousTransition) {
    if (flowRegime == FlowRegimeEnum::Transition) {
      lattice.get(side).profiles.skinFrictionCoeff[stationIndex] = 0.03;
    }
    if (flowRegime == FlowRegimeEnum::Turbulent || flowRegime == FlowRegimeEnum::Wake) {
      const double prev =
          (stationIndex >= 1) ? lattice.get(side).profiles.skinFrictionCoeff[stationIndex - 1]
                              : lattice.get(side).profiles.skinFrictionCoeff[stationIndex];
      lattice.get(side).profiles.skinFrictionCoeff[stationIndex] = prev;
    }
    if (flowRegime == FlowRegimeEnum::Transition || flowRegime == FlowRegimeEnum::Turbulent || flowRegime == FlowRegimeEnum::Wake) {
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
  if (stationIndex >= lattice.get(side).profiles.transitionIndex) {
    ctx.dmax = std::max(ctx.dmax,
                        std::fabs(blc.rhs[0] / (10.0 * ctx.cti)));
  }

  xfoil.rlx = 1.0;
  if (ctx.dmax > 0.3) {
    xfoil.rlx = 0.3 / ctx.dmax;
  }

  if (stationIndex < lattice.get(side).profiles.transitionIndex) {
    ami += xfoil.rlx * blc.rhs[0];
    ctx.ami = ami;
  }
  if (stationIndex >= lattice.get(side).profiles.transitionIndex) {
    ctx.cti += xfoil.rlx * blc.rhs[0];
  }
  ctx.thi += xfoil.rlx * blc.rhs[1];
  ctx.dsi += xfoil.rlx * blc.rhs[2];
  ctx.uei += xfoil.rlx * blc.rhs[3];

  if (stationIndex >= lattice.get(side).profiles.transitionIndex) {
    ctx.cti = std::clamp(ctx.cti, 0.0000001, 0.30);
  }

  const double hklim =
      (stationIndex <= lattice.get(side).trailingEdgeIndex) ? 1.02 : 1.00005;
  const double uei_sq = ctx.uei * ctx.uei;
  const double msq =
      uei_sq * blCompressibility.hstinv /
      (blCompressibility.gm1bl *
       (1.0 - 0.5 * uei_sq * blCompressibility.hstinv));
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
      uei * (1.0 - blCompressibility.tkbl) /
      (1.0 - blCompressibility.tkbl * (uei / blCompressibility.qinfbl) *
                 (uei / blCompressibility.qinfbl));
  data.param.uz_uei =
      (1.0 + blCompressibility.tkbl *
                (2.0 * data.param.uz * uei / blCompressibility.qinfbl /
                     blCompressibility.qinfbl -
                 1.0)) /
      (1.0 - blCompressibility.tkbl *
                 (uei / blCompressibility.qinfbl) *
                 (uei / blCompressibility.qinfbl));
  data.param.uz_ms =
      (data.param.uz * (uei / blCompressibility.qinfbl) *
           (uei / blCompressibility.qinfbl) -
       uei) *
      blCompressibility.tkbl_ms /
      (1.0 - blCompressibility.tkbl *
                 (uei / blCompressibility.qinfbl) *
                 (uei / blCompressibility.qinfbl));
  return data;
}

bool BoundaryLayerWorkflow::blsys(XFoil& xfoil) {
  blData& previous = state.previous();
  blData& current = state.current();

  SkinFrictionCoefficients skinFriction = blmid(flowRegime);
  current = blvar(current, flowRegime);

  if (flowRegime == FlowRegimeEnum::Similarity) {
    state.stepbl();
  }

  if (flowRegime == FlowRegimeEnum::Transition) {
    trdif(xfoil);
  } else {
    blc = blDiffSolver.solve(flowRegime, state, skinFriction,
                             blTransition.amcrit);
  }

  if (flowRegime == FlowRegimeEnum::Similarity) {
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

  state.station2.param.amplz = blTransition.amcrit;
  state.station2.param.sz = 0.0;

  //---- calculate laminar secondary "t" variables
  blkin(state);
  state.station2 = blvar(state.station2, FlowRegimeEnum::Laminar);

  //---- calculate x1-xt midpoint cfm value
  SkinFrictionCoefficients laminarSkinFriction =
      blmid(FlowRegimeEnum::Laminar);

  //=    at this point, all "2" variables are really "t" variables at xt

  //---- set up newton system for dam, dth, dds, due, dxi  at  x1 and xt
  blc = blDiffSolver.solve(FlowRegimeEnum::Laminar, state, laminarSkinFriction,
                           blTransition.amcrit);

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
                           blTransition.amcrit);

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
  if (ibl < lattice.get(is).profiles.transitionIndex) {
    return state.station1.hkz.scalar +
           0.03 * (state.station2.param.xz - state.station1.param.xz) / state.station1.param.tz;
  }

  if (ibl == lattice.get(is).profiles.transitionIndex) {
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
            blTransition.amcrit);

  //---- set initial guess for iterate n2 (ampl2) at x2
  state.station2.param.amplz = state.station1.param.amplz +
                        ax_result.ax * (state.station2.param.xz - state.station1.param.xz);
  //---- solve implicit system for amplification ampl2
  auto iterateAmplification = [&]() -> bool {
    for (int itam = 0; itam < 30; itam++) {
      //---- define weighting factors wf1,wf2 for defining "t" quantities
      if (state.station2.param.amplz <= blTransition.amcrit) {
        //------ there is no transition yet,  "t" is the same as "2"
        amplt = state.station2.param.amplz;
        amplt_a2 = 1.0;
        sfa = 1.0;
        sfa_a1 = 0.0;
        sfa_a2 = 0.0;
      } else {
        //------ there is transition in x1..x2, "t" is set from n1, n2
        amplt = blTransition.amcrit;
        amplt_a2 = 0.0;
        sfa = (amplt - state.station1.param.amplz) /
              (state.station2.param.amplz - state.station1.param.amplz);
        sfa_a1 = (sfa - 1.0) / (state.station2.param.amplz - state.station1.param.amplz);
        sfa_a2 = (-sfa) / (state.station2.param.amplz - state.station1.param.amplz);
      }

      if (blTransition.xiforc < state.station2.param.xz) {
        sfx = (blTransition.xiforc - state.station1.param.xz) /
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
      blkin(state);

      blData::blVector hkt = state.station2.hkz;
      blData::blVector rtt = state.station2.rtz;

      //---- restore clobbered "2" variables, except for ampl2
      amsave = state.station2.param.amplz;

      state.station2 = boundaryLayerStore.restoreblData(2);

      state.station2.param.amplz = amsave;

      //---- calculate amplification rate ax over current x1-xt interval
      ax_result = BoundaryLayerUtil::axset(state.station1.hkz.scalar, state.station1.param.tz,
                        state.station1.rtz.scalar, state.station1.param.amplz, hkt.scalar, tt, rtt.scalar,
                        amplt, blTransition.amcrit);

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

      if ((state.station2.param.amplz > blTransition.amcrit &&
           state.station2.param.amplz + xfoil.rlx * da2 < blTransition.amcrit) ||
          (state.station2.param.amplz < blTransition.amcrit &&
           state.station2.param.amplz + xfoil.rlx * da2 > blTransition.amcrit)) {
        //------ limited newton step so ampl2 doesn't step across amcrit either
        // way
        state.station2.param.amplz = blTransition.amcrit;
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
  xfoil.trfree = (state.station2.param.amplz >= blTransition.amcrit);
  xfoil.trforc =
      (blTransition.xiforc > state.station1.param.xz) &&
      (blTransition.xiforc <= state.station2.param.xz);

  //---- set transition interval flag
  const bool transitionDetected = (xfoil.trforc || xfoil.trfree);
  flowRegime = transitionDetected ? FlowRegimeEnum::Transition
                                   : FlowRegimeEnum::Laminar;

  if (!transitionDetected)
    return false;

  //---- resolve if both forced and free transition
  if (xfoil.trfree && xfoil.trforc) {
    xfoil.trforc = blTransition.xiforc < xt.scalar;
    xfoil.trfree = blTransition.xiforc >= xt.scalar;
  }

  if (xfoil.trforc) {
    //----- if forced transition, then xt is prescribed,
    //-     no sense calculating the sensitivities, since we know them...
    xt.scalar = blTransition.xiforc;
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
  if (error_message) {
    error_message->clear();
  }
  const int lattice_size = point_count + wake_point_count;
  lattice.top.resize(lattice_size);
  lattice.bottom.resize(lattice_size);

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
  const int system_size = xfoil.nsys + 1;
  xfoil.va.resize(system_size, XFoil::Matrix3x2d::Zero());
  xfoil.vb.resize(system_size, XFoil::Matrix3x2d::Zero());
  xfoil.vdel.resize(system_size, XFoil::Matrix3x2d::Zero());
  xfoil.vm.resize(system_size);

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

      lattice.top.profiles.transitionIndex += delta;
      lattice.bottom.profiles.transitionIndex -= delta;

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

      lattice.top.profiles.transitionIndex -= delta;
      lattice.bottom.profiles.transitionIndex += delta;

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

SetblOutputView SetblOutputView::fromXFoil(const XFoil& xfoil) {
  SetblOutputView output;
  output.lblini = xfoil.lblini;
  output.blCompressibility = xfoil.boundaryLayerWorkflow.blCompressibility;
  output.blReynolds = xfoil.boundaryLayerWorkflow.blReynolds;
  output.profiles.top = xfoil.boundaryLayerWorkflow.lattice.top.profiles;
  output.profiles.bottom = xfoil.boundaryLayerWorkflow.lattice.bottom.profiles;
  output.va = xfoil.va;
  output.vb = xfoil.vb;
  output.vdel = xfoil.vdel;
  output.vm = xfoil.vm;
  output.vz = xfoil.vz;
  output.flowRegime = xfoil.boundaryLayerWorkflow.flowRegime;
  output.blTransition = xfoil.boundaryLayerWorkflow.blTransition;
  return output;
}

void SetblOutputView::applyToXFoil(XFoil& xfoil) {
  xfoil.lblini = lblini;
  xfoil.boundaryLayerWorkflow.blCompressibility = blCompressibility;
  xfoil.boundaryLayerWorkflow.blReynolds = blReynolds;
  xfoil.boundaryLayerWorkflow.lattice.top.profiles = std::move(profiles.top);
  xfoil.boundaryLayerWorkflow.lattice.bottom.profiles =
      std::move(profiles.bottom);
  xfoil.va = std::move(va);
  xfoil.vb = std::move(vb);
  xfoil.vdel = std::move(vdel);
  xfoil.vm = std::move(vm);
  xfoil.vz = vz;
  xfoil.boundaryLayerWorkflow.flowRegime = flowRegime;
  xfoil.boundaryLayerWorkflow.blTransition = blTransition;
}

void XFoil::checkTransitionIfNeeded(int side, int ibl, bool skipCheck,
                                    int laminarAdvance, double& ami) {
  FlowRegimeEnum& flowRegime = boundaryLayerWorkflow.flowRegime;
  if (skipCheck || flowRegime == FlowRegimeEnum::Turbulent || flowRegime == FlowRegimeEnum::Wake) {
    return;
  }

  boundaryLayerWorkflow.trchek(*this);
  ami = boundaryLayerWorkflow.state.station2.param.amplz;
  if (flowRegime == FlowRegimeEnum::Transition) {
    boundaryLayerWorkflow.lattice.get(side).profiles.transitionIndex = ibl;
  } else {
    boundaryLayerWorkflow.lattice.get(side).profiles.transitionIndex = ibl + laminarAdvance;
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
    boundaryLayerWorkflow.blkin(boundaryLayerWorkflow.state);

    checkTransitionIfNeeded(side, ibl, ctx.simi, 1, ami);

    const bool startOfWake =
        boundaryLayerWorkflow.isStartOfWake(side, ibl);
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
  boundaryLayerWorkflow.blkin(boundaryLayerWorkflow.state);

  checkTransitionIfNeeded(side, ibl, ctx.simi, 2, ami);

  boundaryLayerWorkflow.syncStationRegimeStates(side, ibl, ctx.wake, *this);

  ctx.ami = ami;
}

XFoil::BlReferenceParams XFoil::computeBlReferenceParams() const {
  BlReferenceParams params;
  double clmr = 0.0;
  double ma_clmr = 0.0;
  double herat = 0.0;
  double herat_ms = 0.0;
  //---- set the cl used to define mach, reynolds numbers
  if (analysis_state_.controlByAlpha)
    clmr = aero_coeffs_.cl;
  else
    clmr = analysis_state_.clspec;

  //---- set current minf(cl)
  const double cla = std::max(clmr, 0.000001);
  switch (analysis_state_.machType) {
  case MachType::CONSTANT:
    params.currentMach = analysis_state_.referenceMach;
    ma_clmr = 0.0;
    break;
  case MachType::FIXED_LIFT:
    params.currentMach = analysis_state_.referenceMach / std::sqrt(cla);
    ma_clmr = -0.5 * params.currentMach / cla;
    break;
  case MachType::FIXED_LIFT_AND_DYNAMIC_PRESSURE:
    params.currentMach = analysis_state_.referenceMach;
    ma_clmr = 0.0;
    break;
  default:
    params.currentMach = analysis_state_.currentMach;
    ma_clmr = 0.0;
    break;
  }

  switch (analysis_state_.reynoldsType) {
  case ReynoldsType::CONSTANT:
    params.currentRe = analysis_state_.referenceRe;
    params.re_clmr = 0.0;
    break;
  case ReynoldsType::FIXED_LIFT:
    params.currentRe = analysis_state_.referenceRe / std::sqrt(cla);
    params.re_clmr = -0.5 * params.currentRe / cla;
    break;
  case ReynoldsType::FIXED_LIFT_AND_DYNAMIC_PRESSURE:
    params.currentRe = analysis_state_.referenceRe / cla;
    params.re_clmr = -params.currentRe / cla;
    break;
  default:
    params.currentRe = analysis_state_.currentRe;
    params.re_clmr = 0.0;
    break;
  }

  params.msq_clmr = 2.0 * params.currentMach * ma_clmr;

  //---- set compressibility parameter tklam and derivative tk_msq
  const double beta = std::sqrt(1.0 - params.currentMach * params.currentMach);
  const double beta_msq = -0.5 / beta;
  const double karmanTsienFactor =
      (params.currentMach / (1.0 + beta)) * (params.currentMach / (1.0 + beta));
  const double karmanTsienFactor_msq =
      1.0 / ((1.0 + beta) * (1.0 + beta)) -
      2.0 * karmanTsienFactor / (1.0 + beta) * beta_msq;
  params.tklam = karmanTsienFactor;
  params.tkl_msq = karmanTsienFactor_msq;

  //---- set gas constant (= cp/cv)
  params.blCompressibility.gm1bl = gamm1;

  //---- set parameters for compressibility correction
  params.blCompressibility.qinfbl = analysis_state_.qinf;
  params.blCompressibility.tkbl = params.tklam;
  params.blCompressibility.tkbl_ms = params.tkl_msq;

  //---- stagnation density and 1/enthalpy
  params.blCompressibility.rstbl =
      pow((1.0 + 0.5 * params.blCompressibility.gm1bl *
                       params.currentMach * params.currentMach),
          (1.0 / params.blCompressibility.gm1bl));
  params.blCompressibility.rstbl_ms =
      0.5 * params.blCompressibility.rstbl /
      (1.0 + 0.5 * params.blCompressibility.gm1bl *
                 params.currentMach * params.currentMach);
  params.blCompressibility.hstinv =
      params.blCompressibility.gm1bl *
      MathUtil::pow(params.currentMach / params.blCompressibility.qinfbl,
                    2) /
      (1.0 + 0.5 * params.blCompressibility.gm1bl *
                 params.currentMach * params.currentMach);
  params.blCompressibility.hstinv_ms =
      params.blCompressibility.gm1bl *
          MathUtil::pow(1.0 / params.blCompressibility.qinfbl, 2) /
          (1.0 + 0.5 * params.blCompressibility.gm1bl *
                     params.currentMach * params.currentMach) -
      0.5 * params.blCompressibility.gm1bl *
          params.blCompressibility.hstinv /
          (1.0 + 0.5 * params.blCompressibility.gm1bl *
                     params.currentMach * params.currentMach);

  //---- set reynolds number based on freestream density, velocity, viscosity
  herat = 1.0 - 0.5 * params.blCompressibility.qinfbl *
                     params.blCompressibility.qinfbl *
                     params.blCompressibility.hstinv;
  herat_ms = -0.5 * params.blCompressibility.qinfbl *
             params.blCompressibility.qinfbl *
             params.blCompressibility.hstinv_ms;

  params.blReynolds.reybl =
      params.currentRe * sqrt(herat * herat * herat) *
      (1.0 + BoundaryLayerWorkflow::kHvrat) /
      (herat + BoundaryLayerWorkflow::kHvrat);
  params.blReynolds.reybl_re = sqrt(herat * herat * herat) *
                               (1.0 + BoundaryLayerWorkflow::kHvrat) /
                               (herat + BoundaryLayerWorkflow::kHvrat);
  params.blReynolds.reybl_ms =
      params.blReynolds.reybl *
      (1.5 / herat - 1.0 / (herat + BoundaryLayerWorkflow::kHvrat)) *
      herat_ms;

  params.amcrit = acrit;
  return params;
}

XFoil::BlInitializationPlan XFoil::computeBlInitializationPlan(
    bool lblini) const {
  BlInitializationPlan plan;
  if (!lblini) {
    //----- initialize bl by marching with ue (fudge at separation)
    plan.needsInitialization = true;
    plan.message = "   Initializing bl ...\n";
  }
  return plan;
}

XFoil::EdgeVelocitySensitivityResult XFoil::prepareEdgeVelocityAndSensitivities(
    SidePairRef<const BoundaryLayerSideProfiles> profiles) const {
  EdgeVelocitySensitivityResult result;
  result.usav.top = profiles.top.edgeVelocity;
  result.usav.bottom = profiles.bottom.edgeVelocity;

  result.edgeVelocity = boundaryLayerWorkflow.ueset(aerodynamicCache.dij);
  result.outputEdgeVelocity = result.edgeVelocity;
  for (int is = 1; is <= 2; ++is) {
    for (int ibl = 0;
         ibl < boundaryLayerWorkflow.lattice.get(is).stationCount - 1; ++ibl) {
      result.usav.get(is)[ibl] = result.edgeVelocity.get(is)[ibl];
      result.outputEdgeVelocity.get(is)[ibl] = profiles.get(is).edgeVelocity[ibl];
    }
  }
  result.jvte.top = boundaryLayerWorkflow.lattice.top
                        .stationToSystem
                            [boundaryLayerWorkflow.lattice.top.trailingEdgeIndex];
  result.jvte.bottom = boundaryLayerWorkflow.lattice.bottom
                           .stationToSystem
                               [boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex];

  result.dule.top = result.outputEdgeVelocity.top[0] - result.usav.top[0];
  result.dule.bottom =
      result.outputEdgeVelocity.bottom[0] - result.usav.bottom[0];

  //---- set le and te ue sensitivities wrt all m values
  const auto le_te_sensitivities = computeLeTeSensitivities(
      boundaryLayerWorkflow.lattice.get(1).stationToPanel[0],
      boundaryLayerWorkflow.lattice.get(2).stationToPanel[0],
      boundaryLayerWorkflow.lattice.get(1).stationToPanel
          [boundaryLayerWorkflow.lattice.top.trailingEdgeIndex],
      boundaryLayerWorkflow.lattice.get(2).stationToPanel
          [boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex]);
  result.ule_m.top = le_te_sensitivities.ule1_m;
  result.ule_m.bottom = le_te_sensitivities.ule2_m;
  result.ute_m.top = le_te_sensitivities.ute1_m;
  result.ute_m.bottom = le_te_sensitivities.ute2_m;

  result.ule_a.top =
      boundaryLayerWorkflow.lattice.get(1).inviscidEdgeVelocityMatrix(1, 0);
  result.ule_a.bottom =
      boundaryLayerWorkflow.lattice.get(2).inviscidEdgeVelocityMatrix(1, 0);
  return result;
}

void XFoil::assembleBlJacobianForStation(
    int is, int iv, int nsys, const SidePairRef<const VectorXd>& d_m,
    const SidePairRef<const VectorXd>& u_m,
    const SidePairRef<const double>& xi_ule,
    const SidePairRef<const VectorXd>& ule_m,
    const SidePairRef<const double>& ule_a,
    const SidePairRef<const double>& u_a,
    const SidePairRef<const double>& d_a,
    const SidePairRef<const double>& due,
    const SidePairRef<const double>& dds,
    const SidePairRef<const double>& dule, double re_clmr, double msq_clmr,
    SetblOutputView& output) {
  for (int jv = 1; jv <= nsys; jv++) {
    output.vm.at(0, jv, iv) =
        boundaryLayerWorkflow.blc.a1(0, 2) * d_m.get(1)[jv] +
        boundaryLayerWorkflow.blc.a1(0, 3) * u_m.get(1)[jv] +
        boundaryLayerWorkflow.blc.a2(0, 2) * d_m.get(2)[jv] +
        boundaryLayerWorkflow.blc.a2(0, 3) * u_m.get(2)[jv] +
        (boundaryLayerWorkflow.blc.a1(0, 4) +
         boundaryLayerWorkflow.blc.a2(0, 4) +
         boundaryLayerWorkflow.blc.d_xi[0]) *
            (xi_ule.get(1) * ule_m.get(1)[jv] +
             xi_ule.get(2) * ule_m.get(2)[jv]);
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
        (boundaryLayerWorkflow.blc.a1(0, 3) * u_a.get(1) +
         boundaryLayerWorkflow.blc.a1(0, 2) * d_a.get(1)) +
        (boundaryLayerWorkflow.blc.a2(0, 3) * u_a.get(2) +
         boundaryLayerWorkflow.blc.a2(0, 2) * d_a.get(2)) +
        (boundaryLayerWorkflow.blc.a1(0, 4) +
         boundaryLayerWorkflow.blc.a2(0, 4) +
         boundaryLayerWorkflow.blc.d_xi[0]) *
            (xi_ule.get(1) * ule_a.get(1) +
             xi_ule.get(2) * ule_a.get(2));

  output.vdel[iv](0, 0) =
      boundaryLayerWorkflow.blc.rhs[0] +
      (boundaryLayerWorkflow.blc.a1(0, 3) * due.get(1) +
       boundaryLayerWorkflow.blc.a1(0, 2) * dds.get(1)) +
      (boundaryLayerWorkflow.blc.a2(0, 3) * due.get(2) +
       boundaryLayerWorkflow.blc.a2(0, 2) * dds.get(2)) +
      (boundaryLayerWorkflow.blc.a1(0, 4) +
       boundaryLayerWorkflow.blc.a2(0, 4) +
       boundaryLayerWorkflow.blc.d_xi[0]) *
          (xi_ule.get(1) * dule.get(1) +
           xi_ule.get(2) * dule.get(2));

  for (int jv = 1; jv <= nsys; jv++) {
    output.vm.at(1, jv, iv) =
        boundaryLayerWorkflow.blc.a1(1, 2) * d_m.get(1)[jv] +
        boundaryLayerWorkflow.blc.a1(1, 3) * u_m.get(1)[jv] +
        boundaryLayerWorkflow.blc.a2(1, 2) * d_m.get(2)[jv] +
        boundaryLayerWorkflow.blc.a2(1, 3) * u_m.get(2)[jv] +
        (boundaryLayerWorkflow.blc.a1(1, 4) +
         boundaryLayerWorkflow.blc.a2(1, 4) +
         boundaryLayerWorkflow.blc.d_xi[1]) *
            (xi_ule.get(1) * ule_m.get(1)[jv] +
             xi_ule.get(2) * ule_m.get(2)[jv]);
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
        (boundaryLayerWorkflow.blc.a1(1, 3) * u_a.get(1) +
         boundaryLayerWorkflow.blc.a1(1, 2) * d_a.get(1)) +
        (boundaryLayerWorkflow.blc.a2(1, 3) * u_a.get(2) +
         boundaryLayerWorkflow.blc.a2(1, 2) * d_a.get(2)) +
        (boundaryLayerWorkflow.blc.a1(1, 4) +
         boundaryLayerWorkflow.blc.a2(1, 4) +
         boundaryLayerWorkflow.blc.d_xi[1]) *
            (xi_ule.get(1) * ule_a.get(1) +
             xi_ule.get(2) * ule_a.get(2));

  output.vdel[iv](1, 0) =
      boundaryLayerWorkflow.blc.rhs[1] +
      (boundaryLayerWorkflow.blc.a1(1, 3) * due.get(1) +
       boundaryLayerWorkflow.blc.a1(1, 2) * dds.get(1)) +
      (boundaryLayerWorkflow.blc.a2(1, 3) * due.get(2) +
       boundaryLayerWorkflow.blc.a2(1, 2) * dds.get(2)) +
      (boundaryLayerWorkflow.blc.a1(1, 4) +
       boundaryLayerWorkflow.blc.a2(1, 4) +
       boundaryLayerWorkflow.blc.d_xi[1]) *
          (xi_ule.get(1) * dule.get(1) +
           xi_ule.get(2) * dule.get(2));

  // memory overlap problem
  for (int jv = 1; jv <= nsys; jv++) {
    output.vm.at(2, jv, iv) =
        boundaryLayerWorkflow.blc.a1(2, 2) * d_m.get(1)[jv] +
        boundaryLayerWorkflow.blc.a1(2, 3) * u_m.get(1)[jv] +
        boundaryLayerWorkflow.blc.a2(2, 2) * d_m.get(2)[jv] +
        boundaryLayerWorkflow.blc.a2(2, 3) * u_m.get(2)[jv] +
        (boundaryLayerWorkflow.blc.a1(2, 4) +
         boundaryLayerWorkflow.blc.a2(2, 4) +
         boundaryLayerWorkflow.blc.d_xi[2]) *
            (xi_ule.get(1) * ule_m.get(1)[jv] +
             xi_ule.get(2) * ule_m.get(2)[jv]);
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
        (boundaryLayerWorkflow.blc.a1(2, 3) * u_a.get(1) +
         boundaryLayerWorkflow.blc.a1(2, 2) * d_a.get(1)) +
        (boundaryLayerWorkflow.blc.a2(2, 3) * u_a.get(2) +
         boundaryLayerWorkflow.blc.a2(2, 2) * d_a.get(2)) +
        (boundaryLayerWorkflow.blc.a1(2, 4) +
         boundaryLayerWorkflow.blc.a2(2, 4) +
         boundaryLayerWorkflow.blc.d_xi[2]) *
            (xi_ule.get(1) * ule_a.get(1) +
             xi_ule.get(2) * ule_a.get(2));

  output.vdel[iv](2, 0) =
      boundaryLayerWorkflow.blc.rhs[2] +
      (boundaryLayerWorkflow.blc.a1(2, 3) * due.get(1) +
       boundaryLayerWorkflow.blc.a1(2, 2) * dds.get(1)) +
      (boundaryLayerWorkflow.blc.a2(2, 3) * due.get(2) +
       boundaryLayerWorkflow.blc.a2(2, 2) * dds.get(2)) +
      (boundaryLayerWorkflow.blc.a1(2, 4) +
       boundaryLayerWorkflow.blc.a2(2, 4) +
       boundaryLayerWorkflow.blc.d_xi[2]) *
          (xi_ule.get(1) * dule.get(1) +
           xi_ule.get(2) * dule.get(2));
}

XFoil::SimilarityStationCoefficients
XFoil::resetSimilarityStationCoefficients(const VectorXd& u_m1,
                                          const VectorXd& d_m1) const {
  SimilarityStationCoefficients result;
  result.u_m1 = u_m1;
  result.d_m1 = d_m1;
  for (int js = 1; js <= 2; ++js) {
    for (int jbl = 0; jbl < boundaryLayerWorkflow.lattice.get(js).stationCount - 1;
         ++jbl) {
      const int jv = boundaryLayerWorkflow.lattice.get(js).stationToSystem[jbl];
      result.u_m1[jv] = 0.0;
      result.d_m1[jv] = 0.0;
    }
  }
  return result;
}

XFoil::SideSweepInitResult XFoil::initializeSideSweepState(int is) const {
  SideSweepInitResult result;
  result.u_a1 = 0.0;
  result.d_a1 = 0.0;
  result.due1 = 0.0;
  result.dds1 = 0.0;
  result.xiforc = boundaryLayerWorkflow.xifset(foil, stagnation, is);
  return result;
}

XFoil::StationPrimaryVars XFoil::loadStationPrimaryVars(
    int is, int ibl, bool stationIsWake, const SetblOutputView& output,
    double ami, double cti) const {
  StationPrimaryVars vars;
  vars.xsi = boundaryLayerWorkflow.lattice.get(is).arcLengthCoordinates[ibl];

  if (ibl < output.profiles.get(is).transitionIndex)
    vars.ami = output.profiles.get(is).skinFrictionCoeff[ibl];
  else
    vars.cti = output.profiles.get(is).skinFrictionCoeff[ibl];
  if (ibl < output.profiles.get(is).transitionIndex) {
    vars.cti = cti;
  } else {
    vars.ami = ami;
  }

  vars.uei = output.profiles.get(is).edgeVelocity[ibl];
  vars.thi = output.profiles.get(is).momentumThickness[ibl];
  vars.mdi = output.profiles.get(is).massFlux[ibl];
  vars.dsi = vars.mdi / vars.uei;

  if (stationIsWake) {
    int iw = ibl - boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex;
    vars.dswaki = boundaryLayerWorkflow.wgap[iw - 1];
  } else {
    vars.dswaki = 0.0;
  }

  return vars;
}

XFoil::StationUpdateResult XFoil::updateStationMatricesAndState(
    int is, int ibl, int iv, const StationPrimaryVars& vars,
    const SidePair<VectorXd>& usav, const SetblOutputView& output,
    const BoundaryLayerState& base_state, int system_size) {
  StationUpdateResult result;
  const double d2_m2 = 1.0 / vars.uei;
  const double d2_u2 = -vars.dsi / vars.uei;

  result.u_m2 = VectorXd::Zero(system_size);
  result.d_m2 = VectorXd::Zero(system_size);
  for (int js = 1; js <= 2; js++) {
    for (int jbl = 0;
         jbl < boundaryLayerWorkflow.lattice.get(js).stationCount - 1; ++jbl) {
      int jv = boundaryLayerWorkflow.lattice.get(js).stationToSystem[jbl];
      result.u_m2[jv] =
          -boundaryLayerWorkflow.lattice.get(is).panelInfluenceFactor[ibl] *
          boundaryLayerWorkflow.lattice.get(js).panelInfluenceFactor[jbl] *
          aerodynamicCache.dij(
              boundaryLayerWorkflow.lattice.get(is).stationToPanel[ibl],
              boundaryLayerWorkflow.lattice.get(js).stationToPanel[jbl]);
      result.d_m2[jv] = d2_u2 * result.u_m2[jv];
    }
  }
  result.d_m2[iv] = result.d_m2[iv] + d2_m2;

  result.u_a2 =
      boundaryLayerWorkflow.lattice.get(is).inviscidEdgeVelocityMatrix(
      1, ibl);
  result.d_a2 = d2_u2 * result.u_a2;

  // "forced" changes from mismatch between edge velocities and usav
  result.due2 = output.profiles.get(is).edgeVelocity[ibl] - usav.get(is)[ibl];
  result.dds2 = d2_u2 * result.due2;

  result.state = base_state;
  result.state.current() =
      boundaryLayerWorkflow.blprv(*this, result.state.current(), vars.xsi,
                                  vars.ami, vars.cti, vars.thi, vars.dsi,
                                  vars.dswaki, vars.uei);
  boundaryLayerWorkflow.blkin(result.state);
  return result;
}

XFoil::TransitionLogResult XFoil::buildTransitionLog(
    bool stationIsTransitionCandidate, FlowRegimeEnum flowRegime) const {
  TransitionLogResult result;
  if (stationIsTransitionCandidate &&
      flowRegime != FlowRegimeEnum::Transition) {
    std::stringstream ss;
    ss << "setbl: xtr???  n1="
       << boundaryLayerWorkflow.state.station1.param.amplz
       << " n2=" << boundaryLayerWorkflow.state.station2.param.amplz << ":\n";
    result.message = ss.str();
  }
  return result;
}

XFoil::TeWakeUpdateResult XFoil::computeTeWakeCoefficients(
    int is, int ibl, const SidePair<VectorXd>& usav,
    const SidePair<VectorXd>& ute_m, const SidePair<int>& jvte,
    const VectorXd& d_m1_template, const SetblOutputView& output) const {
  TeWakeUpdateResult result;
  if (ibl != boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex + 1) {
    return result;
  }
  result.isStartOfWake = true;

  result.coeffs.tte =
      output.profiles.get(1).momentumThickness[
          boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] +
      output.profiles.get(2).momentumThickness[
          boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex];
  result.coeffs.dte =
      output.profiles.get(1).displacementThickness[
          boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] +
      output.profiles.get(2).displacementThickness[
          boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex] +
      foil.edge.ante;
  result.coeffs.cte =
      (output.profiles.get(1).skinFrictionCoeff[
           boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] *
           output.profiles.get(1).momentumThickness[
               boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] +
       output.profiles.get(2).skinFrictionCoeff[
           boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex] *
           output.profiles.get(2).momentumThickness[
               boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex]) /
      result.coeffs.tte;

  result.coeffs.tte_tte1 = 1.0;
  result.coeffs.tte_tte2 = 1.0;
  result.coeffs.dte_mte1 =
      1.0 /
      output.profiles.top.edgeVelocity[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex];
  result.coeffs.dte_ute1 =
      -output.profiles.get(1).displacementThickness[
           boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] /
      output.profiles.top.edgeVelocity[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex];
  result.coeffs.dte_mte2 =
      1.0 /
      output.profiles.bottom.edgeVelocity[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex];
  result.coeffs.dte_ute2 =
      -output.profiles.get(2).displacementThickness[
           boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex] /
      output.profiles.bottom.edgeVelocity[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex];
  result.coeffs.cte_cte1 =
      output.profiles.get(1).momentumThickness[
          boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] /
      result.coeffs.tte;
  result.coeffs.cte_cte2 =
      output.profiles.get(2).momentumThickness[
          boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex] /
      result.coeffs.tte;
  result.coeffs.cte_tte1 =
      (output.profiles.get(1).skinFrictionCoeff[
           boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] -
       result.coeffs.cte) /
      result.coeffs.tte;
  result.coeffs.cte_tte2 =
      (output.profiles.get(2).skinFrictionCoeff[
           boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex] -
       result.coeffs.cte) /
      result.coeffs.tte;

  // Re-define d1 sensitivities wrt m since d1 depends on both te ds values.
  result.d_m1 = d_m1_template;
  for (int js = 1; js <= 2; js++) {
    for (int jbl = 0;
         jbl < boundaryLayerWorkflow.lattice.get(js).stationCount - 1; ++jbl) {
      int jv = boundaryLayerWorkflow.lattice.get(js).stationToSystem[jbl];
      result.d_m1[jv] =
          result.coeffs.dte_ute1 * ute_m.get(1)[jv] +
          result.coeffs.dte_ute2 * ute_m.get(2)[jv];
    }
  }
  result.d_m1[jvte.get(1)] =
      result.d_m1[jvte.get(1)] + result.coeffs.dte_mte1;
  result.d_m1[jvte.get(2)] =
      result.d_m1[jvte.get(2)] + result.coeffs.dte_mte2;

  // "forced" changes from edge velocity mismatch
  result.due1 = 0.0;
  result.dds1 =
      result.coeffs.dte_ute1 *
          (output.profiles.top.edgeVelocity[
               boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] -
           usav.top[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex]) +
      result.coeffs.dte_ute2 *
          (output.profiles.bottom.edgeVelocity[
               boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex] -
           usav.bottom[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex]);

  return result;
}

XFoil::TeWakeJacobianAdjustments XFoil::computeTeWakeJacobianAdjustments(
    const TeWakeCoefficients& coeffs) const {
  TeWakeJacobianAdjustments result;
  result.vz[0][0] = boundaryLayerWorkflow.blc.a1(0, 0) * coeffs.cte_cte1;
  result.vz[0][1] = boundaryLayerWorkflow.blc.a1(0, 0) * coeffs.cte_tte1 +
                    boundaryLayerWorkflow.blc.a1(0, 1) * coeffs.tte_tte1;
  result.vb(0, 0) = boundaryLayerWorkflow.blc.a1(0, 0) * coeffs.cte_cte2;
  result.vb(0, 1) = boundaryLayerWorkflow.blc.a1(0, 0) * coeffs.cte_tte2 +
                    boundaryLayerWorkflow.blc.a1(0, 1) * coeffs.tte_tte2;

  result.vz[1][0] = boundaryLayerWorkflow.blc.a1(1, 0) * coeffs.cte_cte1;
  result.vz[1][1] = boundaryLayerWorkflow.blc.a1(1, 0) * coeffs.cte_tte1 +
                    boundaryLayerWorkflow.blc.a1(1, 1) * coeffs.tte_tte1;
  result.vb(1, 0) = boundaryLayerWorkflow.blc.a1(1, 0) * coeffs.cte_cte2;
  result.vb(1, 1) = boundaryLayerWorkflow.blc.a1(1, 0) * coeffs.cte_tte2 +
                    boundaryLayerWorkflow.blc.a1(1, 1) * coeffs.tte_tte2;

  result.vz[2][0] = boundaryLayerWorkflow.blc.a1(2, 0) * coeffs.cte_cte1;
  result.vz[2][1] = boundaryLayerWorkflow.blc.a1(2, 0) * coeffs.cte_tte1 +
                    boundaryLayerWorkflow.blc.a1(2, 1) * coeffs.tte_tte1;
  result.vb(2, 0) = boundaryLayerWorkflow.blc.a1(2, 0) * coeffs.cte_cte2;
  result.vb(2, 1) = boundaryLayerWorkflow.blc.a1(2, 0) * coeffs.cte_tte2 +
                    boundaryLayerWorkflow.blc.a1(2, 1) * coeffs.tte_tte2;
  return result;
}

XFoil::StationArraysAdvanceResult XFoil::advanceStationArrays(
    const VectorXd& u_m2, const VectorXd& d_m2, double u_a2, double d_a2,
    double due2, double dds2) const {
  StationArraysAdvanceResult result;
  result.u_m1 = u_m2;
  result.d_m1 = d_m2;
  result.u_a1 = u_a2;
  result.d_a1 = d_a2;
  result.due1 = due2;
  result.dds1 = dds2;
  return result;
}
namespace {
struct SetblStation {
  VectorXd u_m;
  VectorXd d_m;
  double u_a = 0.0;
  double d_a = 0.0;
  double due = 0.0;
  double dds = 0.0;
  double xi_ule = 0.0;

  void resizeSystem(int system_size) {
    u_m = VectorXd::Zero(system_size);
    d_m = VectorXd::Zero(system_size);
  }
};

struct SetblSideData {
  SidePair<int> jvte{0, 0};
  SidePair<VectorXd> usav;
  SidePair<VectorXd> ule_m;
  SidePair<VectorXd> ute_m;
  SidePair<double> ule_a{0.0, 0.0};
  SidePair<double> dule{0.0, 0.0};

  void resizeSystem(int system_size) {
    usav.top = VectorXd::Zero(system_size);
    usav.bottom = VectorXd::Zero(system_size);
    ule_m.top = VectorXd::Zero(system_size);
    ule_m.bottom = VectorXd::Zero(system_size);
    ute_m.top = VectorXd::Zero(system_size);
    ute_m.bottom = VectorXd::Zero(system_size);
  }
};
}  // namespace

SetblOutputView XFoil::setbl(
    SidePairRef<const BoundaryLayerSideProfiles> profiles) {
  //-------------------------------------------------
  //	   sets up the bl newton system coefficients for the current bl
  // variables
  //     and the edge velocities received from setup. the local bl system
  //     coefficients are then incorporated into the global newton system.
  //-------------------------------------------------
  SetblOutputView output{};
  FlowRegimeEnum& flowRegime = boundaryLayerWorkflow.flowRegime;
  BlCompressibilityParams& blCompressibility = boundaryLayerWorkflow.blCompressibility;
  BlReynoldsParams& blReynolds = boundaryLayerWorkflow.blReynolds;
  BlTransitionParams& blTransition = boundaryLayerWorkflow.blTransition;
  const int system_size = nsys + 1;
  if (output.vm.size < system_size) {
    output.vm.resize(system_size);
  }
  if (static_cast<int>(output.va.size()) < system_size) {
    output.va.resize(system_size, Matrix3x2d::Zero());
  }
  if (static_cast<int>(output.vb.size()) < system_size) {
    output.vb.resize(system_size, Matrix3x2d::Zero());
  }
  if (static_cast<int>(output.vdel.size()) < system_size) {
    output.vdel.resize(system_size, Matrix3x2d::Zero());
  }

  std::array<SetblStation, 2> setblStations{};
  setblStations[0].resizeSystem(system_size);
  setblStations[1].resizeSystem(system_size);

  SetblSideData setblSides;
  setblSides.resizeSystem(system_size);

  double msq_clmr = 0.0;
  double re_clmr = 0.0;
  double cti = 0.0;
  double ami = 0.0;

  cti = 0.0; // techwinder added, otherwise variable is not initialized

  BlReferenceParams reference_params = computeBlReferenceParams();
  analysis_state_.currentMach = reference_params.currentMach;
  analysis_state_.currentRe = reference_params.currentRe;
  re_clmr = reference_params.re_clmr;
  msq_clmr = reference_params.msq_clmr;
  tklam = reference_params.tklam;
  tkl_msq = reference_params.tkl_msq;
  output.blCompressibility = reference_params.blCompressibility;
  output.blReynolds = reference_params.blReynolds;
  output.blTransition.amcrit = reference_params.amcrit;
  blCompressibility = output.blCompressibility;
  blReynolds = output.blReynolds;
  blTransition.amcrit = output.blTransition.amcrit;

  const bool current_lblini = lblini;
  output.lblini = current_lblini;
  BlInitializationPlan bl_init = computeBlInitializationPlan(current_lblini);
  if (bl_init.needsInitialization) {
    //----- initialize bl by marching with ue (fudge at separation)
    writeString(bl_init.message);
    boundaryLayerWorkflow.mrchue(*this);
    output.lblini = true;
    lblini = true;
  }

  //---- march bl with current ue and ds to establish transition
  boundaryLayerWorkflow.mrchdu(*this);
  output.profiles.top = boundaryLayerWorkflow.lattice.top.profiles;
  output.profiles.bottom = boundaryLayerWorkflow.lattice.bottom.profiles;

  EdgeVelocitySensitivityResult edge_result =
      prepareEdgeVelocityAndSensitivities(profiles);
  setblSides.usav = edge_result.usav;
  setblSides.jvte = edge_result.jvte;
  setblSides.dule = edge_result.dule;
  setblSides.ule_m = edge_result.ule_m;
  setblSides.ute_m = edge_result.ute_m;
  setblSides.ule_a = edge_result.ule_a;
  boundaryLayerWorkflow.lattice.top.profiles.edgeVelocity =
      edge_result.edgeVelocity.top;
  boundaryLayerWorkflow.lattice.bottom.profiles.edgeVelocity =
      edge_result.edgeVelocity.bottom;
  output.profiles.top.edgeVelocity = edge_result.outputEdgeVelocity.top;
  output.profiles.bottom.edgeVelocity = edge_result.outputEdgeVelocity.bottom;
  boundaryLayerWorkflow.lattice.top.profiles.edgeVelocity =
      output.profiles.top.edgeVelocity;
  boundaryLayerWorkflow.lattice.bottom.profiles.edgeVelocity =
      output.profiles.bottom.edgeVelocity;

  //*** process each boundary layer side
  for (int is = 1; is <= 2; is++) {
    //---- there is no station "1" at similarity, so zero everything out
    SimilarityStationCoefficients similarity_coeffs =
        resetSimilarityStationCoefficients(setblStations[0].u_m,
                                            setblStations[0].d_m);
    setblStations[0].u_m = similarity_coeffs.u_m1;
    setblStations[0].d_m = similarity_coeffs.d_m1;

    SideSweepInitResult sweep_init = initializeSideSweepState(is);
    setblStations[0].u_a = sweep_init.u_a1;
    setblStations[0].d_a = sweep_init.d_a1;
    setblStations[0].due = sweep_init.due1;
    setblStations[0].dds = sweep_init.dds1;
    output.blTransition.xiforc = sweep_init.xiforc;
    blTransition.xiforc = output.blTransition.xiforc;

    //**** sweep downstream setting up bl equation linearizations
    for (int ibl = 0; ibl < boundaryLayerWorkflow.lattice.get(is).stationCount - 1; ++ibl) {
      
      int iv = boundaryLayerWorkflow.lattice.get(is).stationToSystem[ibl];

      const bool stationIsSimilarity = (ibl == 0);
      const bool stationIsWake = (ibl > boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex);
      const bool stationIsTransitionCandidate = (ibl == output.profiles.get(is).transitionIndex);
      
      output.flowRegime =
          boundaryLayerWorkflow.determineRegimeForStation(is, ibl,
                                                          stationIsSimilarity,
                                                          stationIsWake);
      flowRegime = output.flowRegime;

      //---- set primary variables for current station
      StationPrimaryVars vars = loadStationPrimaryVars(
          is, ibl, stationIsWake, output, ami, cti);
      ami = vars.ami;
      cti = vars.cti;

      StationUpdateResult station_update = updateStationMatricesAndState(
          is, ibl, iv, vars, setblSides.usav, output,
          boundaryLayerWorkflow.state,
          setblStations[1].u_m.size());
      setblStations[1].u_m = station_update.u_m2;
      setblStations[1].d_m = station_update.d_m2;
      setblStations[1].u_a = station_update.u_a2;
      setblStations[1].d_a = station_update.d_a2;
      setblStations[1].due = station_update.due2;
      setblStations[1].dds = station_update.dds2;
      boundaryLayerWorkflow.state = station_update.state;

      //---- check for transition and set xt, etc. if found
      if (stationIsTransitionCandidate) {
        boundaryLayerWorkflow.trchek(*this);
        ami = boundaryLayerWorkflow.state.station2.param.amplz;
      }
      TransitionLogResult transition_log =
          buildTransitionLog(stationIsTransitionCandidate, output.flowRegime);
      if (!transition_log.message.empty()) {
        writeString(transition_log.message);
      }

      //---- assemble 10x4 linearized system for dskinFrictionCoeff, dth, dds, due, dxi
      //	   at the previous "1" station and the current "2" station

      TeWakeUpdateResult te_update = computeTeWakeCoefficients(
          is, ibl, setblSides.usav, setblSides.ute_m, setblSides.jvte,
          setblStations[0].d_m, output);
      if (te_update.isStartOfWake) {
        boundaryLayerWorkflow.tesys(boundaryLayerWorkflow.lattice.top.profiles,
                                    boundaryLayerWorkflow.lattice.bottom.profiles,
                                    foil.edge);
        setblStations[0].d_m = te_update.d_m1;
        setblStations[0].due = te_update.due1;
        setblStations[0].dds = te_update.dds1;
      } else {
        boundaryLayerWorkflow.blsys(*this);
      }

      //---- save wall shear and equil. max shear coefficient for plotting
      // output
      output.profiles.get(is).skinFrictionCoeffHistory[ibl] = boundaryLayerWorkflow.state.station2.cqz.scalar;

      //---- set xi sensitivities wrt le ue changes
      if (is == 1) {
        setblStations[0].xi_ule = stagnation.sst_go;
        setblStations[1].xi_ule = -stagnation.sst_gp;
      } else {
        setblStations[0].xi_ule = -stagnation.sst_go;
        setblStations[1].xi_ule = stagnation.sst_gp;
      }

      //---- stuff bl system coefficients into main jacobian matrix
      assembleBlJacobianForStation(
          is, iv, nsys,
          SidePairRef<const VectorXd>{setblStations[0].d_m, setblStations[1].d_m},
          SidePairRef<const VectorXd>{setblStations[0].u_m, setblStations[1].u_m},
          SidePairRef<const double>{setblStations[0].xi_ule, setblStations[1].xi_ule},
          SidePairRef<const VectorXd>{setblSides.ule_m.get(1), setblSides.ule_m.get(2)},
          SidePairRef<const double>{setblSides.ule_a.get(1), setblSides.ule_a.get(2)},
          SidePairRef<const double>{setblStations[0].u_a, setblStations[1].u_a},
          SidePairRef<const double>{setblStations[0].d_a, setblStations[1].d_a},
          SidePairRef<const double>{setblStations[0].due, setblStations[1].due},
          SidePairRef<const double>{setblStations[0].dds, setblStations[1].dds},
          SidePairRef<const double>{setblSides.dule.get(1), setblSides.dule.get(2)},
          re_clmr, msq_clmr, output);

      if (te_update.isStartOfWake) {
        //----- redefine coefficients for tte, dte, etc
        TeWakeJacobianAdjustments te_jacobian =
            computeTeWakeJacobianAdjustments(te_update.coeffs);
        for (int row = 0; row < 3; ++row) {
          output.vz[row][0] = te_jacobian.vz[row][0];
          output.vz[row][1] = te_jacobian.vz[row][1];
        }
        output.vb[iv] = te_jacobian.vb;
      }

      //---- turbulent intervals will follow if currently at transition interval
      if (output.flowRegime == FlowRegimeEnum::Transition) {
        //------ save transition location
        output.profiles.get(is).transitionIndex = ibl;
        boundaryLayerWorkflow.lattice.get(is).profiles.transitionIndex = ibl;
        output.flowRegime = FlowRegimeEnum::Turbulent;
        flowRegime = output.flowRegime;
      }

      if (ibl == boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex) {
        //----- set "2" variables at te to output.wake correlations for next station
        output.flowRegime = FlowRegimeEnum::Wake;
        flowRegime = output.flowRegime;
        boundaryLayerWorkflow.state.station2 =
            boundaryLayerWorkflow.blvar(
                boundaryLayerWorkflow.state.station2,
                FlowRegimeEnum::Wake);
        boundaryLayerWorkflow.blmid(FlowRegimeEnum::Wake);
      }

      StationArraysAdvanceResult advance = advanceStationArrays(
          setblStations[1].u_m, setblStations[1].d_m,
          setblStations[1].u_a, setblStations[1].d_a,
          setblStations[1].due, setblStations[1].dds);
      setblStations[0].u_m = advance.u_m1;
      setblStations[0].d_m = advance.d_m1;
      setblStations[0].u_a = advance.u_a1;
      setblStations[0].d_a = advance.d_a1;
      setblStations[0].due = advance.due1;
      setblStations[0].dds = advance.dds1;

      //---- set bl variables for next station
      boundaryLayerWorkflow.state.stepbl();
    }
  }

  return output;
}

SidePair<Eigen::VectorXd> BoundaryLayerWorkflow::ueset(
    const Eigen::MatrixXd& dij) const {
  //---------------------------------------------------------
  //     sets ue from inviscid ue plus all source influence
  //---------------------------------------------------------
  SidePair<Eigen::VectorXd> edge_velocity;
  edge_velocity.top = lattice.top.profiles.edgeVelocity;
  edge_velocity.bottom = lattice.bottom.profiles.edgeVelocity;
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 0; ibl < lattice.get(is).stationCount - 1; ++ibl) {
      double dui = 0.0;
      for (int js = 1; js <= 2; js++) {
        for (int jbl = 0; jbl < lattice.get(js).stationCount - 1; ++jbl) {
          double ue_m = -lattice.get(is).panelInfluenceFactor[ibl] * lattice.get(js).panelInfluenceFactor[jbl] *
                        dij(lattice.get(is).stationToPanel[ibl],
                            lattice.get(js).stationToPanel[jbl]);
          dui += ue_m * lattice.get(js).profiles.massFlux[jbl];
        }
      }
      edge_velocity.get(is)[ibl] =
          lattice.get(is).inviscidEdgeVelocityMatrix(0, ibl) + dui;
    }
  }
  return edge_velocity;
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
double BoundaryLayerWorkflow::xifset(const Foil& foil,
                                     const StagnationResult& stagnation,
                                     int is) const {
  std::stringstream ss;
  VectorXd w1 = VectorXd::Zero(foil.foil_shape.n);
  double str;

  if (lattice.get(is).transitionLocation >= 1.0) {
    return lattice.get(is).arcLengthCoordinates[lattice.get(is).trailingEdgeIndex];
  }

  Vector2d point_chord = foil.edge.point_te - foil.edge.point_le;

  //---- calculate chord-based x/c, y/c
  for (int i = 0; i < foil.foil_shape.n; i++) {
    w1[i] =
        (foil.foil_shape.points.col(i) - foil.edge.point_le)
            .dot(point_chord.normalized());
  }

  VectorXd w3 = spline::splind(w1, foil.foil_shape.spline_length);
  if (is == 1) {
    str = foil.edge.sle +
          (foil.foil_shape.spline_length[0] - foil.edge.sle) *
              lattice.top.transitionLocation;
  } else {
    str = foil.edge.sle +
          (foil.foil_shape.spline_length[foil.foil_shape.n - 1] -
           foil.edge.sle) *
              lattice.bottom.transitionLocation;
  }
  str = spline::sinvrt(str, lattice.get(is).transitionLocation, w1, w3,
                       foil.foil_shape.spline_length,
                       foil.foil_shape.n);
  double xiforc = std::min((str - stagnation.sst),
                           lattice.get(is).arcLengthCoordinates[lattice.get(is).trailingEdgeIndex]);
  if (xiforc < 0.0) {
    std::cout << " ***  stagnation point is past trip on side " << is << std::endl;
    return lattice.get(is).arcLengthCoordinates[lattice.get(is).trailingEdgeIndex];
  }

  return xiforc;
}
