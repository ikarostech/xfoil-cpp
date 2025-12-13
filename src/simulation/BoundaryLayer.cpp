#include "BoundaryLayer.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

#include "XFoil.h"
#include "domain/coefficient/skin_friction.hpp"
#include "core/boundary_layer_util.hpp"

using BoundaryContext = BoundaryLayerWorkflow::MixedModeStationContext;
using Eigen::Matrix;
using Eigen::Vector;

bool BoundaryLayerWorkflow::saveblData(int icom) {
  if (icom == 1) {
    blsav[icom] = state.station1;
  } else {
    blsav[icom] = state.station2;
  }
  return true;
}

bool BoundaryLayerWorkflow::restoreblData(int icom) {
  if (icom == 1) {
    state.station1 = blsav[icom];
  } else if (icom == 2) {
    state.station2 = blsav[icom];
  }
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

  SkinFrictionCoefficients skinFriction = blmid(xfoil, xfoil.flowRegime);
  current = blvar(current, xfoil.flowRegime);

  if (xfoil.flowRegime == FlowRegimeEnum::Similarity) {
    state.stepbl();
  }

  if (xfoil.flowRegime == FlowRegimeEnum::Transition) {
    trdif(xfoil);
  } else {
    blc = blDiffSolver.solve(xfoil.flowRegime, state, skinFriction, xfoil.amcrit);
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

  saveblData(1);
  saveblData(2);

  //---- weighting factors for linear interpolation to transition point
  blDiff wf;
  wf.scalar = (xfoil.xt.scalar - state.station1.param.xz) /
                  (state.station2.param.xz - state.station1.param.xz);
  
  double wf2 = (xfoil.xt.scalar - state.station1.param.xz) / (state.station2.param.xz - state.station1.param.xz);
  double wf_xt = 1.0 / (state.station2.param.xz - state.station1.param.xz);
  wf.vector = xfoil.xt.vector * wf_xt;
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
  state.station2.param.xz = xfoil.xt.scalar;
  state.station2.param.tz = tt.scalar;
  state.station2.param.dz = dt.scalar;
  state.station2.param.uz = ut.scalar;

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
  blc = blDiffSolver.solve(FlowRegimeEnum::Laminar, state, laminarSkinFriction, xfoil.amcrit);

  //---- the current newton system is in terms of "1" and "t" variables,
  //-    so calculate its equivalent in terms of "1" and "2" variables.
  //-    in other words, convert residual sensitivities wrt "t" variables
  //-    into sensitivities wrt "1" and "2" variables.  the amplification
  //-    equation is unnecessary here, so the k=1 row is left empty.
  blrez = blc.rhs;
  for (int k = 1; k < 3; k++) {
    blm[k] = blc.d_msq[k] + blc.a2(k, 1) * tt.ms() + blc.a2(k, 2) * dt.ms() +
             blc.a2(k, 3) * ut.ms() + blc.a2(k, 4) * xfoil.xt.ms();
    blr[k] = blc.d_re[k] + blc.a2(k, 1) * tt.re() + blc.a2(k, 2) * dt.re() +
             blc.a2(k, 3) * ut.re() + blc.a2(k, 4) * xfoil.xt.re();
    blx[k] = blc.d_xi[k] + blc.a2(k, 1) * tt.xf() + blc.a2(k, 2) * dt.xf() +
             blc.a2(k, 3) * ut.xf() + blc.a2(k, 4) * xfoil.xt.xf();
  }
  const Eigen::Matrix<double, 4, 5> bl1_transform{
    {tt.a(), tt.t1(), tt.d1(), tt.u1(), tt.x1()},
    {dt.a(), dt.t1(), dt.d1(), dt.u1(), dt.x1()},
    {ut.a(), ut.t1(), ut.d1(), ut.u1(), ut.x1()},
    {xfoil.xt.a(), xfoil.xt.t1(), xfoil.xt.d1(), xfoil.xt.u1(), xfoil.xt.x1()}
  };
  bl1.block<2, 5>(1, 0) =
      blc.a1.middleRows<2>(1) + blc.a2.block<2, 4>(1, 1) * bl1_transform;

  const Eigen::Matrix<double, 4, 4> bl2_transform{
  {tt.t2(), tt.d2(), tt.u2(), tt.x2()},
  {dt.t2(), dt.d2(), dt.u2(), dt.x2()},
  {ut.t2(), ut.d2(), ut.u2(), ut.x2()},
  {xfoil.xt.t2(), xfoil.xt.d2(), xfoil.xt.u2(), xfoil.xt.x2()}
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
  restoreblData(2);

  //---- calculate xt-x2 midpoint cfm value
  SkinFrictionCoefficients turbulentSkinFriction =
      blmid(xfoil, FlowRegimeEnum::Turbulent);

  //---- set up newton system for dct, dth, dds, due, dxi  at  xt and x2
  blc = blDiffSolver.solve(FlowRegimeEnum::Turbulent, state, turbulentSkinFriction, xfoil.amcrit);

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
             blc.a1(k, 2) * dt.ms() + blc.a1(k, 3) * ut.ms() + blc.a1(k, 4) * xfoil.xt.ms();
    btr[k] = blc.d_re[k] + blc.a1(k, 0) * st_re + blc.a1(k, 1) * tt.re() +
             blc.a1(k, 2) * dt.re() + blc.a1(k, 3) * ut.re() + blc.a1(k, 4) * xfoil.xt.re();
    btx[k] = blc.d_xi[k] + blc.a1(k, 0) * st_xf + blc.a1(k, 1) * tt.xf() +
             blc.a1(k, 2) * dt.xf() + blc.a1(k, 3) * ut.xf() + blc.a1(k, 4) * xfoil.xt.xf();
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
  restoreblData(1);

  return true;
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
  double z_ax = 0.0, z_a1 = 0.0, z_t1 = 0.0, z_d1 = 0.0, z_u1 = 0.0, z_x1 = 0.0,
         z_a2 = 0.0, z_t2 = 0.0, z_d2 = 0.0, z_u2 = 0.0, z_x2 = 0.0, z_ms = 0.0,
         z_re = 0.0;
  double amplt_a2, wf, wf_a1, wf_a2, wf_xf, wf_x1, wf_x2;
  double xt_a2, dt_a2, tt_a2;
  double ut_a2;
  double daeps = 0.00005;

  amplt_a2 = 0.0;
  xt_a2 = dt_a2 = tt_a2 = 0.0;
  ut_a2 = 0.0;

  //---- save variables and sensitivities at ibl ("2") for future restoration
  saveblData(2);

  //---- calculate average amplification rate ax over x1..x2 interval
  BoundaryLayerUtil::AxResult ax_result =
      BoundaryLayerUtil::axset(state.station1.hkz.scalar, state.station1.param.tz, state.station1.rtz.scalar,
            state.station1.param.amplz, state.station2.hkz.scalar, state.station2.param.tz,
            state.station2.rtz.scalar, state.station2.param.amplz, xfoil.amcrit);

  //---- set initial guess for iterate n2 (ampl2) at x2
  state.station2.param.amplz = state.station1.param.amplz +
                        ax_result.ax * (state.station2.param.xz - state.station1.param.xz);
  //---- solve implicit system for amplification ampl2
  auto iterateAmplification = [&]() -> bool {
    for (int itam = 0; itam < 30; itam++) {
      //---- define weighting factors wf1,wf2 for defining "t" quantities
      if (state.station2.param.amplz <= xfoil.amcrit) {
        //------ there is no transition yet,  "t" is the same as "2"
        amplt = state.station2.param.amplz;
        amplt_a2 = 1.0;
        sfa = 1.0;
        sfa_a1 = 0.0;
        sfa_a2 = 0.0;
      } else {
        //------ there is transition in x1..x2, "t" is set from n1, n2
        amplt = xfoil.amcrit;
        amplt_a2 = 0.0;
        sfa = (amplt - state.station1.param.amplz) /
              (state.station2.param.amplz - state.station1.param.amplz);
        sfa_a1 = (sfa - 1.0) / (state.station2.param.amplz - state.station1.param.amplz);
        sfa_a2 = (-sfa) / (state.station2.param.amplz - state.station1.param.amplz);
      }

      if (xfoil.xiforc < state.station2.param.xz) {
        sfx = (xfoil.xiforc - state.station1.param.xz) /
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
      xfoil.xt.scalar = state.station1.param.xz * (1 - wf) + state.station2.param.xz * wf;
      tt = state.station1.param.tz * (1 - wf) + state.station2.param.tz * wf;
      dt = state.station1.param.dz * (1 - wf) + state.station2.param.dz * wf;
      ut = state.station1.param.uz * (1 - wf) + state.station2.param.uz * wf;

      xt_a2 = (state.station2.param.xz - state.station1.param.xz) * wf_a2;
      tt_a2 = (state.station2.param.tz - state.station1.param.tz) * wf_a2;
      dt_a2 = (state.station2.param.dz - state.station1.param.dz) * wf_a2;
      ut_a2 = (state.station2.param.uz - state.station1.param.uz) * wf_a2;

      //---- temporarily set "2" variables from "t" for blkin
      state.station2.param.xz = xfoil.xt.scalar;
      state.station2.param.tz = tt;
      state.station2.param.dz = dt;
      state.station2.param.uz = ut;

      //---- calculate laminar secondary "t" variables hkt, rtt
      xfoil.blkin(state);

      blData::blVector hkt = state.station2.hkz;
      blData::blVector rtt = state.station2.rtz;

      //---- restore clobbered "2" variables, except for ampl2
      amsave = state.station2.param.amplz;

      restoreblData(2);

      state.station2.param.amplz = amsave;

      //---- calculate amplification rate ax over current x1-xt interval
      ax_result = BoundaryLayerUtil::axset(state.station1.hkz.scalar, state.station1.param.tz,
                        state.station1.rtz.scalar, state.station1.param.amplz, hkt.scalar, tt, rtt.scalar,
                        amplt, xfoil.amcrit);

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

      if ((state.station2.param.amplz > xfoil.amcrit &&
           state.station2.param.amplz + xfoil.rlx * da2 < xfoil.amcrit) ||
          (state.station2.param.amplz < xfoil.amcrit &&
           state.station2.param.amplz + xfoil.rlx * da2 > xfoil.amcrit)) {
        //------ limited newton step so ampl2 doesn't step across amcrit either
        // way
        state.station2.param.amplz = xfoil.amcrit;
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
  xfoil.trfree = (state.station2.param.amplz >= xfoil.amcrit);
  xfoil.trforc = (xfoil.xiforc > state.station1.param.xz) && (xfoil.xiforc <= state.station2.param.xz);

  //---- set transition interval flag
  const bool transitionDetected = (xfoil.trforc || xfoil.trfree);
  xfoil.flowRegime = transitionDetected ? FlowRegimeEnum::Transition
                                   : FlowRegimeEnum::Laminar;

  if (!transitionDetected)
    return false;

  //---- resolve if both forced and free transition
  if (xfoil.trfree && xfoil.trforc) {
    xfoil.trforc = xfoil.xiforc < xfoil.xt.scalar;
    xfoil.trfree = xfoil.xiforc >= xfoil.xt.scalar;
  }

  if (xfoil.trforc) {
    //----- if forced transition, then xt is prescribed,
    //-     no sense calculating the sensitivities, since we know them...
    xfoil.xt.scalar = xfoil.xiforc;
    xfoil.xt.a() = 0.0;
    xfoil.xt.x1() = 0.0;
    xfoil.xt.t1() = 0.0;
    xfoil.xt.d1() = 0.0;
    xfoil.xt.u1() = 0.0;
    xfoil.xt.x2() = 0.0;
    xfoil.xt.t2() = 0.0;
    xfoil.xt.d2() = 0.0;
    xfoil.xt.u2() = 0.0;
    xfoil.xt.ms() = 0.0;
    xfoil.xt.re() = 0.0;
    xfoil.xt.xf() = 1.0;
    return true;
  }

  //---- free transition ... set sensitivities of xt

  xfoil.xt.x1() = (1 - wf);
  tt_t1 = (1 - wf);
  dt_d1 = (1 - wf);
  ut_u1 = (1 - wf);

  xfoil.xt.x2() = wf;
  tt_t2 = wf;
  dt_d2 = wf;
  ut_u2 = wf;

  xfoil.xt.a() = (state.station2.param.xz - state.station1.param.xz) * wf_a1;
  tt_a1 = (state.station2.param.tz - state.station1.param.tz) * wf_a1;
  dt_a1 = (state.station2.param.dz - state.station1.param.dz) * wf_a1;
  ut_a1 = (state.station2.param.uz - state.station1.param.uz) * wf_a1;

  xfoil.xt.x1() += (state.station2.param.xz - state.station1.param.xz) * wf_x1;
  tt_x1 = (state.station2.param.tz - state.station1.param.tz) * wf_x1;
  dt_x1 = (state.station2.param.dz - state.station1.param.dz) * wf_x1;
  ut_x1 = (state.station2.param.uz - state.station1.param.uz) * wf_x1;

  xfoil.xt.x2() += (state.station2.param.xz - state.station1.param.xz) * wf_x2;
  tt_x2 = (state.station2.param.tz - state.station1.param.tz) * wf_x2;
  dt_x2 = (state.station2.param.dz - state.station1.param.dz) * wf_x2;
  ut_x2 = (state.station2.param.uz - state.station1.param.uz) * wf_x2;

  xfoil.xt.xf() = (state.station2.param.xz - state.station1.param.xz) * wf_xf;

  //---- at this point, ax = ax( hk1, t1, rt1, a1, hkt, tt, rtt, at )
  blData::blVector hkt = state.station2.hkz;
  blData::blVector rtt = state.station2.rtz;

  //---- set sensitivities of ax( t1 d1 u1 a1 t2 d2 u2 a2 ms re )
  double ax_t1 = ax_result.ax_hk1 * state.station1.hkz.t() + ax_result.ax_t1 +
                 ax_result.ax_rt1 * state.station1.rtz.t() +
                 (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
                  ax_result.ax_rt2 * rtt.t()) *
                     tt_t1;
  double ax_d1 =
      ax_result.ax_hk1 * state.station1.hkz.d() + (ax_result.ax_hk2 * hkt.d()) * dt_d1;
  double ax_u1 =
      ax_result.ax_hk1 * state.station1.hkz.u() + ax_result.ax_rt1 * state.station1.rtz.u() +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_u1;
  double ax_a1 =
      ax_result.ax_a1 +
      (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
       ax_result.ax_rt2 * rtt.t()) *
          tt_a1 +
      (ax_result.ax_hk2 * hkt.d()) * dt_a1 +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_a1;
  double ax_x1 =
      (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
       ax_result.ax_rt2 * rtt.t()) *
          tt_x1 +
      (ax_result.ax_hk2 * hkt.d()) * dt_x1 +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_x1;

  double ax_t2 = (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
                  ax_result.ax_rt2 * rtt.t()) *
                 tt_t2;
  double ax_d2 = (ax_result.ax_hk2 * hkt.d()) * dt_d2;
  double ax_u2 =
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_u2;
  double ax_a2 =
      ax_result.ax_a2 * amplt_a2 +
      (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
       ax_result.ax_rt2 * rtt.t()) *
          tt_a2 +
      (ax_result.ax_hk2 * hkt.d()) * dt_a2 +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_a2;
  double ax_x2 =
      (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
       ax_result.ax_rt2 * rtt.t()) *
          tt_x2 +
      (ax_result.ax_hk2 * hkt.d()) * dt_x2 +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_x2;

  double ax_ms = ax_result.ax_hk2 * hkt.ms() + ax_result.ax_rt2 * rtt.ms() +
                 ax_result.ax_hk1 * state.station1.hkz.ms() +
                 ax_result.ax_rt1 * state.station1.rtz.ms();
  double ax_re =
      ax_result.ax_rt2 * rtt.re() + ax_result.ax_rt1 * state.station1.rtz.re();

  //---- set sensitivities of residual res
  z_ax = -(state.station2.param.xz - state.station1.param.xz);

  z_a1 = z_ax * ax_a1 - 1.0;
  z_t1 = z_ax * ax_t1;
  z_d1 = z_ax * ax_d1;
  z_u1 = z_ax * ax_u1;
  z_x1 = z_ax * ax_x1 + ax_result.ax;

  z_a2 = z_ax * ax_a2 + 1.0;
  z_t2 = z_ax * ax_t2;
  z_d2 = z_ax * ax_d2;
  z_u2 = z_ax * ax_u2;
  z_x2 = z_ax * ax_x2 - ax_result.ax;

  z_ms = z_ax * ax_ms;
  z_re = z_ax * ax_re;

  //---- set sensitivities of xt, with res being stationary for a2 constraint
  xfoil.xt.a() = xfoil.xt.a() - (xt_a2 / z_a2) * z_a1;
  xfoil.xt.t1() = -(xt_a2 / z_a2) * z_t1;
  xfoil.xt.d1() = -(xt_a2 / z_a2) * z_d1;
  xfoil.xt.u1() = -(xt_a2 / z_a2) * z_u1;
  xfoil.xt.x1() = xfoil.xt.x1() - (xt_a2 / z_a2) * z_x1;
  xfoil.xt.t2() = -(xt_a2 / z_a2) * z_t2;
  xfoil.xt.d2() = -(xt_a2 / z_a2) * z_d2;
  xfoil.xt.u2() = -(xt_a2 / z_a2) * z_u2;
  xfoil.xt.x2() = xfoil.xt.x2() - (xt_a2 / z_a2) * z_x2;
  xfoil.xt.ms() = -(xt_a2 / z_a2) * z_ms;
  xfoil.xt.re() = -(xt_a2 / z_a2) * z_re;
  xfoil.xt.xf() = 0.0;

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
