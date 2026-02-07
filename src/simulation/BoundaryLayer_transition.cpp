#include "simulation/BoundaryLayer_transition.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "XFoil.h"
#include "core/boundary_layer_util.hpp"
#include "infrastructure/logger.hpp"
#include "simulation/BoundaryLayer.hpp"

using Eigen::Matrix;
using Eigen::Vector;

BoundaryLayerTransitionSolver::BoundaryLayerTransitionSolver(
    BoundaryLayerWorkflow& workflow)
    : workflow_(&workflow) {}

BoundaryLayerWorkflow::BoundaryLayerWorkflow()
    : transitionSolver(*this),
      geometry(lattice, wgap, stagnationIndex, stagnationSst) {}

double BoundaryLayerTransitionSolver::computeTransitionLocation(
    double weightingFactor) const {
  const double upstreamLocation = workflow_->state.station1.param.xz;
  const double downstreamLocation = workflow_->state.station2.param.xz;
  return upstreamLocation * (1.0 - weightingFactor) +
         downstreamLocation * weightingFactor;
}

bool BoundaryLayerTransitionSolver::trdif() {
  //-----------------------------------------------
  //     sets up the newton system governing the
  //     transition interval.  equations governing
  //     the  laminar  part  x1 < xi < xt  and
  //     the turbulent part  xt < xi < x2
  //     are simply summed.
  //-----------------------------------------------
  auto& state = workflow_->state;
  auto& blTransition = workflow_->blTransition;
  auto& blc = workflow_->blc;
  auto& xt = workflow_->xt;

  Matrix<double, 4, 5> bl1, bl2, bt1, bt2;
  Vector<double, 4> blrez, blm, blr, blx, btrez, btm, btr, btx;

  double ctr, ctr_hk2;

  boundaryLayerStore.saveblData(state.station1, 1);
  boundaryLayerStore.saveblData(state.station2, 2);

  //---- weighting factors for linear interpolation to transition point
  blDiff wf;
  wf.scalar = (xt.scalar - state.station1.param.xz) /
              (state.station2.param.xz - state.station1.param.xz);

  double wf2 =
      (xt.scalar - state.station1.param.xz) /
      (state.station2.param.xz - state.station1.param.xz);
  double wf_xt =
      1.0 / (state.station2.param.xz - state.station1.param.xz);
  wf.vector = xt.vector * wf_xt;
  wf.x1() +=
      (wf2 - 1.0) / (state.station2.param.xz - state.station1.param.xz);
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
  workflow_->blkin(state);
  state.station2 =
      boundaryLayerVariablesSolver.solve(state.station2, FlowRegimeEnum::Laminar);

  //---- calculate x1-xt midpoint cfm value
  SkinFrictionCoefficients laminarSkinFriction =
      workflow_->blmid(FlowRegimeEnum::Laminar);

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
  const Eigen::Matrix<double, 4, 5> bl1_transform{{tt.a(), tt.t1(), tt.d1(),
                                                   tt.u1(), tt.x1()},
                                                  {dt.a(), dt.t1(), dt.d1(),
                                                   dt.u1(), dt.x1()},
                                                  {ut.a(), ut.t1(), ut.d1(),
                                                   ut.u1(), ut.x1()},
                                                  {xt.a(), xt.t1(), xt.d1(),
                                                   xt.u1(), xt.x1()}};
  bl1.block<2, 5>(1, 0) =
      blc.a1.middleRows<2>(1) + blc.a2.block<2, 4>(1, 1) * bl1_transform;

  const Eigen::Matrix<double, 4, 4> bl2_transform{
      {tt.t2(), tt.d2(), tt.u2(), tt.x2()},
      {dt.t2(), dt.d2(), dt.u2(), dt.x2()},
      {ut.t2(), ut.d2(), ut.u2(), ut.x2()},
      {xt.t2(), xt.d2(), xt.u2(), xt.x2()}};

  bl2.block<2, 1>(1, 0).setZero();
  bl2.block<2, 4>(1, 1) = blc.a2.block<2, 4>(1, 1) * bl2_transform;

  //**** second, set up turbulent part between xt and x2  ****

  //---- calculate equilibrium shear coefficient cqt at transition point
  state.station2 = boundaryLayerVariablesSolver.solve(state.station2,
                                                      FlowRegimeEnum::Turbulent);

  //---- set initial shear coefficient value st at transition point
  //-    ( note that cq2, cq2_t2, etc. are really "cqt", "cqt_tt", etc.)

  ctr = 1.8 * exp(-3.3 / (state.station2.hkz.scalar - 1.0));
  ctr_hk2 = ctr * 3.3 / (state.station2.hkz.scalar - 1.0) /
            (state.station2.hkz.scalar - 1.0);

  double st = ctr * state.station2.cqz.scalar;
  double st_tt = ctr * state.station2.cqz.t() +
                 state.station2.cqz.scalar * ctr_hk2 * state.station2.hkz.t();
  double st_dt = ctr * state.station2.cqz.d() +
                 state.station2.cqz.scalar * ctr_hk2 * state.station2.hkz.d();
  double st_ut = ctr * state.station2.cqz.u() +
                 state.station2.cqz.scalar * ctr_hk2 * state.station2.hkz.u();
  double st_ms = ctr * state.station2.cqz.ms() +
                 state.station2.cqz.scalar * ctr_hk2 * state.station2.hkz.ms();
  double st_re = ctr * state.station2.cqz.re();

  state.station2.param.amplz = 0.0;
  state.station2.param.sz = st;

  //---- recalculate turbulent secondary "t" variables using proper cti
  state.station2 = boundaryLayerVariablesSolver.solve(state.station2,
                                                      FlowRegimeEnum::Turbulent);

  state.stepbl();
  state.station2 = boundaryLayerStore.restoreblData(2);

  //---- calculate xt-x2 midpoint cfm value
  SkinFrictionCoefficients turbulentSkinFriction =
      workflow_->blmid(FlowRegimeEnum::Turbulent);

  //---- set up newton system for dct, dth, dds, due, dxi  at  xt and x2
  blc =
      blDiffSolver.solve(FlowRegimeEnum::Turbulent, state, turbulentSkinFriction,
                         blTransition.amcrit);

  //---- convert sensitivities wrt "t" variables into sensitivities
  //-    wrt "1" and "2" variables as done before for the laminar part
  Eigen::VectorXd common_st = Eigen::Vector3d{st_tt, st_dt, st_ut};
  Eigen::VectorXd st1 =
      Eigen::Matrix<double, 5, 3>{{tt.a(), dt.a(), ut.a()},
                                  {tt.t1(), dt.t1(), ut.t1()},
                                  {tt.d1(), dt.d1(), ut.d1()},
                                  {tt.u1(), dt.u1(), ut.u1()},
                                  {tt.x1(), dt.x1(), ut.x1()}} *
      common_st;
  Eigen::VectorXd st2 =
      Eigen::Matrix<double, 5, 3>{{0, 0, 0},
                                  {tt.t2(), dt.t2(), ut.t2()},
                                  {tt.d2(), dt.d2(), ut.d2()},
                                  {tt.u2(), dt.u2(), ut.u2()},
                                  {tt.x2(), dt.x2(), ut.x2()}} *
      common_st;

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
    btm[k] = blc.d_msq[k] + blc.a1(k, 0) * st_ms +
             blc.a1(k, 1) * tt.ms() + blc.a1(k, 2) * dt.ms() +
             blc.a1(k, 3) * ut.ms() + blc.a1(k, 4) * xt.ms();
    btr[k] = blc.d_re[k] + blc.a1(k, 0) * st_re +
             blc.a1(k, 1) * tt.re() + blc.a1(k, 2) * dt.re() +
             blc.a1(k, 3) * ut.re() + blc.a1(k, 4) * xt.re();
    btx[k] = blc.d_xi[k] + blc.a1(k, 0) * st_xf +
             blc.a1(k, 1) * tt.xf() + blc.a1(k, 2) * dt.xf() +
             blc.a1(k, 3) * ut.xf() + blc.a1(k, 4) * xt.xf();
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

bool BoundaryLayerTransitionSolver::trchek() {
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
  auto& state = workflow_->state;
  auto& blTransition = workflow_->blTransition;
  auto& flowRegime = workflow_->flowRegime;
  auto& xt = workflow_->xt;

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
  BoundaryLayerUtil::AxResult ax_result = BoundaryLayerUtil::axset(
      state.station1.hkz.scalar, state.station1.param.tz,
      state.station1.rtz.scalar, state.station1.param.amplz,
      state.station2.hkz.scalar, state.station2.param.tz,
      state.station2.rtz.scalar, state.station2.param.amplz,
      blTransition.amcrit);

  //---- set initial guess for iterate n2 (ampl2) at x2
  state.station2.param.amplz =
      state.station1.param.amplz +
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
        sfa_a1 = (sfa - 1.0) /
                 (state.station2.param.amplz - state.station1.param.amplz);
        sfa_a2 = (-sfa) /
                 (state.station2.param.amplz - state.station1.param.amplz);
      }

      if (blTransition.xiforc < state.station2.param.xz) {
        sfx = (blTransition.xiforc - state.station1.param.xz) /
              (state.station2.param.xz - state.station1.param.xz);
        sfx_x1 =
            (sfx - 1.0) / (state.station2.param.xz - state.station1.param.xz);
        sfx_x2 =
            (-sfx) / (state.station2.param.xz - state.station1.param.xz);
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
      workflow_->blkin(state);

      blData::blVector hkt = state.station2.hkz;
      blData::blVector rtt = state.station2.rtz;

      //---- restore clobbered "2" variables, except for ampl2
      amsave = state.station2.param.amplz;

      state.station2 = boundaryLayerStore.restoreblData(2);

      state.station2.param.amplz = amsave;

      //---- calculate amplification rate ax over current x1-xt interval
      ax_result = BoundaryLayerUtil::axset(
          state.station1.hkz.scalar, state.station1.param.tz,
          state.station1.rtz.scalar, state.station1.param.amplz, hkt.scalar,
          tt, rtt.scalar, amplt, blTransition.amcrit);

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
          (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) *
              ut_a2 +
          ax_result.ax_a2 * amplt_a2;

      //---- residual for implicit ampl2 definition (amplification equation)
      res = state.station2.param.amplz - state.station1.param.amplz -
            ax_result.ax * (state.station2.param.xz - state.station1.param.xz);
      res_a2 = 1.0 - ax_result.ax_a2 *
                         (state.station2.param.xz - state.station1.param.xz);

      da2 = -res / res_a2;

      double rlx = 1.0;
      dxt = xt_a2 * da2;

      if (rlx *
              fabs(dxt /
                   (state.station2.param.xz - state.station1.param.xz)) >
          0.05) {
        rlx =
            0.05 * fabs((state.station2.param.xz - state.station1.param.xz) /
                        dxt);
      }

      if (rlx * fabs(da2) > 1.0) {
        rlx = 1.0 * fabs(1.0 / da2);
      }

      //---- check if converged
      if (fabs(da2) < daeps) {
        return true;
      }

      if ((state.station2.param.amplz > blTransition.amcrit &&
           state.station2.param.amplz + rlx * da2 <
               blTransition.amcrit) ||
          (state.station2.param.amplz < blTransition.amcrit &&
           state.station2.param.amplz + rlx * da2 >
               blTransition.amcrit)) {
        //------ limited newton step so ampl2 doesn't step across amcrit either
        // way
        state.station2.param.amplz = blTransition.amcrit;
      } else {
        //------ regular newton step
        state.station2.param.amplz = state.station2.param.amplz + rlx * da2;
      }
    }
    return false;
  };

  if (!iterateAmplification()) {
    // TRACE("trchek2 - n2 convergence failed\n");
    Logger::instance().write("trchek2 - n2 convergence failed\n");
  }

  //---- test for free or forced transition
  bool trfree = (state.station2.param.amplz >= blTransition.amcrit);
  bool trforc = (blTransition.xiforc > state.station1.param.xz) &&
                (blTransition.xiforc <= state.station2.param.xz);

  //---- set transition interval flag
  const bool transitionDetected = (trforc || trfree);
  flowRegime =
      transitionDetected ? FlowRegimeEnum::Transition : FlowRegimeEnum::Laminar;

  if (!transitionDetected)
    return false;

  //---- resolve if both forced and free transition
  if (trfree && trforc) {
    trforc = blTransition.xiforc < xt.scalar;
    trfree = blTransition.xiforc >= xt.scalar;
  }

  if (trforc) {
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
      ax_result.ax_hk1 * state.station1.hkz.d() +
      (ax_result.ax_hk2 * hkt.d()) * dt_d1;
  ax.u1() = ax_result.ax_hk1 * state.station1.hkz.u() +
            ax_result.ax_rt1 * state.station1.rtz.u() +
            (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) *
                ut_u1;
  ax.a() = ax_result.ax_a1 +
           (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
            ax_result.ax_rt2 * rtt.t()) *
               tt_a1 +
           (ax_result.ax_hk2 * hkt.d()) * dt_a1 +
           (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) *
               ut_a1;
  ax.x1() =
      (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
       ax_result.ax_rt2 * rtt.t()) *
          tt_x1 +
      (ax_result.ax_hk2 * hkt.d()) * dt_x1 +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_x1;

  ax.t2() =
      (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
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
  ax.re() =
      ax_result.ax_rt2 * rtt.re() + ax_result.ax_rt1 * state.station1.rtz.re();

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
