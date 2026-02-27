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

struct BoundaryLayerTransitionSolver::TrchekData {
  explicit TrchekData(BoundaryLayerWorkflow& workflow)
      : state(workflow.state),
        blTransition(workflow.blTransition),
        flowRegime(workflow.flowRegime),
        xt(workflow.xt) {
    tt_sens.vector.setZero();
    dt_sens.vector.setZero();
    ut_sens.vector.setZero();
    tt_sens.scalar = 0.0;
    dt_sens.scalar = 0.0;
    ut_sens.scalar = 0.0;
  }

  BoundaryLayerState& state;
  BlTransitionParams& blTransition;
  FlowRegimeEnum& flowRegime;
  blDiff& xt;

  double amplt = 0.0;
  double sfa = 0.0;
  double sfa_a1 = 0.0;
  double sfa_a2 = 0.0;
  double sfx = 0.0;
  double sfx_x1 = 0.0;
  double sfx_x2 = 0.0;
  double sfx_xf = 0.0;
  double tt = 0.0;
  double dt = 0.0;
  double ut = 0.0;
  double amsave = 0.0;
  double res = 0.0;
  double res_a2 = 0.0;
  double da2 = 0.0;
  double dxt = 0.0;
  blDiff tt_sens;
  blDiff dt_sens;
  blDiff ut_sens;
  double amplt_a2 = 0.0;
  double wf = 0.0;
  double wf_a1 = 0.0;
  double wf_a2 = 0.0;
  double wf_xf = 0.0;
  double wf_x1 = 0.0;
  double wf_x2 = 0.0;
  double xt_a2 = 0.0;
  BoundaryLayerUtil::AxResult ax_result{};
};

struct BoundaryLayerTransitionSolver::TrdifData {
  explicit TrdifData(BoundaryLayerWorkflow& workflow)
      : state(workflow.state),
        blTransition(workflow.blTransition),
        blc(workflow.blc),
        xt(workflow.xt) {}

  BoundaryLayerState& state;
  BlTransitionParams& blTransition;
  BlSystemCoeffs& blc;
  blDiff& xt;

  Matrix<double, 4, 5> bl1;
  Matrix<double, 4, 5> bl2;
  Matrix<double, 4, 5> bt1;
  Matrix<double, 4, 5> bt2;
  Vector<double, 4> blrez;
  Vector<double, 4> blm;
  Vector<double, 4> blr;
  Vector<double, 4> blx;
  Vector<double, 4> btrez;
  Vector<double, 4> btm;
  Vector<double, 4> btr;
  Vector<double, 4> btx;

  double ctr = 0.0;
  double ctr_hk2 = 0.0;

  blDiff wf;
  double wf1 = 0.0;
  double wf2 = 0.0;
  double wf_xt = 0.0;
  blDiff tt;
  blDiff dt;
  blDiff ut;

  Eigen::Matrix<double, 4, 5> bl1_transform;
  Eigen::Matrix<double, 4, 4> bl2_transform;

  double st = 0.0;
  double st_tt = 0.0;
  double st_dt = 0.0;
  double st_ut = 0.0;
  double st_ms = 0.0;
  double st_re = 0.0;
  double st_xf = 0.0;
};

double BoundaryLayerTransitionSolver::computeTransitionLocation(
    double weightingFactor) const {
  const double upstreamLocation = workflow_->state.station1.param.xz;
  const double downstreamLocation = workflow_->state.station2.param.xz;
  return upstreamLocation * (1.0 - weightingFactor) +
         downstreamLocation * weightingFactor;
}

bool BoundaryLayerTransitionSolver::iterateAmplification(TrchekData& data) {
  static constexpr double kAmplificationEps = 0.00005;

  data.ax_result = BoundaryLayerUtil::axset(
      data.state.station1.hkz.scalar, data.state.station1.param.tz,
      data.state.station1.rtz.scalar, data.state.station1.param.amplz,
      data.state.station2.hkz.scalar, data.state.station2.param.tz,
      data.state.station2.rtz.scalar, data.state.station2.param.amplz,
      data.blTransition.amcrit);

  data.state.station2.param.amplz =
      data.state.station1.param.amplz +
      data.ax_result.ax *
          (data.state.station2.param.xz - data.state.station1.param.xz);

  for (int itam = 0; itam < 30; ++itam) {
    if (data.state.station2.param.amplz <= data.blTransition.amcrit) {
      data.amplt = data.state.station2.param.amplz;
      data.amplt_a2 = 1.0;
      data.sfa = 1.0;
      data.sfa_a1 = 0.0;
      data.sfa_a2 = 0.0;
    } else {
      data.amplt = data.blTransition.amcrit;
      data.amplt_a2 = 0.0;
      data.sfa = (data.amplt - data.state.station1.param.amplz) /
                 (data.state.station2.param.amplz -
                  data.state.station1.param.amplz);
      data.sfa_a1 =
          (data.sfa - 1.0) /
          (data.state.station2.param.amplz - data.state.station1.param.amplz);
      data.sfa_a2 =
          (-data.sfa) /
          (data.state.station2.param.amplz - data.state.station1.param.amplz);
    }

    if (data.blTransition.xiforc < data.state.station2.param.xz) {
      data.sfx = (data.blTransition.xiforc - data.state.station1.param.xz) /
                 (data.state.station2.param.xz - data.state.station1.param.xz);
      data.sfx_x1 =
          (data.sfx - 1.0) /
          (data.state.station2.param.xz - data.state.station1.param.xz);
      data.sfx_x2 =
          (-data.sfx) /
          (data.state.station2.param.xz - data.state.station1.param.xz);
      data.sfx_xf =
          1.0 / (data.state.station2.param.xz - data.state.station1.param.xz);
    } else {
      data.sfx = 1.0;
      data.sfx_x1 = 0.0;
      data.sfx_x2 = 0.0;
      data.sfx_xf = 0.0;
    }

    if (data.sfa < data.sfx) {
      data.wf = data.sfa;
      data.wf_a1 = data.sfa_a1;
      data.wf_a2 = data.sfa_a2;
      data.wf_x1 = 0.0;
      data.wf_x2 = 0.0;
      data.wf_xf = 0.0;
    } else {
      data.wf = data.sfx;
      data.wf_a1 = 0.0;
      data.wf_a2 = 0.0;
      data.wf_x1 = data.sfx_x1;
      data.wf_x2 = data.sfx_x2;
      data.wf_xf = data.sfx_xf;
    }

    data.xt.scalar = computeTransitionLocation(data.wf);
    data.tt = data.state.station1.param.tz * (1.0 - data.wf) +
              data.state.station2.param.tz * data.wf;
    data.dt = data.state.station1.param.dz * (1.0 - data.wf) +
              data.state.station2.param.dz * data.wf;
    data.ut = data.state.station1.param.uz * (1.0 - data.wf) +
              data.state.station2.param.uz * data.wf;

    data.xt_a2 =
        (data.state.station2.param.xz - data.state.station1.param.xz) *
        data.wf_a2;
    data.tt_sens.a2() =
        (data.state.station2.param.tz - data.state.station1.param.tz) *
        data.wf_a2;
    data.dt_sens.a2() =
        (data.state.station2.param.dz - data.state.station1.param.dz) *
        data.wf_a2;
    data.ut_sens.a2() =
        (data.state.station2.param.uz - data.state.station1.param.uz) *
        data.wf_a2;

    data.state.station2.param.xz = data.xt.scalar;
    data.state.station2.param.tz = data.tt;
    data.state.station2.param.dz = data.dt;
    data.state.station2.param.uz = data.ut;

    workflow_->blkin(data.state);

    blData::blVector hkt = data.state.station2.hkz;
    blData::blVector rtt = data.state.station2.rtz;

    data.amsave = data.state.station2.param.amplz;
    data.state.station2 = boundaryLayerStore.restoreblData(2);
    data.state.station2.param.amplz = data.amsave;

    data.ax_result = BoundaryLayerUtil::axset(
        data.state.station1.hkz.scalar, data.state.station1.param.tz,
        data.state.station1.rtz.scalar, data.state.station1.param.amplz,
        hkt.scalar, data.tt, rtt.scalar, data.amplt, data.blTransition.amcrit);

    if (data.ax_result.ax <= 0.0) {
      return true;
    }

    data.ax_result.ax_a2 =
        (data.ax_result.ax_hk2 * hkt.t() + data.ax_result.ax_t2 +
         data.ax_result.ax_rt2 * rtt.t()) *
            data.tt_sens.a2() +
        (data.ax_result.ax_hk2 * hkt.d()) * data.dt_sens.a2() +
        (data.ax_result.ax_hk2 * hkt.u() + data.ax_result.ax_rt2 * rtt.u()) *
            data.ut_sens.a2() +
        data.ax_result.ax_a2 * data.amplt_a2;

    data.res = data.state.station2.param.amplz - data.state.station1.param.amplz -
               data.ax_result.ax *
                   (data.state.station2.param.xz - data.state.station1.param.xz);
    data.res_a2 =
        1.0 -
        data.ax_result.ax_a2 *
            (data.state.station2.param.xz - data.state.station1.param.xz);

    data.da2 = -data.res / data.res_a2;

    double rlx = 1.0;
    data.dxt = data.xt_a2 * data.da2;

    if (rlx * std::fabs(data.dxt /
                        (data.state.station2.param.xz -
                         data.state.station1.param.xz)) >
        0.05) {
      rlx =
          0.05 *
          std::fabs((data.state.station2.param.xz - data.state.station1.param.xz) /
                    data.dxt);
    }

    if (rlx * std::fabs(data.da2) > 1.0) {
      rlx = std::fabs(1.0 / data.da2);
    }

    if (std::fabs(data.da2) < kAmplificationEps) {
      return true;
    }

    if ((data.state.station2.param.amplz > data.blTransition.amcrit &&
         data.state.station2.param.amplz + rlx * data.da2 <
             data.blTransition.amcrit) ||
        (data.state.station2.param.amplz < data.blTransition.amcrit &&
         data.state.station2.param.amplz + rlx * data.da2 >
             data.blTransition.amcrit)) {
      data.state.station2.param.amplz = data.blTransition.amcrit;
    } else {
      data.state.station2.param.amplz += rlx * data.da2;
    }
  }

  return false;
}

bool BoundaryLayerTransitionSolver::resolveTransitionLocationAndSensitivities(
    TrchekData& data) {
  bool trfree = (data.state.station2.param.amplz >= data.blTransition.amcrit);
  bool trforc = (data.blTransition.xiforc > data.state.station1.param.xz) &&
                (data.blTransition.xiforc <= data.state.station2.param.xz);

  const bool transitionDetected = trforc || trfree;
  data.flowRegime =
      transitionDetected ? FlowRegimeEnum::Transition : FlowRegimeEnum::Laminar;

  if (!transitionDetected) {
    return false;
  }

  if (trfree && trforc) {
    trforc = data.blTransition.xiforc < data.xt.scalar;
  }

  if (trforc) {
    data.xt.scalar = data.blTransition.xiforc;
    data.xt.a1() = 0.0;
    data.xt.x1() = 0.0;
    data.xt.t1() = 0.0;
    data.xt.d1() = 0.0;
    data.xt.u1() = 0.0;
    data.xt.x2() = 0.0;
    data.xt.t2() = 0.0;
    data.xt.d2() = 0.0;
    data.xt.u2() = 0.0;
    data.xt.a2() = 0.0;
    data.xt.ms() = 0.0;
    data.xt.re() = 0.0;
    data.xt.xf() = 1.0;
    return true;
  }

  data.xt.x1() = (1.0 - data.wf);
  data.tt_sens.t1() = (1.0 - data.wf);
  data.dt_sens.d1() = (1.0 - data.wf);
  data.ut_sens.u1() = (1.0 - data.wf);

  data.xt.x2() = data.wf;
  data.tt_sens.t2() = data.wf;
  data.dt_sens.d2() = data.wf;
  data.ut_sens.u2() = data.wf;

  data.xt.a1() = (data.state.station2.param.xz - data.state.station1.param.xz) *
                 data.wf_a1;
  data.tt_sens.a1() =
      (data.state.station2.param.tz - data.state.station1.param.tz) *
      data.wf_a1;
  data.dt_sens.a1() =
      (data.state.station2.param.dz - data.state.station1.param.dz) *
      data.wf_a1;
  data.ut_sens.a1() =
      (data.state.station2.param.uz - data.state.station1.param.uz) *
      data.wf_a1;

  data.xt.x1() += (data.state.station2.param.xz - data.state.station1.param.xz) *
                  data.wf_x1;
  data.tt_sens.x1() =
      (data.state.station2.param.tz - data.state.station1.param.tz) *
      data.wf_x1;
  data.dt_sens.x1() =
      (data.state.station2.param.dz - data.state.station1.param.dz) *
      data.wf_x1;
  data.ut_sens.x1() =
      (data.state.station2.param.uz - data.state.station1.param.uz) *
      data.wf_x1;

  data.xt.x2() += (data.state.station2.param.xz - data.state.station1.param.xz) *
                  data.wf_x2;
  data.tt_sens.x2() =
      (data.state.station2.param.tz - data.state.station1.param.tz) *
      data.wf_x2;
  data.dt_sens.x2() =
      (data.state.station2.param.dz - data.state.station1.param.dz) *
      data.wf_x2;
  data.ut_sens.x2() =
      (data.state.station2.param.uz - data.state.station1.param.uz) *
      data.wf_x2;

  data.xt.a2() = data.xt_a2;
  data.xt.xf() = (data.state.station2.param.xz - data.state.station1.param.xz) *
                 data.wf_xf;

  blData::blVector hkt = data.state.station2.hkz;
  blData::blVector rtt = data.state.station2.rtz;

  blDiff ax;
  ax.vector.setZero();
  ax.scalar = 0.0;

  ax.t1() = data.ax_result.ax_hk1 * data.state.station1.hkz.t() +
            data.ax_result.ax_t1 +
            data.ax_result.ax_rt1 * data.state.station1.rtz.t() +
            (data.ax_result.ax_hk2 * hkt.t() + data.ax_result.ax_t2 +
             data.ax_result.ax_rt2 * rtt.t()) *
                data.tt_sens.t1();
  ax.d1() = data.ax_result.ax_hk1 * data.state.station1.hkz.d() +
            (data.ax_result.ax_hk2 * hkt.d()) * data.dt_sens.d1();
  ax.u1() = data.ax_result.ax_hk1 * data.state.station1.hkz.u() +
            data.ax_result.ax_rt1 * data.state.station1.rtz.u() +
            (data.ax_result.ax_hk2 * hkt.u() + data.ax_result.ax_rt2 * rtt.u()) *
                data.ut_sens.u1();
  ax.a1() = data.ax_result.ax_a1 +
            (data.ax_result.ax_hk2 * hkt.t() + data.ax_result.ax_t2 +
             data.ax_result.ax_rt2 * rtt.t()) *
                data.tt_sens.a1() +
            (data.ax_result.ax_hk2 * hkt.d()) * data.dt_sens.a1() +
            (data.ax_result.ax_hk2 * hkt.u() + data.ax_result.ax_rt2 * rtt.u()) *
                data.ut_sens.a1();
  ax.x1() =
      (data.ax_result.ax_hk2 * hkt.t() + data.ax_result.ax_t2 +
       data.ax_result.ax_rt2 * rtt.t()) *
          data.tt_sens.x1() +
      (data.ax_result.ax_hk2 * hkt.d()) * data.dt_sens.x1() +
      (data.ax_result.ax_hk2 * hkt.u() + data.ax_result.ax_rt2 * rtt.u()) *
          data.ut_sens.x1();

  ax.t2() =
      (data.ax_result.ax_hk2 * hkt.t() + data.ax_result.ax_t2 +
       data.ax_result.ax_rt2 * rtt.t()) *
      data.tt_sens.t2();
  ax.d2() = (data.ax_result.ax_hk2 * hkt.d()) * data.dt_sens.d2();
  ax.u2() =
      (data.ax_result.ax_hk2 * hkt.u() + data.ax_result.ax_rt2 * rtt.u()) *
      data.ut_sens.u2();
  ax.a2() = data.ax_result.ax_a2 * data.amplt_a2 +
            (data.ax_result.ax_hk2 * hkt.t() + data.ax_result.ax_t2 +
             data.ax_result.ax_rt2 * rtt.t()) *
                data.tt_sens.a2() +
            (data.ax_result.ax_hk2 * hkt.d()) * data.dt_sens.a2() +
            (data.ax_result.ax_hk2 * hkt.u() + data.ax_result.ax_rt2 * rtt.u()) *
                data.ut_sens.a2();
  ax.x2() =
      (data.ax_result.ax_hk2 * hkt.t() + data.ax_result.ax_t2 +
       data.ax_result.ax_rt2 * rtt.t()) *
          data.tt_sens.x2() +
      (data.ax_result.ax_hk2 * hkt.d()) * data.dt_sens.x2() +
      (data.ax_result.ax_hk2 * hkt.u() + data.ax_result.ax_rt2 * rtt.u()) *
          data.ut_sens.x2();

  ax.ms() = data.ax_result.ax_hk2 * hkt.ms() + data.ax_result.ax_rt2 * rtt.ms() +
            data.ax_result.ax_hk1 * data.state.station1.hkz.ms() +
            data.ax_result.ax_rt1 * data.state.station1.rtz.ms();
  ax.re() =
      data.ax_result.ax_rt2 * rtt.re() + data.ax_result.ax_rt1 * data.state.station1.rtz.re();

  blDiff z;
  z.vector.setZero();
  z.scalar = -(data.state.station2.param.xz - data.state.station1.param.xz);
  z.vector = z.scalar * ax.vector;
  z.a1() -= 1.0;
  z.a2() += 1.0;
  z.x1() += data.ax_result.ax;
  z.x2() -= data.ax_result.ax;

  data.xt.vector -= (data.xt_a2 / z.a2()) * z.vector;

  return true;
}

void BoundaryLayerTransitionSolver::setupLaminarTransitionSystem(TrdifData& data) {
  data.wf.scalar = (data.xt.scalar - data.state.station1.param.xz) /
                   (data.state.station2.param.xz - data.state.station1.param.xz);

  data.wf2 = (data.xt.scalar - data.state.station1.param.xz) /
             (data.state.station2.param.xz - data.state.station1.param.xz);
  data.wf_xt = 1.0 / (data.state.station2.param.xz - data.state.station1.param.xz);
  data.wf.vector = data.xt.vector * data.wf_xt;
  data.wf.x1() +=
      (data.wf2 - 1.0) / (data.state.station2.param.xz - data.state.station1.param.xz);
  data.wf.x2() -=
      data.wf2 / (data.state.station2.param.xz - data.state.station1.param.xz);

  data.wf1 = 1.0 - data.wf2;

  data.tt.scalar =
      data.state.station1.param.tz * data.wf1 + data.state.station2.param.tz * data.wf2;
  data.tt.vector =
      (data.state.station2.param.tz - data.state.station1.param.tz) * data.wf.vector;
  data.tt.t1() += data.wf1;
  data.tt.t2() += data.wf2;

  data.dt.scalar =
      data.state.station1.param.dz * data.wf1 + data.state.station2.param.dz * data.wf2;
  data.dt.vector =
      (data.state.station2.param.dz - data.state.station1.param.dz) * data.wf.vector;
  data.dt.d1() += data.wf1;
  data.dt.d2() += data.wf2;

  data.ut.scalar =
      data.state.station1.param.uz * data.wf1 + data.state.station2.param.uz * data.wf2;
  data.ut.vector =
      (data.state.station2.param.uz - data.state.station1.param.uz) * data.wf.vector;
  data.ut.u1() += data.wf1;
  data.ut.u2() += data.wf2;

  data.state.station2.param.xz = data.xt.scalar;
  data.state.station2.param.tz = data.tt.scalar;
  data.state.station2.param.dz = data.dt.scalar;
  data.state.station2.param.uz = data.ut.scalar;

  data.state.station2.param.amplz = data.blTransition.amcrit;
  data.state.station2.param.sz = 0.0;

  workflow_->blkin(data.state);
  data.state.station2 =
      boundaryLayerVariablesSolver.solve(data.state.station2, FlowRegimeEnum::Laminar);

  const SkinFrictionCoefficients laminarSkinFriction =
      workflow_->blmid(FlowRegimeEnum::Laminar);

  data.blc = blDiffSolver.solve(FlowRegimeEnum::Laminar, data.state,
                                laminarSkinFriction, data.blTransition.amcrit);

  data.blrez = data.blc.rhs;
  for (int k = 1; k < 3; ++k) {
    data.blm[k] = data.blc.d_msq[k] + data.blc.a2(k, 1) * data.tt.ms() +
                  data.blc.a2(k, 2) * data.dt.ms() +
                  data.blc.a2(k, 3) * data.ut.ms() +
                  data.blc.a2(k, 4) * data.xt.ms();
    data.blr[k] = data.blc.d_re[k] + data.blc.a2(k, 1) * data.tt.re() +
                  data.blc.a2(k, 2) * data.dt.re() +
                  data.blc.a2(k, 3) * data.ut.re() +
                  data.blc.a2(k, 4) * data.xt.re();
    data.blx[k] = data.blc.d_xi[k] + data.blc.a2(k, 1) * data.tt.xf() +
                  data.blc.a2(k, 2) * data.dt.xf() +
                  data.blc.a2(k, 3) * data.ut.xf() +
                  data.blc.a2(k, 4) * data.xt.xf();
  }

  data.bl1_transform = Eigen::Matrix<double, 4, 5>{
      {data.tt.a1(), data.tt.t1(), data.tt.d1(), data.tt.u1(), data.tt.x1()},
      {data.dt.a1(), data.dt.t1(), data.dt.d1(), data.dt.u1(), data.dt.x1()},
      {data.ut.a1(), data.ut.t1(), data.ut.d1(), data.ut.u1(), data.ut.x1()},
      {data.xt.a1(), data.xt.t1(), data.xt.d1(), data.xt.u1(), data.xt.x1()}};

  data.bl1.block<2, 5>(1, 0) =
      data.blc.a1.middleRows<2>(1) + data.blc.a2.block<2, 4>(1, 1) * data.bl1_transform;

  data.bl2_transform = Eigen::Matrix<double, 4, 4>{
      {data.tt.t2(), data.tt.d2(), data.tt.u2(), data.tt.x2()},
      {data.dt.t2(), data.dt.d2(), data.dt.u2(), data.dt.x2()},
      {data.ut.t2(), data.ut.d2(), data.ut.u2(), data.ut.x2()},
      {data.xt.t2(), data.xt.d2(), data.xt.u2(), data.xt.x2()}};

  data.bl2.block<2, 1>(1, 0).setZero();
  data.bl2.block<2, 4>(1, 1) =
      data.blc.a2.block<2, 4>(1, 1) * data.bl2_transform;
}

void BoundaryLayerTransitionSolver::setupTurbulentTransitionSystem(TrdifData& data) {
  data.state.station2 =
      boundaryLayerVariablesSolver.solve(data.state.station2, FlowRegimeEnum::Turbulent);

  data.ctr = 1.8 * std::exp(-3.3 / (data.state.station2.hkz.scalar - 1.0));
  data.ctr_hk2 = data.ctr * 3.3 / (data.state.station2.hkz.scalar - 1.0) /
                 (data.state.station2.hkz.scalar - 1.0);

  data.st = data.ctr * data.state.station2.cqz.scalar;
  data.st_tt = data.ctr * data.state.station2.cqz.t() +
               data.state.station2.cqz.scalar * data.ctr_hk2 *
                   data.state.station2.hkz.t();
  data.st_dt = data.ctr * data.state.station2.cqz.d() +
               data.state.station2.cqz.scalar * data.ctr_hk2 *
                   data.state.station2.hkz.d();
  data.st_ut = data.ctr * data.state.station2.cqz.u() +
               data.state.station2.cqz.scalar * data.ctr_hk2 *
                   data.state.station2.hkz.u();
  data.st_ms = data.ctr * data.state.station2.cqz.ms() +
               data.state.station2.cqz.scalar * data.ctr_hk2 *
                   data.state.station2.hkz.ms();
  data.st_re = data.ctr * data.state.station2.cqz.re();

  data.state.station2.param.amplz = 0.0;
  data.state.station2.param.sz = data.st;

  data.state.station2 =
      boundaryLayerVariablesSolver.solve(data.state.station2, FlowRegimeEnum::Turbulent);

  data.state.stepbl();
  data.state.station2 = boundaryLayerStore.restoreblData(2);

  const SkinFrictionCoefficients turbulentSkinFriction =
      workflow_->blmid(FlowRegimeEnum::Turbulent);

  data.blc = blDiffSolver.solve(FlowRegimeEnum::Turbulent, data.state,
                                turbulentSkinFriction, data.blTransition.amcrit);

  const Eigen::VectorXd common_st = Eigen::Vector3d{data.st_tt, data.st_dt, data.st_ut};
  const Eigen::VectorXd st1 =
      Eigen::Matrix<double, 5, 3>{{data.tt.a1(), data.dt.a1(), data.ut.a1()},
                                  {data.tt.t1(), data.dt.t1(), data.ut.t1()},
                                  {data.tt.d1(), data.dt.d1(), data.ut.d1()},
                                  {data.tt.u1(), data.dt.u1(), data.ut.u1()},
                                  {data.tt.x1(), data.dt.x1(), data.ut.x1()}} *
      common_st;
  const Eigen::VectorXd st2 =
      Eigen::Matrix<double, 5, 3>{{0, 0, 0},
                                  {data.tt.t2(), data.dt.t2(), data.ut.t2()},
                                  {data.tt.d2(), data.dt.d2(), data.ut.d2()},
                                  {data.tt.u2(), data.dt.u2(), data.ut.u2()},
                                  {data.tt.x2(), data.dt.x2(), data.ut.x2()}} *
      common_st;

  data.st_ms = data.st_tt * data.tt.ms() + data.st_dt * data.dt.ms() +
               data.st_ut * data.ut.ms() + data.st_ms;
  data.st_re = data.st_tt * data.tt.re() + data.st_dt * data.dt.re() +
               data.st_ut * data.ut.re() + data.st_re;
  data.st_xf = data.st_tt * data.tt.xf() + data.st_dt * data.dt.xf() +
               data.st_ut * data.ut.xf();

  Matrix<double, 5, 5> bt1_right = Matrix<double, 5, 5>::Zero();
  bt1_right.block<4, 5>(1, 0) = data.bl1_transform;
  bt1_right.row(0) = st1.transpose();

  Matrix<double, 5, 5> bt2_right = Matrix<double, 5, 5>::Zero();
  bt2_right.block<4, 4>(1, 1) = data.bl2_transform;
  bt2_right.row(0) = st2.transpose();

  data.bt1.block(0, 0, 3, 5) = data.blc.a1.block(0, 0, 3, 5) * bt1_right;
  data.bt2.block(0, 0, 3, 5) = data.blc.a1.block(0, 0, 3, 5) * bt2_right;
  data.bt2 += data.blc.a2;

  for (int k = 0; k < 3; ++k) {
    data.btrez[k] = data.blc.rhs[k];
    data.btm[k] = data.blc.d_msq[k] + data.blc.a1(k, 0) * data.st_ms +
                  data.blc.a1(k, 1) * data.tt.ms() +
                  data.blc.a1(k, 2) * data.dt.ms() +
                  data.blc.a1(k, 3) * data.ut.ms() +
                  data.blc.a1(k, 4) * data.xt.ms();
    data.btr[k] = data.blc.d_re[k] + data.blc.a1(k, 0) * data.st_re +
                  data.blc.a1(k, 1) * data.tt.re() +
                  data.blc.a1(k, 2) * data.dt.re() +
                  data.blc.a1(k, 3) * data.ut.re() +
                  data.blc.a1(k, 4) * data.xt.re();
    data.btx[k] = data.blc.d_xi[k] + data.blc.a1(k, 0) * data.st_xf +
                  data.blc.a1(k, 1) * data.tt.xf() +
                  data.blc.a1(k, 2) * data.dt.xf() +
                  data.blc.a1(k, 3) * data.ut.xf() +
                  data.blc.a1(k, 4) * data.xt.xf();
  }
}

void BoundaryLayerTransitionSolver::mergeTransitionSystems(TrdifData& data) {
  data.blc.rhs[0] = data.btrez[0];
  data.blc.rhs[1] = data.blrez[1] + data.btrez[1];
  data.blc.rhs[2] = data.blrez[2] + data.btrez[2];
  data.blc.d_msq[0] = data.btm[0];
  data.blc.d_msq[1] = data.blm[1] + data.btm[1];
  data.blc.d_msq[2] = data.blm[2] + data.btm[2];
  data.blc.d_re[0] = data.btr[0];
  data.blc.d_re[1] = data.blr[1] + data.btr[1];
  data.blc.d_re[2] = data.blr[2] + data.btr[2];
  data.blc.d_xi[0] = data.btx[0];
  data.blc.d_xi[1] = data.blx[1] + data.btx[1];
  data.blc.d_xi[2] = data.blx[2] + data.btx[2];
  data.blc.a1.row(0) = data.bt1.row(0);
  data.blc.a2.row(0) = data.bt2.row(0);
  data.blc.a1.middleRows(1, 2) =
      data.bl1.middleRows(1, 2) + data.bt1.middleRows(1, 2);
  data.blc.a2.middleRows(1, 2) =
      data.bl2.middleRows(1, 2) + data.bt2.middleRows(1, 2);

  data.state.station1 = boundaryLayerStore.restoreblData(1);
}

bool BoundaryLayerTransitionSolver::trdif() {
  //-----------------------------------------------
  //     sets up the newton system governing the
  //     transition interval.  equations governing
  //     the  laminar  part  x1 < xi < xt  and
  //     the turbulent part  xt < xi < x2
  //     are simply summed.
  //-----------------------------------------------
  TrdifData data(*workflow_);

  boundaryLayerStore.saveblData(data.state.station1, 1);
  boundaryLayerStore.saveblData(data.state.station2, 2);

  setupLaminarTransitionSystem(data);
  setupTurbulentTransitionSystem(data);
  mergeTransitionSystems(data);

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
  TrchekData data(*workflow_);

  boundaryLayerStore.saveblData(data.state.station2, 2);

  if (!iterateAmplification(data)) {
    Logger::instance().write("trchek2 - n2 convergence failed\n");
  }

  return resolveTransitionLocationAndSensitivities(data);
}
