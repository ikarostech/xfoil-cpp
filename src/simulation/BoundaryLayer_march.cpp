#include "BoundaryLayer.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "XFoil.h"

using BoundaryContext = BoundaryLayerWorkflow::MixedModeStationContext;

int BoundaryLayerWorkflow::resetSideState(int side, XFoil& xfoil) {
  const int previousTransition = lattice.transitionIndex.get(side);
  xfoil.xiforc = xfoil.xifset(side);
  xfoil.tran = false;
  xfoil.turb = false;
  lattice.transitionIndex.get(side) = lattice.trailingEdgeIndex.get(side);
  return previousTransition;
}

void BoundaryLayerWorkflow::storeStationStateCommon(
    int side, int stationIndex, double ami, double cti, double thi,
    double dsi, double uei, double xsi, double dswaki, XFoil& xfoil) {
  if (stationIndex < lattice.transitionIndex.get(side)) {
    lattice.ctau.get(side)[stationIndex] = ami;
  } else {
    lattice.ctau.get(side)[stationIndex] = cti;
  }
  lattice.thet.get(side)[stationIndex] = thi;
  lattice.dstr.get(side)[stationIndex] = dsi;
  lattice.uedg.get(side)[stationIndex] = uei;
  lattice.mass.get(side)[stationIndex] = dsi * uei;
  lattice.ctq.get(side)[stationIndex] = state.station2.cqz.scalar;

  {
    blData updatedCurrent =
        blprv(xfoil, state.current(), xsi, ami, cti, thi, dsi, dswaki, uei);
    state.current() = updatedCurrent;
  }
  xfoil.blkin(state);
  state.stepbl();

  if (xfoil.tran ||
      stationIndex == lattice.trailingEdgeIndex.get(side)) {
    xfoil.turb = true;
  }

  xfoil.tran = false;
}

double BoundaryLayerWorkflow::fallbackEdgeVelocity(
    int side, int stationIndex,
    EdgeVelocityFallbackMode edgeMode) const {
  switch (edgeMode) {
    case EdgeVelocityFallbackMode::UsePreviousStation:
      return lattice.uedg.get(side)[stationIndex - 1];
    case EdgeVelocityFallbackMode::AverageNeighbors: {
      double uei = lattice.uedg.get(side)[stationIndex];
      if (stationIndex < lattice.stationCount.get(side) - 1) {
        uei = 0.5 * (lattice.uedg.get(side)[stationIndex - 1] +
                     lattice.uedg.get(side)[stationIndex + 1]);
      }
      return uei;
    }
  }
  return lattice.uedg.get(side)[stationIndex];
}

void BoundaryLayerWorkflow::syncStationRegimeStates(int side,
                                                    int stationIndex,
                                                    bool wake,
                                                    XFoil& xfoil) {
  if (stationIndex < lattice.transitionIndex.get(side)) {
    state.station2 = blvar(state.station2, FlowRegimeEnum::Laminar);
    blmid(xfoil, FlowRegimeEnum::Laminar);
  }
  if (stationIndex >= lattice.transitionIndex.get(side)) {
    state.station2 = blvar(state.station2, FlowRegimeEnum::Turbulent);
    blmid(xfoil, FlowRegimeEnum::Turbulent);
  }
  if (wake) {
    state.station2 = blvar(state.station2, FlowRegimeEnum::Wake);
    blmid(xfoil, FlowRegimeEnum::Wake);
  }
}

bool BoundaryLayerWorkflow::mrchdu(XFoil& xfoil) {
  return mrchdu(state, lattice, xfoil);
}

bool BoundaryLayerWorkflow::mrchdu(BoundaryLayerState& state,
                                   [[maybe_unused]] BoundaryLayerLattice& lattice,
                                   XFoil& xfoil) {
  const double deps = 0.000005;
  const double senswt = 1000.0;

  double sens = 0.0;
  double sennew = 0.0;
  double ami = 0.0;

  for (int side = 1; side <= 2; ++side) {
    if (!marchBoundaryLayerSide(state, side, deps, senswt, sens, sennew, ami,
                                xfoil)) {
      return false;
    }
  }
  return true;
}

bool BoundaryLayerWorkflow::marchBoundaryLayerSide(
    BoundaryLayerState& state, int side, double deps, double senswt,
    double& sens, double& sennew, double& ami, XFoil& xfoil) {
  const int previousTransition = resetSideState(side, xfoil);

  for (int stationIndex = 0;
       stationIndex < lattice.stationCount.get(side) - 1; ++stationIndex) {
    if (!processBoundaryLayerStation(state, side, stationIndex,
                                     previousTransition, deps, senswt, sens,
                                     sennew, ami, xfoil)) {
      return false;
    }
  }

  return true;
}

bool BoundaryLayerWorkflow::processBoundaryLayerStation(
    BoundaryLayerState& state, int side, int stationIndex,
    int previousTransition, double deps, double senswt, double& sens,
    double& sennew, double& ami, XFoil& xfoil) {
  BoundaryContext ctx =
      xfoil.prepareMixedModeStation(side, stationIndex, previousTransition,
                                    ami);

  bool converged =
      xfoil.performMixedModeNewtonIteration(side, stationIndex,
                                            previousTransition, ctx, deps,
                                            senswt, sens, sennew, ami);
  if (!converged) {
    xfoil.handleMixedModeNonConvergence(side, stationIndex, ctx, ami);
  }

  sens = sennew;
  storeStationStateCommon(side, stationIndex, ctx.ami, ctx.cti, ctx.thi,
                          ctx.dsi, ctx.uei, ctx.xsi, ctx.dswaki, xfoil);

  if (XFoil::isCancelled()) {
    return false;
  }

  return true;
}

bool BoundaryLayerWorkflow::mrchue(XFoil& xfoil) {
  return mrchue(state, lattice, xfoil);
}

bool BoundaryLayerWorkflow::mrchue(BoundaryLayerState& state,
                                   [[maybe_unused]] BoundaryLayerLattice& lattice,
                                   XFoil& xfoil) {
  std::stringstream ss;
  for (int side = 1; side <= 2; ++side) {
    if (!marchMrchueSide(state, side, xfoil, ss)) {
      return false;
    }
  }
  return true;
}

bool BoundaryLayerWorkflow::marchMrchueSide(BoundaryLayerState& state,
                                            int side, XFoil& xfoil,
                                            std::stringstream& ss) {
  ss << "    Side " << side << " ...\n";
  xfoil.writeString(ss.str());
  ss.str("");

  resetSideState(side, xfoil);

  double thi = 0.0;
  double dsi = 0.0;
  double ami = 0.0;
  double cti = 0.0;
  initializeMrchueSide(side, thi, dsi, ami, cti, xfoil);

  for (int stationIndex = 0;
       stationIndex < lattice.stationCount.get(side) - 1; ++stationIndex) {
    MrchueStationContext ctx;
    prepareMrchueStationContext(side, stationIndex, ctx, thi, dsi, ami, cti,
                                xfoil);
    bool converged =
        performMrchueNewtonLoop(side, stationIndex, ctx, xfoil, ss);
    if (!converged) {
      handleMrchueStationFailure(side, stationIndex, ctx, xfoil, ss);
    }
    storeMrchueStationState(side, stationIndex, ctx, xfoil);

    ami = ctx.ami;
    cti = ctx.cti;
    thi = ctx.thi;
    dsi = ctx.dsi;

    if (stationIndex == lattice.trailingEdgeIndex.get(side)) {
      thi = lattice.thet.get(1)[lattice.trailingEdgeIndex.top] +
            lattice.thet.get(2)[lattice.trailingEdgeIndex.bottom];
      dsi = lattice.dstr.get(1)[lattice.trailingEdgeIndex.top] +
            lattice.dstr.get(2)[lattice.trailingEdgeIndex.bottom] +
            xfoil.foil.edge.ante;
    }

    if (XFoil::isCancelled()) {
      return false;
    }
  }

  return true;
}

void BoundaryLayerWorkflow::initializeMrchueSide(int side, double& thi,
                                                 double& dsi, double& ami,
                                                 double& cti, XFoil& xfoil) {
  const double xsi = lattice.xssi.get(side)[0];
  const double uei = lattice.uedg.get(side)[0];
  const double ucon = uei / xsi;
  const double tsq = 0.45 / (ucon * 6.0 * xfoil.reybl);
  thi = std::sqrt(tsq);
  dsi = 2.2 * thi;
  ami = 0.0;
  cti = 0.03;
}

void BoundaryLayerWorkflow::prepareMrchueStationContext(
    int side, int stationIndex, MrchueStationContext& ctx, double thi,
    double dsi, double ami, double cti, XFoil& xfoil) {
  ctx.simi = (stationIndex == 0);
  ctx.wake = stationIndex > lattice.trailingEdgeIndex.get(side);
  ctx.xsi = lattice.xssi.get(side)[stationIndex];
  ctx.uei = lattice.uedg.get(side)[stationIndex];
  ctx.thi = thi;
  ctx.dsi = dsi;
  ctx.ami = ami;
  ctx.cti = cti;
  ctx.direct = true;
  ctx.dmax = 0.0;
  ctx.hmax = 0.0;
  ctx.htarg = 0.0;
  if (ctx.wake) {
    const int iw = stationIndex - lattice.trailingEdgeIndex.get(side);
    ctx.dswaki = xfoil.wgap[iw - 1];
  } else {
    ctx.dswaki = 0.0;
  }
  xfoil.simi = ctx.simi;
  xfoil.wake = ctx.wake;
}

bool BoundaryLayerWorkflow::performMrchueNewtonLoop(
    int side, int stationIndex, MrchueStationContext& ctx, XFoil& xfoil,
    std::stringstream& ss) {
  constexpr double kHlmax = 3.8;
  constexpr double kHtmax = 2.5;

  bool converged = false;
  bool direct = true;
  double htarg = 0.0;
  double dmax_local = 0.0;
  xfoil.simi = ctx.simi;
  xfoil.wake = ctx.wake;

  for (int itbl = 1; itbl <= 25; ++itbl) {
    {
      blData updatedCurrent =
          blprv(xfoil, state.current(), ctx.xsi, ctx.ami, ctx.cti, ctx.thi,
                ctx.dsi, ctx.dswaki, ctx.uei);
      state.current() = updatedCurrent;
    }
    xfoil.blkin(state);

    if ((!ctx.simi) && (!xfoil.turb)) {
      xfoil.trchek();
      ctx.ami = state.station2.param.amplz;

      if (xfoil.tran) {
        lattice.transitionIndex.get(side) = stationIndex;
        if (ctx.cti <= 0.0) {
          ctx.cti = 0.03;
          state.station2.param.sz = ctx.cti;
        }
      } else {
        lattice.transitionIndex.get(side) = stationIndex + 2;
      }
    }

    if (stationIndex ==
        lattice.trailingEdgeIndex.get(side) + 1) {
      ctx.tte = lattice.thet.get(1)[lattice.trailingEdgeIndex.top] +
                lattice.thet.get(2)[lattice.trailingEdgeIndex.bottom];
      ctx.dte = lattice.dstr.get(1)[lattice.trailingEdgeIndex.top] +
                lattice.dstr.get(2)[lattice.trailingEdgeIndex.bottom] +
                xfoil.foil.edge.ante;
      ctx.cte =
          (lattice.ctau.get(1)[lattice.trailingEdgeIndex.top] *
               lattice.thet.get(1)[lattice.trailingEdgeIndex.top] +
           lattice.ctau.get(2)[lattice.trailingEdgeIndex.bottom] *
               lattice.thet.get(2)[lattice.trailingEdgeIndex.bottom]) /
          ctx.tte;
      tesys(xfoil, ctx.cte, ctx.tte, ctx.dte);
    } else {
      blsys(xfoil);
    }

    if (direct) {
      blc.a2(3, 0) = 0.0;
      blc.a2(3, 1) = 0.0;
      blc.a2(3, 2) = 0.0;
      blc.a2(3, 3) = 1.0;
      blc.rhs[3] = 0.0;
      blc.rhs =
          blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);

      dmax_local =
          std::max(std::fabs(blc.rhs[1] / ctx.thi),
                   std::fabs(blc.rhs[2] / ctx.dsi));
      if (stationIndex < lattice.transitionIndex.get(side)) {
        dmax_local =
            std::max(dmax_local, std::fabs(blc.rhs[0] / 10.0));
      }
      if (stationIndex >= lattice.transitionIndex.get(side)) {
        dmax_local =
            std::max(dmax_local, std::fabs(blc.rhs[0] / ctx.cti));
      }

      double rlx = 1.0;
      if (dmax_local > 0.3) {
        rlx = 0.3 / dmax_local;
      }

      if (stationIndex != lattice.trailingEdgeIndex.get(side) + 1) {
        const double msq =
            ctx.uei * ctx.uei * xfoil.hstinv /
            (xfoil.gm1bl *
             (1.0 - 0.5 * ctx.uei * ctx.uei * xfoil.hstinv));
        const double htest =
            (ctx.dsi + rlx * blc.rhs[2]) /
            (ctx.thi + rlx * blc.rhs[1]);
        const auto hkin_result = boundary_layer::hkin(htest, msq);
        const double hktest = hkin_result.hk;

        const double hmax =
            (stationIndex < lattice.transitionIndex.get(side))
                ? kHlmax
                : kHtmax;
        direct = (hktest < hmax);
        if (!direct) {
          htarg = xfoil.calcHtarg(stationIndex, side, ctx.wake);
          if (ctx.wake) {
            htarg = std::max(htarg, 1.01);
          } else {
            htarg = std::max(htarg, hmax);
          }

          ss << "     mrchue: inverse mode at " << stationIndex
             << "    hk=" << std::fixed << std::setprecision(3) << htarg
             << "\n";
          xfoil.writeString(ss.str());
          ss.str("");
          continue;
        }
      }

      if (stationIndex >= lattice.transitionIndex.get(side)) {
        ctx.cti += rlx * blc.rhs[0];
      }
      ctx.thi += rlx * blc.rhs[1];
      ctx.dsi += rlx * blc.rhs[2];
      ctx.uei += rlx * blc.rhs[3];
    } else {
      blc.a2(3, 0) = 0.0;
      blc.a2(3, 1) = state.station2.hkz.t();
      blc.a2(3, 2) = state.station2.hkz.d();
      blc.a2(3, 3) = state.station2.hkz.u();
      blc.rhs[3] = ctx.htarg - state.station2.hkz.scalar;
      blc.rhs =
          blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);

      dmax_local =
          std::max(std::fabs(blc.rhs[1] / ctx.thi),
                   std::fabs(blc.rhs[2] / ctx.dsi));
      if (stationIndex >= lattice.transitionIndex.get(side)) {
        dmax_local =
            std::max(dmax_local, std::fabs(blc.rhs[0] / ctx.cti));
      }

      double rlx = 1.0;
      if (dmax_local > 0.3) {
        rlx = 0.3 / dmax_local;
      }

      if (stationIndex >= lattice.transitionIndex.get(side)) {
        ctx.cti += rlx * blc.rhs[0];
      }
      ctx.thi += rlx * blc.rhs[1];
      ctx.dsi += rlx * blc.rhs[2];
      ctx.uei += rlx * blc.rhs[3];
    }

    if (stationIndex >= lattice.transitionIndex.get(side)) {
      ctx.cti = std::clamp(ctx.cti, 0.0000001, 0.30);
    }

    const double hklim =
        (stationIndex <= lattice.trailingEdgeIndex.get(side)) ? 1.02
                                                              : 1.00005;
    const double msq =
        ctx.uei * ctx.uei * xfoil.hstinv /
        (xfoil.gm1bl *
         (1.0 - 0.5 * ctx.uei * ctx.uei * xfoil.hstinv));
    double dsw = ctx.dsi - ctx.dswaki;
    xfoil.dslim(dsw, ctx.thi, msq, hklim);
    ctx.dsi = dsw + ctx.dswaki;

    if (dmax_local <= 0.00001) {
      converged = true;
      break;
    }
  }

  ctx.dmax = dmax_local;
  ctx.htarg = htarg;
  ctx.direct = direct;
  return converged;
}

void BoundaryLayerWorkflow::handleMrchueStationFailure(
    int side, int stationIndex, MrchueStationContext& ctx, XFoil& xfoil,
    std::stringstream& ss) {
  xfoil.simi = ctx.simi;
  xfoil.wake = ctx.wake;

  ss << "     mrchue: convergence failed at " << stationIndex << ",  side "
     << side << ", res =" << std::fixed << std::setprecision(3) << ctx.dmax
     << "\n";
  xfoil.writeString(ss.str());
  ss.str("");

  resetStationKinematicsAfterFailure(
      side, stationIndex, ctx,
      EdgeVelocityFallbackMode::AverageNeighbors);

  {
    blData updatedCurrent =
        blprv(xfoil, state.current(), ctx.xsi, ctx.ami, ctx.cti, ctx.thi,
              ctx.dsi, ctx.dswaki, ctx.uei);
    state.current() = updatedCurrent;
  }
  xfoil.blkin(state);

  xfoil.checkTransitionIfNeeded(side, stationIndex, ctx.simi, 2, ctx.ami);
  syncStationRegimeStates(side, stationIndex, ctx.wake, xfoil);
}

void BoundaryLayerWorkflow::storeMrchueStationState(
    int side, int stationIndex, const MrchueStationContext& ctx, XFoil& xfoil) {
  storeStationStateCommon(side, stationIndex, ctx.ami, ctx.cti, ctx.thi,
                          ctx.dsi, ctx.uei, ctx.xsi, ctx.dswaki, xfoil);
}
