#include "simulation/march/ue.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "XFoil.h"
#include "core/math_util.hpp"
#include "infrastructure/logger.hpp"

using EdgeVelocityFallbackMode =
    BoundaryLayerWorkflow::EdgeVelocityFallbackMode;
using MrchueStationContext = MarcherUe::MrchueStationContext;

namespace {
double adjustDisplacementForHkLimit(double displacementThickness,
                                    double momentumThickness, double msq,
                                    double hklim) {
  const double h = displacementThickness / momentumThickness;
  const auto hkin_result = boundary_layer::hkin(h, msq);
  const double dh = std::max(0.0, hklim - hkin_result.hk) / hkin_result.hk_h;
  return displacementThickness + dh * momentumThickness;
}

double computeMrchueDmax(const BoundaryLayerWorkflow &workflow, int side,
                         int stationIndex, const MrchueStationContext &ctx,
                         bool includeLaminarAmpTerm) {
  double dmax_local = std::max(std::fabs(workflow.blc.rhs[1] / ctx.thi),
                               std::fabs(workflow.blc.rhs[2] / ctx.dsi));

  if (includeLaminarAmpTerm &&
      stationIndex < workflow.lattice.get(side).profiles.transitionIndex) {
    dmax_local = std::max(dmax_local, std::fabs(workflow.blc.rhs[0] / 10.0));
  }
  if (stationIndex >= workflow.lattice.get(side).profiles.transitionIndex) {
    dmax_local = std::max(dmax_local, std::fabs(workflow.blc.rhs[0] / ctx.cti));
  }
  return dmax_local;
}

double computeMrchueRelaxation(double dmax_local) {
  if (dmax_local > 0.3) {
    return 0.3 / dmax_local;
  }
  return 1.0;
}

void applyMrchueStateDelta(const BoundaryLayerWorkflow &workflow, int side,
                           int stationIndex, MrchueStationContext &ctx,
                           double relaxation) {
  if (stationIndex >= workflow.lattice.get(side).profiles.transitionIndex) {
    ctx.cti += relaxation * workflow.blc.rhs[0];
  }
  ctx.thi += relaxation * workflow.blc.rhs[1];
  ctx.dsi += relaxation * workflow.blc.rhs[2];
  ctx.uei += relaxation * workflow.blc.rhs[3];
}

bool maybeSwitchToInverseMode(BoundaryLayerWorkflow &workflow, int side,
                              int stationIndex, MrchueStationContext &ctx,
                              double relaxation, double kHlmax, double kHtmax,
                              double &htarg, bool &direct,
                              std::stringstream &ss) {
  if (stationIndex == workflow.lattice.get(side).trailingEdgeIndex + 1) {
    return false;
  }

  const double msq =
      ctx.uei * ctx.uei * workflow.blCompressibility.hstinv /
      (workflow.blCompressibility.gm1bl *
       (1.0 - 0.5 * ctx.uei * ctx.uei * workflow.blCompressibility.hstinv));
  const double htest = (ctx.dsi + relaxation * workflow.blc.rhs[2]) /
                       (ctx.thi + relaxation * workflow.blc.rhs[1]);
  const auto hkin_result = boundary_layer::hkin(htest, msq);
  const double hktest = hkin_result.hk;
  const double hmax =
      (stationIndex < workflow.lattice.get(side).profiles.transitionIndex)
          ? kHlmax
          : kHtmax;

  direct = (hktest < hmax);
  if (direct) {
    return false;
  }

  htarg = workflow.calcHtarg(stationIndex, side, ctx.wake);
  if (ctx.wake) {
    htarg = std::max(htarg, 1.01);
  } else {
    htarg = std::max(htarg, hmax);
  }

  ss << "     mrchue: inverse mode at " << stationIndex
     << "    hk=" << std::fixed << std::setprecision(3) << htarg << "\n";
  Logger::instance().write(ss.str());
  ss.str("");
  return true;
}

void solveMrchueDirectSystem(BoundaryLayerWorkflow &workflow) {
  workflow.blc.a2(3, 0) = 0.0;
  workflow.blc.a2(3, 1) = 0.0;
  workflow.blc.a2(3, 2) = 0.0;
  workflow.blc.a2(3, 3) = 1.0;
  workflow.blc.rhs[3] = 0.0;
  workflow.blc.rhs =
      workflow.blc.a2.block(0, 0, 4, 4).fullPivLu().solve(workflow.blc.rhs);
}

void solveMrchueInverseSystem(BoundaryLayerWorkflow &workflow,
                              const MrchueStationContext &ctx) {
  workflow.blc.a2(3, 0) = 0.0;
  workflow.blc.a2(3, 1) = workflow.state.station2.hkz.t();
  workflow.blc.a2(3, 2) = workflow.state.station2.hkz.d();
  workflow.blc.a2(3, 3) = workflow.state.station2.hkz.u();
  workflow.blc.rhs[3] = ctx.htarg - workflow.state.station2.hkz.scalar;
  workflow.blc.rhs =
      workflow.blc.a2.block(0, 0, 4, 4).fullPivLu().solve(workflow.blc.rhs);
}

void clampAndAdjustMrchueStation(BoundaryLayerWorkflow &workflow, int side,
                                 int stationIndex, MrchueStationContext &ctx) {
  if (stationIndex >= workflow.lattice.get(side).profiles.transitionIndex) {
    ctx.cti = std::clamp(ctx.cti, 0.0000001, 0.30);
  }

  const double hklim =
      (stationIndex <= workflow.lattice.get(side).trailingEdgeIndex) ? 1.02
                                                                     : 1.00005;
  const double msq =
      ctx.uei * ctx.uei * workflow.blCompressibility.hstinv /
      (workflow.blCompressibility.gm1bl *
       (1.0 - 0.5 * ctx.uei * ctx.uei * workflow.blCompressibility.hstinv));
  double dsw = ctx.dsi - ctx.dswaki;
  dsw = adjustDisplacementForHkLimit(dsw, ctx.thi, msq, hklim);
  ctx.dsi = dsw + ctx.dswaki;
}
} // namespace

bool MarcherUe::mrchue(BoundaryLayerWorkflow &workflow, const Foil &foil,
                       const StagnationResult &stagnation) {
  return mrchue(workflow, workflow.state, foil, stagnation);
}

bool MarcherUe::mrchue(BoundaryLayerWorkflow &workflow,
                       BoundaryLayerState &state, const Foil &foil,
                       const StagnationResult &stagnation) {
  std::stringstream ss;
  for (int side = 1; side <= 2; ++side) {
    if (!marchMrchueSide(workflow, state, side, foil, stagnation, ss)) {
      return false;
    }
  }
  return true;
}

bool MarcherUe::marchMrchueSide(BoundaryLayerWorkflow &workflow,
                                BoundaryLayerState &state, int side,
                                const Foil &foil,
                                const StagnationResult &stagnation,
                                std::stringstream &ss) {
  ss << "    Side " << side << " ...\n";
  Logger::instance().write(ss.str());
  ss.str("");

  resetSideState(workflow, side, foil, stagnation);

  double thi = 0.0;
  double dsi = 0.0;
  double ami = 0.0;
  double cti = 0.0;
  initializeMrchueSide(workflow, side, thi, dsi, ami, cti);

  for (int stationIndex = 0;
       stationIndex < workflow.lattice.get(side).stationCount - 1;
       ++stationIndex) {
    MrchueStationContext ctx;
    prepareMrchueStationContext(workflow, side, stationIndex, ctx, thi, dsi,
                                ami, cti);
    bool converged = performMrchueNewtonLoop(workflow, side, stationIndex, ctx,
                                             foil.edge, ss);
    if (!converged) {
      handleMrchueStationFailure(workflow, side, stationIndex, ctx, ss);
    }
    storeMrchueStationState(workflow, side, stationIndex, ctx);

    ami = ctx.ami;
    cti = ctx.cti;
    thi = ctx.thi;
    dsi = ctx.dsi;

    if (stationIndex == workflow.lattice.get(side).trailingEdgeIndex) {
      thi = workflow.lattice.get(1)
                .profiles
                .momentumThickness[workflow.lattice.top.trailingEdgeIndex] +
            workflow.lattice.get(2)
                .profiles
                .momentumThickness[workflow.lattice.bottom.trailingEdgeIndex];
      dsi = workflow.lattice.get(1)
                .profiles
                .displacementThickness[workflow.lattice.top.trailingEdgeIndex] +
            workflow.lattice.get(2).profiles.displacementThickness
                [workflow.lattice.bottom.trailingEdgeIndex] +
            foil.edge.ante;
    }
  }

  return true;
}

void MarcherUe::initializeMrchueSide(BoundaryLayerWorkflow &workflow, int side,
                                     double &thi, double &dsi, double &ami,
                                     double &cti) {
  const double xsi = workflow.lattice.get(side).arcLengthCoordinates[0];
  const double uei = workflow.lattice.get(side).profiles.edgeVelocity[0];
  const double ucon = uei / xsi;
  const double tsq = 0.45 / (ucon * 6.0 * workflow.blReynolds.reybl);
  thi = std::sqrt(tsq);
  dsi = 2.2 * thi;
  ami = 0.0;
  cti = 0.03;
}

void MarcherUe::prepareMrchueStationContext(BoundaryLayerWorkflow &workflow,
                                            int side, int stationIndex,
                                            MrchueStationContext &ctx,
                                            double thi, double dsi, double ami,
                                            double cti) {
  ctx.simi = (stationIndex == 0);
  ctx.wake = stationIndex > workflow.lattice.get(side).trailingEdgeIndex;
  ctx.xsi = workflow.lattice.get(side).arcLengthCoordinates[stationIndex];
  ctx.uei = workflow.lattice.get(side).profiles.edgeVelocity[stationIndex];
  ctx.thi = thi;
  ctx.dsi = dsi;
  ctx.ami = ami;
  ctx.cti = cti;
  ctx.direct = true;
  ctx.dmax = 0.0;
  ctx.hmax = 0.0;
  ctx.htarg = 0.0;
  if (ctx.wake) {
    const int iw = stationIndex - workflow.lattice.get(side).trailingEdgeIndex;
    ctx.dswaki = workflow.wgap[iw - 1];
  } else {
    ctx.dswaki = 0.0;
  }
  workflow.flowRegime = determineRegimeForStation(workflow, side, stationIndex,
                                                  ctx.simi, ctx.wake);
}

bool MarcherUe::performMrchueNewtonLoop(BoundaryLayerWorkflow &workflow,
                                        int side, int stationIndex,
                                        MrchueStationContext &ctx,
                                        const Edge &edge,
                                        std::stringstream &ss) {
  constexpr double kHlmax = 3.8;
  constexpr double kHtmax = 2.5;

  bool converged = false;
  bool direct = true;
  double htarg = 0.0;
  double dmax_local = 0.0;
  workflow.flowRegime = determineRegimeForStation(workflow, side, stationIndex,
                                                  ctx.simi, ctx.wake);

  for (int itbl = 1; itbl <= 25; ++itbl) {
    workflow.state.current() =
        workflow.blprv(workflow.state.current(), ctx.xsi, ctx.ami, ctx.cti,
                       ctx.thi, ctx.dsi, ctx.dswaki, ctx.uei);

    workflow.blkin(workflow.state);

    if ((!ctx.simi) && (!(workflow.flowRegime == FlowRegimeEnum::Turbulent ||
                          workflow.flowRegime == FlowRegimeEnum::Wake))) {
      workflow.transitionSolver.trchek();
      ctx.ami = workflow.state.station2.param.amplz;

      if (workflow.flowRegime == FlowRegimeEnum::Transition) {
        workflow.lattice.get(side).profiles.transitionIndex = stationIndex;
        if (ctx.cti <= 0.0) {
          ctx.cti = 0.03;
          workflow.state.station2.param.sz = ctx.cti;
        }
      } else {
        workflow.lattice.get(side).profiles.transitionIndex = stationIndex + 2;
      }
    }

    if (stationIndex == workflow.lattice.get(side).trailingEdgeIndex + 1) {
      ctx.tte =
          workflow.lattice.get(1)
              .profiles
              .momentumThickness[workflow.lattice.top.trailingEdgeIndex] +
          workflow.lattice.get(2)
              .profiles
              .momentumThickness[workflow.lattice.bottom.trailingEdgeIndex];
      ctx.dte =
          workflow.lattice.get(1)
              .profiles
              .displacementThickness[workflow.lattice.top.trailingEdgeIndex] +
          workflow.lattice.get(2).profiles.displacementThickness
              [workflow.lattice.bottom.trailingEdgeIndex] +
          edge.ante;
      ctx.cte =
          (workflow.lattice.get(1)
                   .profiles
                   .skinFrictionCoeff[workflow.lattice.top.trailingEdgeIndex] *
               workflow.lattice.get(1)
                   .profiles
                   .momentumThickness[workflow.lattice.top.trailingEdgeIndex] +
           workflow.lattice.get(2).profiles.skinFrictionCoeff
                   [workflow.lattice.bottom.trailingEdgeIndex] *
               workflow.lattice.get(2).profiles.momentumThickness
                   [workflow.lattice.bottom.trailingEdgeIndex]) /
          ctx.tte;
      workflow.tesys(workflow.lattice.top.profiles,
                     workflow.lattice.bottom.profiles, edge);
    } else {
      workflow.blsys();
    }

    if (direct) {
      solveMrchueDirectSystem(workflow);
      dmax_local = computeMrchueDmax(workflow, side, stationIndex, ctx, true);
      const double rlx = computeMrchueRelaxation(dmax_local);

      if (maybeSwitchToInverseMode(workflow, side, stationIndex, ctx, rlx,
                                   kHlmax, kHtmax, htarg, direct, ss)) {
        continue;
      }

      applyMrchueStateDelta(workflow, side, stationIndex, ctx, rlx);
    } else {
      solveMrchueInverseSystem(workflow, ctx);
      dmax_local = computeMrchueDmax(workflow, side, stationIndex, ctx, false);
      const double rlx = computeMrchueRelaxation(dmax_local);
      applyMrchueStateDelta(workflow, side, stationIndex, ctx, rlx);
    }

    clampAndAdjustMrchueStation(workflow, side, stationIndex, ctx);

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

void MarcherUe::handleMrchueStationFailure(BoundaryLayerWorkflow &workflow,
                                           int side, int stationIndex,
                                           MrchueStationContext &ctx,
                                           std::stringstream &ss) {
  workflow.flowRegime = determineRegimeForStation(workflow, side, stationIndex,
                                                  ctx.simi, ctx.wake);

  ss << "     mrchue: convergence failed at " << stationIndex << ",  side "
     << side << ", res =" << std::fixed << std::setprecision(3) << ctx.dmax
     << "\n";
  Logger::instance().write(ss.str());
  ss.str("");

  resetStationKinematicsAfterFailure(
      workflow, side, stationIndex, ctx,
      EdgeVelocityFallbackMode::AverageNeighbors);

  {
    blData updatedCurrent =
        workflow.blprv(workflow.state.current(), ctx.xsi, ctx.ami, ctx.cti,
                       ctx.thi, ctx.dsi, ctx.dswaki, ctx.uei);
    workflow.state.current() = updatedCurrent;
  }
  workflow.blkin(workflow.state);

  workflow.checkTransitionIfNeeded(side, stationIndex, ctx.simi, 2, ctx.ami);
  syncStationRegimeStates(workflow, side, stationIndex, ctx.wake);
}

void MarcherUe::storeMrchueStationState(BoundaryLayerWorkflow &workflow,
                                        int side, int stationIndex,
                                        const MrchueStationContext &ctx) {
  storeStationStateCommon(workflow, side, stationIndex, ctx.ami, ctx.cti,
                          ctx.thi, ctx.dsi, ctx.uei, ctx.xsi, ctx.dswaki);
}
