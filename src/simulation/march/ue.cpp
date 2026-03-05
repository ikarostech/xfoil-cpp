#include "simulation/march/ue.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "XFoil.h"
#include "core/math_util.hpp"

using EdgeVelocityFallbackMode =
    BoundaryLayerWorkflow::EdgeVelocityFallbackMode;
using MrchueStationContext = MarcherUe::MrchueStationContext;

namespace
{
  double adjustDisplacementForHkLimit(double displacementThickness,
                                      double momentumThickness, double msq,
                                      double hklim)
  {
    const double h = displacementThickness / momentumThickness;
    const auto hkin_result = boundary_layer::hkin(h, msq);
    const double dh = std::max(0.0, hklim - hkin_result.hk) / hkin_result.hk_h;
    return displacementThickness + dh * momentumThickness;
  }

  double computeMrchueDmax(const BoundaryLayerWorkflow &workflow, int side,
                           int stationIndex, const MrchueStationContext &ctx,
                           bool includeLaminarAmpTerm)
  {
    const auto station = workflow.readStationModel(side, stationIndex);
    double dmax_local = std::max(std::fabs(workflow.readNewtonRhs(1) / ctx.thi),
                                 std::fabs(workflow.readNewtonRhs(2) / ctx.dsi));

    if (includeLaminarAmpTerm && stationIndex < station.transitionIndex)
    {
      dmax_local = std::max(dmax_local, std::fabs(workflow.readNewtonRhs(0) / 10.0));
    }
    if (stationIndex >= station.transitionIndex)
    {
      dmax_local =
          std::max(dmax_local, std::fabs(workflow.readNewtonRhs(0) / ctx.cti));
    }
    return dmax_local;
  }

  double computeMrchueRelaxation(double dmax_local)
  {
    if (dmax_local > 0.3)
    {
      return 0.3 / dmax_local;
    }
    return 1.0;
  }

  void applyMrchueStateDelta(const BoundaryLayerWorkflow &workflow, int side,
                             int stationIndex, MrchueStationContext &ctx,
                             double relaxation)
  {
    const auto station = workflow.readStationModel(side, stationIndex);
    if (stationIndex >= station.transitionIndex)
    {
      ctx.cti += relaxation * workflow.readNewtonRhs(0);
    }
    ctx.thi += relaxation * workflow.readNewtonRhs(1);
    ctx.dsi += relaxation * workflow.readNewtonRhs(2);
    ctx.uei += relaxation * workflow.readNewtonRhs(3);
  }

  bool maybeSwitchToInverseMode(BoundaryLayerWorkflow &workflow, int side,
                                int stationIndex, MrchueStationContext &ctx,
                                double relaxation, double kHlmax, double kHtmax,
                                double &htarg, bool &direct,
                                std::stringstream &ss)
  {
    const auto station = workflow.readStationModel(side, stationIndex);
    if (stationIndex == station.trailingEdgeIndex + 1)
    {
      return false;
    }

    const double msq =
        ctx.uei * ctx.uei * workflow.blCompressibility.hstinv /
        (workflow.blCompressibility.gm1bl *
         (1.0 - 0.5 * ctx.uei * ctx.uei * workflow.blCompressibility.hstinv));
    const double htest = (ctx.dsi + relaxation * workflow.readNewtonRhs(2)) /
                         (ctx.thi + relaxation * workflow.readNewtonRhs(1));
    const auto hkin_result = boundary_layer::hkin(htest, msq);
    const double hktest = hkin_result.hk;
    const double hmax = (stationIndex < station.transitionIndex) ? kHlmax : kHtmax;

    direct = (hktest < hmax);
    if (direct)
    {
      return false;
    }

    htarg = workflow.calcHtarg(stationIndex, side, ctx.isWake());
    if (ctx.isWake())
    {
      htarg = std::max(htarg, 1.01);
    }
    else
    {
      htarg = std::max(htarg, hmax);
    }

    ss << "     mrchue: inverse mode at " << stationIndex
       << "    hk=" << std::fixed << std::setprecision(3) << htarg << "\n";
    workflow.emitMarchInfoLog(ss.str());
    ss.str("");
    return true;
  }

  void clampAndAdjustMrchueStation(BoundaryLayerWorkflow &workflow, int side,
                                   int stationIndex, MrchueStationContext &ctx)
  {
    const auto station = workflow.readStationModel(side, stationIndex);
    if (stationIndex >= station.transitionIndex)
    {
      ctx.cti = std::clamp(ctx.cti, 0.0000001, 0.30);
    }

    const double hklim = (stationIndex <= station.trailingEdgeIndex) ? 1.02 : 1.00005;
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
                       const StagnationResult &stagnation)
{
  std::stringstream ss;
  for (int side = 1; side <= 2; ++side)
  {
    if (!marchMrchueSide(workflow, side, foil, stagnation, ss))
    {
      return false;
    }
  }
  return true;
}

bool MarcherUe::marchMrchueSide(BoundaryLayerWorkflow &workflow, int side,
                                const Foil &foil,
                                const StagnationResult &stagnation,
                                std::stringstream &ss)
{
  ss << "    Side " << side << " ...\n";
  workflow.emitMarchInfoLog(ss.str());
  ss.str("");

  workflow.resetSideState(side, foil, stagnation);

  initializeMrchueSide(workflow, side);

  const int stationCount = workflow.readSideStationCount(side);
  for (int stationIndex = 0;
       stationIndex < stationCount - 1;
       ++stationIndex)
  {
    MrchueStationContext ctx;
    prepareMrchueStationContext(workflow, side, stationIndex, ctx);
    bool converged = performMrchueNewtonLoop(workflow, side, stationIndex, ctx,
                                             foil.edge, ss);
    if (!converged)
    {
      handleMrchueStationFailure(workflow, side, stationIndex, ctx);
    }
    storeMrchueStationState(workflow, side, stationIndex, ctx);

    ami = ctx.ami;
    cti = ctx.cti;
    thi = ctx.thi;
    dsi = ctx.dsi;

    const auto station = workflow.readStationModel(side, stationIndex);
    if (stationIndex == station.trailingEdgeIndex)
    {
      const auto trailingEdge = workflow.readTrailingEdgeModel();
      thi = trailingEdge.topMomentumThickness + trailingEdge.bottomMomentumThickness;
      dsi = trailingEdge.topDisplacementThickness +
            trailingEdge.bottomDisplacementThickness + foil.edge.ante;
    }
  }

  return true;
}

void MarcherUe::initializeMrchueSide(BoundaryLayerWorkflow &workflow, int side)
{
  const auto station = workflow.readStationModel(side, 0);
  const double xsi = station.arcLength;
  const double uei = station.edgeVelocity;
  const double ucon = uei / xsi;
  const double tsq = 0.45 / (ucon * 6.0 * workflow.blReynolds.reybl);
  thi = std::sqrt(tsq);
  dsi = 2.2 * thi;
  ami = 0.0;
  cti = 0.03;
}

void MarcherUe::prepareMrchueStationContext(BoundaryLayerWorkflow &workflow,
                                            int side, int stationIndex,
                                            MrchueStationContext &ctx)
{
  const auto station = workflow.readStationModel(side, stationIndex);
  ctx.flowRegime = workflow.determineRegimeForStation(side, stationIndex);
  ctx.xsi = station.arcLength;
  ctx.uei = station.edgeVelocity;
  ctx.thi = thi;
  ctx.dsi = dsi;
  ctx.ami = ami;
  ctx.cti = cti;
  ctx.direct = true;
  ctx.dmax = 0.0;
  ctx.hmax = 0.0;
  ctx.htarg = 0.0;
  if (ctx.isWake())
  {
    ctx.dswaki = station.wakeGap;
  }
  else
  {
    ctx.dswaki = 0.0;
  }
  ctx.flowRegime = workflow.applyFlowRegimeCandidate(ctx.flowRegime);
}

bool MarcherUe::performMrchueNewtonLoop(BoundaryLayerWorkflow &workflow,
                                        int side, int stationIndex,
                                        MrchueStationContext &ctx,
                                        const Edge &edge,
                                        std::stringstream &ss)
{
  constexpr double kHlmax = 3.8;
  constexpr double kHtmax = 2.5;

  bool converged = false;
  bool direct = true;
  double htarg = 0.0;
  double dmax_local = 0.0;
  ctx.flowRegime = workflow.applyFlowRegimeCandidate(ctx.flowRegime);

  for (int itbl = 1; itbl <= 25; ++itbl)
  {
    workflow.state.current() =
        workflow.blprv(workflow.state.current(), ctx.xsi, ctx.ami, ctx.cti,
                       ctx.thi, ctx.dsi, ctx.dswaki, ctx.uei);

    workflow.blkin(workflow.state);
    ctx.flowRegime = workflow.currentFlowRegime();

    if ((!ctx.isSimilarity()) &&
        (!(ctx.flowRegime == FlowRegimeEnum::Turbulent ||
           ctx.flowRegime == FlowRegimeEnum::Wake)))
    {
      workflow.runTransitionCheckForMrchue(side, stationIndex, ctx.ami, ctx.cti);
      ctx.flowRegime = workflow.currentFlowRegime();
    }

    const auto station = workflow.readStationModel(side, stationIndex);
    if (stationIndex == station.trailingEdgeIndex + 1)
    {
      const auto trailingEdge = workflow.readTrailingEdgeModel();
      ctx.tte =
          trailingEdge.topMomentumThickness + trailingEdge.bottomMomentumThickness;
      ctx.dte = trailingEdge.topDisplacementThickness +
                trailingEdge.bottomDisplacementThickness + edge.ante;
      ctx.cte = (trailingEdge.topSkinFrictionCoeff *
                     trailingEdge.topMomentumThickness +
                 trailingEdge.bottomSkinFrictionCoeff *
                     trailingEdge.bottomMomentumThickness) /
                ctx.tte;
      workflow.solveTeSystemForCurrentProfiles(edge);
    }
    else
    {
      workflow.blsys();
    }

    if (direct)
    {
      workflow.solveMrchueDirectNewtonSystem();
      dmax_local = computeMrchueDmax(workflow, side, stationIndex, ctx, true);
      const double rlx = computeMrchueRelaxation(dmax_local);

      if (maybeSwitchToInverseMode(workflow, side, stationIndex, ctx, rlx,
                                   kHlmax, kHtmax, htarg, direct, ss))
      {
        continue;
      }

      applyMrchueStateDelta(workflow, side, stationIndex, ctx, rlx);
    }
    else
    {
      workflow.solveMrchueInverseNewtonSystem(ctx.htarg);
      dmax_local = computeMrchueDmax(workflow, side, stationIndex, ctx, false);
      const double rlx = computeMrchueRelaxation(dmax_local);
      applyMrchueStateDelta(workflow, side, stationIndex, ctx, rlx);
    }

    clampAndAdjustMrchueStation(workflow, side, stationIndex, ctx);

    if (dmax_local <= 0.00001)
    {
      converged = true;
      break;
    }
  }

  ctx.dmax = dmax_local;
  ctx.htarg = htarg;
  ctx.direct = direct;
  ctx.flowRegime = workflow.currentFlowRegime();
  return converged;
}

void MarcherUe::handleMrchueStationFailure(BoundaryLayerWorkflow &workflow,
                                           int side, int stationIndex,
                                           MrchueStationContext &ctx)
{
  workflow.emitMarchFailureLog("mrchue", side, stationIndex, ctx.dmax);

  workflow.recoverStationAfterFailure(side, stationIndex, ctx, ctx.ami,
                                      EdgeVelocityFallbackMode::AverageNeighbors,
                                      2);
}

void MarcherUe::storeMrchueStationState(BoundaryLayerWorkflow &workflow,
                                        int side, int stationIndex,
                                        const MrchueStationContext &ctx)
{
  workflow.storeStationStateCommon(side, stationIndex, ctx);
}
