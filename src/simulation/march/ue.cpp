#include "simulation/march/ue.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "XFoil.h"
#include "simulation/march/ue_numerics.hpp"

using EdgeVelocityFallbackMode = MrchueContext::EdgeVelocityFallbackMode;
using MrchueStationContext = MarcherUe::MrchueStationContext;
namespace ue_numerics = march::ue_numerics;

namespace
{
  double computeMrchueDmax(const MrchueContext &context, int side,
                           int stationIndex, const MrchueStationContext &ctx,
                           bool includeLaminarAmpTerm)
  {
    const auto station = context.readStationModel(side, stationIndex);
    return ue_numerics::computeDmax(
        context.readNewtonRhs(0), context.readNewtonRhs(1),
        context.readNewtonRhs(2), ctx.thi, ctx.dsi, ctx.cti,
        includeLaminarAmpTerm, stationIndex >= station.transitionIndex);
  }

  void applyMrchueStateDelta(const MrchueContext &context, int side,
                             int stationIndex, MrchueStationContext &ctx,
                             double relaxation)
  {
    const auto station = context.readStationModel(side, stationIndex);
    if (stationIndex >= station.transitionIndex)
    {
      ctx.cti += relaxation * context.readNewtonRhs(0);
    }
    ctx.thi += relaxation * context.readNewtonRhs(1);
    ctx.dsi += relaxation * context.readNewtonRhs(2);
    ctx.uei += relaxation * context.readNewtonRhs(3);
  }

  bool maybeSwitchToInverseMode(
      MrchueContext &context, int side, int stationIndex,
      MrchueStationContext &ctx, double relaxation, double kHlmax,
      double kHtmax, double &htarg, bool &direct,
      std::vector<Marcher::MarchEvent> &events)
  {
    const auto station = context.readStationModel(side, stationIndex);
    if (stationIndex == station.trailingEdgeIndex + 1)
    {
      return false;
    }

    const double msq = ue_numerics::computeMachSquared(
        ctx.uei, context.readBlCompressibilityHstinv(),
        context.readBlCompressibilityGm1bl());
    const double hktest = ue_numerics::computeProposedHk(
        ctx.dsi, ctx.thi, relaxation, context.readNewtonRhs(2),
        context.readNewtonRhs(1), msq);
    const double hmax = (stationIndex < station.transitionIndex) ? kHlmax : kHtmax;

    ctx.hmax = hmax;
    direct = (hktest < hmax);
    if (direct)
    {
      return false;
    }

    htarg = context.calcHtarg(stationIndex, side, ctx.isWake());
    if (ctx.isWake())
    {
      htarg = std::max(htarg, 1.01);
    }
    else
    {
      htarg = std::max(htarg, hmax);
    }

    std::ostringstream message;
    message << "     mrchue: inverse mode at " << stationIndex
            << "    hk=" << std::fixed << std::setprecision(3) << htarg
            << "\n";
    events.push_back({Marcher::MarchEvent::Kind::Info, message.str()});
    return true;
  }

  void clampAndAdjustMrchueStation(MrchueContext &context, int side,
                                   int stationIndex, MrchueStationContext &ctx)
  {
    const auto station = context.readStationModel(side, stationIndex);
    if (stationIndex >= station.transitionIndex)
    {
      ctx.cti = std::clamp(ctx.cti, 0.0000001, 0.30);
    }

    const double hklim = (stationIndex <= station.trailingEdgeIndex) ? 1.02 : 1.00005;
    const double msq = ue_numerics::computeMachSquared(
        ctx.uei, context.readBlCompressibilityHstinv(),
        context.readBlCompressibilityGm1bl());
    double dsw = ctx.dsi - ctx.dswaki;
    dsw = ue_numerics::adjustDisplacementForHkLimit(dsw, ctx.thi, msq, hklim);
    ctx.dsi = dsw + ctx.dswaki;
  }
} // namespace

bool MarcherUe::mrchue(MrchueContext &context, const Foil &foil,
                       const StagnationResult &stagnation)
{
  for (int side = 1; side <= 2; ++side)
  {
    if (!marchMrchueSide(context, side, foil, stagnation))
    {
      return false;
    }
  }
  return true;
}

bool MarcherUe::marchMrchueSide(MrchueContext &context, int side,
                                const Foil &foil,
                                const StagnationResult &stagnation)
{
  const auto publishEvent = [&](const MarchEvent &event) {
    if (event.kind == MarchEvent::Kind::Failure)
    {
      context.emitMarchInfoLog(event.message);
      return;
    }
    context.emitMarchInfoLog(event.message);
  };
  std::ostringstream message;
  message << "    Side " << side << " ...\n";
  publishEvent({MarchEvent::Kind::Info, message.str()});

  context.resetSideState(side, foil, stagnation);

  SideMarchState sideState = initializeMrchueSide(context, side);

  const int stationCount = context.readSideStationCount(side);
  for (int stationIndex = 0;
       stationIndex < stationCount - 1;
       ++stationIndex)
  {
    MrchueStationContext station;
    prepareMrchueStationContext(context, side, stationIndex, sideState, station);
    const MrchueStationContext resolvedStation = finalizeStationResult(
        performMrchueNewtonLoop(context, side, stationIndex, station, foil.edge),
        [&](StationMarchResult &result) {
          publishEvent(makeFailureEvent(result.recovery.phase, side, stationIndex,
                                        result.station.dmax));
          context.recoverStationAfterFailure(
              side, stationIndex, result.station, result.station.ami,
              result.recovery.edgeMode, result.recovery.laminarAdvance);
        },
        publishEvent,
        [&](const MrchueStationContext &resultStation) {
          storeMrchueStationState(context, side, stationIndex, resultStation);
        });

    updateSideStateFromStation(resolvedStation, sideState);
    updateSideStateForTrailingEdge(context, side, foil, stationIndex, sideState);
  }

  return true;
}

MarcherUe::SideMarchState MarcherUe::initializeMrchueSide(MrchueContext &context,
                                                          int side)
{
  const auto station = context.readStationModel(side, 0);
  const double xsi = station.arcLength;
  const double uei = station.edgeVelocity;
  const double ucon = uei / xsi;
  const double tsq = 0.45 / (ucon * 6.0 * context.readBlReynoldsReybl());
  SideMarchState sideState;
  sideState.thi = std::sqrt(tsq);
  sideState.dsi = 2.2 * sideState.thi;
  sideState.ami = 0.0;
  sideState.cti = 0.03;
  return sideState;
}

void MarcherUe::prepareMrchueStationContext(MrchueContext &context,
                                            int side, int stationIndex,
                                            const SideMarchState &sideState,
                                            MrchueStationContext &ctx)
{
  const auto station = context.readStationModel(side, stationIndex);
  ctx.flowRegime = context.determineRegimeForStation(side, stationIndex);
  ctx.xsi = station.arcLength;
  ctx.uei = station.edgeVelocity;
  ctx.thi = sideState.thi;
  ctx.dsi = sideState.dsi;
  ctx.ami = sideState.ami;
  ctx.cti = sideState.cti;
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
  ctx.flowRegime = context.applyFlowRegimeCandidate(ctx.flowRegime);
}

void MarcherUe::updateSideStateFromStation(const MrchueStationContext &ctx,
                                           SideMarchState &sideState)
{
  sideState.ami = ctx.ami;
  sideState.cti = ctx.cti;
  sideState.thi = ctx.thi;
  sideState.dsi = ctx.dsi;
}

void MarcherUe::updateSideStateForTrailingEdge(MrchueContext &context, int side,
                                               const Foil &foil,
                                               int stationIndex,
                                               SideMarchState &sideState)
{
  const auto station = context.readStationModel(side, stationIndex);
  if (stationIndex != station.trailingEdgeIndex)
  {
    return;
  }

  const auto trailingEdge = context.readTrailingEdgeModel();
  sideState.thi = trailingEdge.topMomentumThickness +
                  trailingEdge.bottomMomentumThickness;
  sideState.dsi = trailingEdge.topDisplacementThickness +
                  trailingEdge.bottomDisplacementThickness + foil.edge.ante;
}

MarcherUe::StationMarchResult
MarcherUe::performMrchueNewtonLoop(MrchueContext &context, int side,
                                   int stationIndex,
                                   MrchueStationContext station,
                                   const Edge &edge)
{
  constexpr double kHlmax = 3.8;
  constexpr double kHtmax = 2.5;

  double dmax_local = 0.0;
  std::vector<MarchEvent> events;
  station.flowRegime = context.applyFlowRegimeCandidate(station.flowRegime);
  auto &state = context.mutableState();

  for (int itbl = 1; itbl <= 25; ++itbl)
  {
    state.current() =
        context.blprv(state.current(), station.xsi, station.ami, station.cti,
                      station.thi, station.dsi, station.dswaki, station.uei);

    context.blkin(state);
    station.flowRegime = context.currentFlowRegime();

    if ((!station.isSimilarity()) &&
        (!(station.flowRegime == FlowRegimeEnum::Turbulent ||
           station.flowRegime == FlowRegimeEnum::Wake)))
    {
      context.runTransitionCheckForMrchue(side, stationIndex, station.ami,
                                          station.cti);
      station.flowRegime = context.currentFlowRegime();
    }

    const auto stationModel = context.readStationModel(side, stationIndex);
    if (stationIndex == stationModel.trailingEdgeIndex + 1)
    {
      const auto trailingEdge = context.readTrailingEdgeModel();
      station.tte =
          trailingEdge.topMomentumThickness + trailingEdge.bottomMomentumThickness;
      station.dte = trailingEdge.topDisplacementThickness +
                    trailingEdge.bottomDisplacementThickness + edge.ante;
      station.cte = (trailingEdge.topSkinFrictionCoeff *
                         trailingEdge.topMomentumThickness +
                     trailingEdge.bottomSkinFrictionCoeff *
                         trailingEdge.bottomMomentumThickness) /
                    station.tte;
      context.solveTeSystemForCurrentProfiles(edge);
    }
    else
    {
      context.blsys();
    }

    if (station.direct)
    {
      context.solveMrchueDirectNewtonSystem();
      dmax_local =
          computeMrchueDmax(context, side, stationIndex, station, true);
      const double rlx = ue_numerics::computeRelaxation(dmax_local);

      if (maybeSwitchToInverseMode(
              context, side, stationIndex, station, rlx, kHlmax, kHtmax,
              station.htarg, station.direct, events))
      {
        continue;
      }

      applyMrchueStateDelta(context, side, stationIndex, station, rlx);
    }
    else
    {
      const double effectiveHtarg = ue_numerics::computeInverseTarget(
          state.station2.hkz.scalar, station.htarg, station.hmax);
      context.solveMrchueInverseNewtonSystem(effectiveHtarg);
      dmax_local =
          computeMrchueDmax(context, side, stationIndex, station, false);
      const double rlx = ue_numerics::computeRelaxation(dmax_local);
      applyMrchueStateDelta(context, side, stationIndex, station, rlx);
    }

    clampAndAdjustMrchueStation(context, side, stationIndex, station);

    if (dmax_local <= 0.00001)
    {
      station.dmax = dmax_local;
      station.flowRegime = context.currentFlowRegime();
      return {station, true, {}, std::move(events)};
    }
  }

  station.dmax = dmax_local;
  station.flowRegime = context.currentFlowRegime();
  return {station,
          false,
          {true, "mrchue", EdgeVelocityFallbackMode::AverageNeighbors, 2},
          std::move(events)};
}

void MarcherUe::storeMrchueStationState(MrchueContext &context,
                                        int side, int stationIndex,
                                        const MrchueStationContext &ctx)
{
  context.storeStationStateCommon(side, stationIndex, ctx);
}
