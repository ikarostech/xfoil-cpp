#include "solver/march/ue.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "solver/march/ue_numerics.hpp"

using EdgeVelocityFallbackMode = MrchueContext::EdgeVelocityFallbackMode;
using MrchueStationContext = MarcherUe::MrchueStationContext;
namespace ue_numerics = march::ue_numerics;

namespace
{
  double computeMrchueDmax(const MrchueContext &context,
                           const MarcherUe::StationInput &input,
                           const MrchueStationContext &ctx,
                           bool includeLaminarAmpTerm)
  {
    return ue_numerics::computeDmax(
        context.readNewtonRhs(0), context.readNewtonRhs(1),
        context.readNewtonRhs(2), ctx.thi, ctx.dsi, ctx.cti,
        includeLaminarAmpTerm,
        input.stationIndex >= input.stationModel.transitionIndex);
  }

  void applyMrchueStateDelta(const MrchueContext &context,
                             const MarcherUe::StationInput &input,
                             MrchueStationContext &ctx,
                             double relaxation)
  {
    if (input.stationIndex >= input.stationModel.transitionIndex)
    {
      ctx.cti += relaxation * context.readNewtonRhs(0);
    }
    ctx.thi += relaxation * context.readNewtonRhs(1);
    ctx.dsi += relaxation * context.readNewtonRhs(2);
    ctx.uei += relaxation * context.readNewtonRhs(3);
  }

  bool maybeSwitchToInverseMode(const MarcherUe::StationInput &input,
                                MrchueContext &context,
                                MrchueStationContext &ctx,
                                double relaxation, double kHlmax,
                                double kHtmax, double &htarg, bool &direct,
                                std::vector<Marcher::MarchEvent> &events)
  {
    if (input.stationIndex == input.stationModel.trailingEdgeIndex + 1)
    {
      return false;
    }

    const double msq = ue_numerics::computeMachSquared(
        ctx.uei, input.hstinv, input.gm1bl);
    const double hktest = ue_numerics::computeProposedHk(
        ctx.dsi, ctx.thi, relaxation, context.readNewtonRhs(2),
        context.readNewtonRhs(1), msq);
    const double hmax =
        (input.stationIndex < input.stationModel.transitionIndex) ? kHlmax
                                                                  : kHtmax;

    ctx.hmax = hmax;
    direct = (hktest < hmax);
    if (direct)
    {
      return false;
    }

    htarg = context.calcHtarg(input.stationIndex, input.side, ctx.isWake());
    if (ctx.isWake())
    {
      htarg = std::max(htarg, 1.01);
    }
    else
    {
      htarg = std::max(htarg, hmax);
    }

    std::ostringstream message;
    message << "     mrchue: inverse mode at " << input.stationIndex
            << "    hk=" << std::fixed << std::setprecision(3) << htarg
            << "\n";
    events.push_back({Marcher::MarchEvent::Kind::Info, message.str()});
    return true;
  }

  void clampAndAdjustMrchueStation(const MarcherUe::StationInput &input,
                                   MrchueStationContext &ctx)
  {
    if (input.stationIndex >= input.stationModel.transitionIndex)
    {
      ctx.cti = std::clamp(ctx.cti, 0.0000001, 0.30);
    }

    const double hklim =
        (input.stationIndex <= input.stationModel.trailingEdgeIndex) ? 1.02
                                                                     : 1.00005;
    const double msq =
        ue_numerics::computeMachSquared(ctx.uei, input.hstinv, input.gm1bl);
    double dsw = ctx.dsi - ctx.dswaki;
    dsw = ue_numerics::adjustDisplacementForHkLimit(dsw, ctx.thi, msq, hklim);
    ctx.dsi = dsw + ctx.dswaki;
  }
} // namespace

bool MarcherUe::mrchue(MrchueContext &context, const Foil &foil,
                       const StagnationFeature &stagnation)
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
                                const StagnationFeature &stagnation)
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

  const SideInput sideInput =
      makeSideInput(context, side, foil, stagnation);
  SideMarchState sideState = initializeMrchueSide(sideInput);

  for (int stationIndex = 0;
       stationIndex < sideInput.stationCount - 1;
       ++stationIndex)
  {
    const StationInput input =
        makeStationInput(context, side, stationIndex, foil);
    MrchueStationContext station;
    prepareMrchueStationContext(context, input, sideState, station);
    const MrchueStationContext resolvedStation = finalizeStationResult(
        performMrchueNewtonLoop(context, input, station, input.edge),
        [&](StationMarchResult &result) {
          publishEvent(makeFailureEvent(result.recovery.phase, input.side,
                                        input.stationIndex, result.station.dmax));
          context.recoverStationAfterFailure(
              input.side, input.stationIndex, result.station, result.station.ami,
              result.recovery.edgeMode, result.recovery.laminarAdvance);
        },
        publishEvent,
        [&](const MrchueStationContext &resultStation) {
          storeMrchueStationState(context, input.side, input.stationIndex,
                                 resultStation);
        });

    updateSideStateFromStation(resolvedStation, sideState);
    updateSideStateForTrailingEdge(input, foil, sideState);
  }

  return true;
}

MarcherUe::SideInput MarcherUe::makeSideInput(MrchueContext &context, int side,
                                              const Foil &foil,
                                              const StagnationFeature &stagnation) const
{
  context.resetSideState(side, foil, stagnation);
  return {side, context.readSideStationCount(side),
          context.readStationModel(side, 0), context.readBlReynoldsReybl()};
}

MarcherUe::SideMarchState
MarcherUe::initializeMrchueSide(const SideInput &input)
{
  const double xsi = input.leadingStationModel.arcLength;
  const double uei = input.leadingStationModel.edgeVelocity;
  const double ucon = uei / xsi;
  const double tsq = 0.45 / (ucon * 6.0 * input.reybl);
  SideMarchState sideState;
  sideState.thi = std::sqrt(tsq);
  sideState.dsi = 2.2 * sideState.thi;
  sideState.ami = 0.0;
  sideState.cti = 0.03;
  return sideState;
}

void MarcherUe::prepareMrchueStationContext(MrchueContext &context,
                                            const StationInput &input,
                                            const SideMarchState &sideState,
                                            MrchueStationContext &ctx)
{
  ctx.flowRegime =
      context.determineRegimeForStation(input.side, input.stationIndex);
  ctx.xsi = input.stationModel.arcLength;
  ctx.uei = input.stationModel.edgeVelocity;
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
    ctx.dswaki = input.stationModel.wakeGap;
  }
  else
  {
    ctx.dswaki = 0.0;
  }
  ctx.flowRegime = context.applyFlowRegimeCandidate(ctx.flowRegime);
}

MarcherUe::StationInput MarcherUe::makeStationInput(MrchueContext &context,
                                                    int side, int stationIndex,
                                                    const Foil &foil) const
{
  return {side,
          stationIndex,
          foil.edge,
          context.readStationModel(side, stationIndex),
          context.readTrailingEdgeModel(),
          context.readBlCompressibilityHstinv(),
          context.readBlCompressibilityGm1bl()};
}

void MarcherUe::updateSideStateFromStation(const MrchueStationContext &ctx,
                                           SideMarchState &sideState)
{
  sideState.ami = ctx.ami;
  sideState.cti = ctx.cti;
  sideState.thi = ctx.thi;
  sideState.dsi = ctx.dsi;
}

void MarcherUe::updateSideStateForTrailingEdge(const StationInput &input,
                                               const Foil &foil,
                                               SideMarchState &sideState)
{
  if (input.stationIndex != input.stationModel.trailingEdgeIndex)
  {
    return;
  }

  sideState.thi = input.trailingEdgeModel.topMomentumThickness +
                  input.trailingEdgeModel.bottomMomentumThickness;
  sideState.dsi = input.trailingEdgeModel.topDisplacementThickness +
                  input.trailingEdgeModel.bottomDisplacementThickness +
                  foil.edge.ante;
}

MarcherUe::StationMarchResult
MarcherUe::performMrchueNewtonLoop(MrchueContext &context,
                                   const StationInput &input,
                                   MrchueStationContext station,
                                   const Edge &edge)
{
  constexpr double kHlmax = 3.8;
  constexpr double kHtmax = 2.5;

  double dmax_local = 0.0;
  std::vector<MarchEvent> events;
  station.flowRegime = context.applyFlowRegimeCandidate(station.flowRegime);

  for (int itbl = 1; itbl <= 25; ++itbl)
  {
    context.refreshCurrentStationState(station.xsi, station.ami, station.cti,
                                       station.thi, station.dsi,
                                       station.dswaki, station.uei);
    station.flowRegime = context.currentFlowRegime();

    if ((!station.isSimilarity()) &&
        (!(station.flowRegime == FlowRegimeEnum::Turbulent ||
           station.flowRegime == FlowRegimeEnum::Wake)))
    {
      context.runTransitionCheckForMrchue(input.side, input.stationIndex,
                                          station.ami, station.cti);
      station.flowRegime = context.currentFlowRegime();
    }

    if (input.stationIndex == input.stationModel.trailingEdgeIndex + 1)
    {
      station.tte =
          input.trailingEdgeModel.topMomentumThickness +
          input.trailingEdgeModel.bottomMomentumThickness;
      station.dte = input.trailingEdgeModel.topDisplacementThickness +
                    input.trailingEdgeModel.bottomDisplacementThickness +
                    edge.ante;
      station.cte = (input.trailingEdgeModel.topSkinFrictionCoeff *
                         input.trailingEdgeModel.topMomentumThickness +
                     input.trailingEdgeModel.bottomSkinFrictionCoeff *
                         input.trailingEdgeModel.bottomMomentumThickness) /
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
      dmax_local = computeMrchueDmax(context, input, station, true);
      const double rlx = ue_numerics::computeRelaxation(dmax_local);

      if (maybeSwitchToInverseMode(
              input, context, station, rlx, kHlmax, kHtmax,
              station.htarg, station.direct, events))
      {
        continue;
      }

      applyMrchueStateDelta(context, input, station, rlx);
    }
    else
    {
      const double effectiveHtarg = ue_numerics::computeInverseTarget(
          context.readCurrentShapeFactor(), station.htarg, station.hmax);
      context.solveMrchueInverseNewtonSystem(effectiveHtarg);
      dmax_local = computeMrchueDmax(context, input, station, false);
      const double rlx = ue_numerics::computeRelaxation(dmax_local);
      applyMrchueStateDelta(context, input, station, rlx);
    }

    clampAndAdjustMrchueStation(input, station);

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
