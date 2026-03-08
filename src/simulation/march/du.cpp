#include "simulation/march/du.hpp"

#include <algorithm>
#include <cmath>

#include "XFoil.h"
#include "core/math_util.hpp"
#include "simulation/march/du_numerics.hpp"

using BoundaryContext = MrchduContext::MixedModeStationContext;
using EdgeVelocityFallbackMode = MrchduContext::EdgeVelocityFallbackMode;
namespace du_numerics = march::du_numerics;

MrchduContext::MixedModeStationContext
MarcherDu::prepareMixedModeStation(MrchduContext &context,
                                   SideMarchState &sideState,
                                   const StationInput &input)
{
  BoundaryContext ctx;
  const auto &station = input.stationModel;

  ctx.flowRegime =
      context.determineRegimeForStation(input.side, input.stationIndex);
  ctx.xsi = station.arcLength;
  ctx.uei = station.edgeVelocity;
  ctx.thi = station.momentumThickness;
  ctx.dsi = station.displacementThickness;

  if (input.stationIndex < input.previousTransition)
  {
    sideState.ami = station.skinFrictionCoeff;
    ctx.cti = 0.03;
  }
  else
  {
    ctx.cti = station.skinFrictionCoeff;
    if (ctx.cti <= 0.0)
    {
      ctx.cti = 0.03;
    }
  }
  ctx.ami = sideState.ami;

  if (ctx.isWake())
  {
    ctx.dswaki = station.wakeGap;
  }
  else
  {
    ctx.dswaki = 0.0;
  }

  ctx.dsi = du_numerics::computeInitialDisplacement(
      ctx.dsi, ctx.dswaki, ctx.thi, input.stationIndex,
      station.trailingEdgeIndex);

  ctx.flowRegime = context.applyFlowRegimeCandidate(ctx.flowRegime);

  return ctx;
}

bool MarcherDu::mrchdu(MrchduContext &context, const Foil &foil,
                       const StagnationResult &stagnation)
{
  return mrchdu(context, context.mutableState(), foil, stagnation);
}

bool MarcherDu::mrchdu(MrchduContext &context,
                       BoundaryLayerState &state, const Foil &foil,
                       const StagnationResult &stagnation)
{
  for (int side = 1; side <= 2; ++side)
  {
    SideMarchState sideState;
    if (!marchBoundaryLayerSide(context, state, side, sideState, foil,
                                stagnation))
    {
      return false;
    }
  }
  return true;
}

MarcherDu::SideInput MarcherDu::makeSideInput(MrchduContext &context, int side,
                                              const Foil &foil,
                                              const StagnationResult &stagnation) const
{
  return {side, context.readSideStationCount(side),
          context.resetSideState(side, foil, stagnation)};
}

bool MarcherDu::marchBoundaryLayerSide(MrchduContext &context,
                                       BoundaryLayerState &state, int side,
                                       SideMarchState &sideState,
                                       const Foil &foil,
                                       const StagnationResult &stagnation)
{
  const SideInput sideInput =
      makeSideInput(context, side, foil, stagnation);

  for (int stationIndex = 0;
       stationIndex < sideInput.stationCount - 1;
       ++stationIndex)
  {
    if (!processBoundaryLayerStation(context, state, sideState, side,
                                     stationIndex,
                                     sideInput.previousTransition, foil))
    {
      return false;
    }
  }

  return true;
}

bool MarcherDu::processBoundaryLayerStation(
    MrchduContext &context, BoundaryLayerState &state,
    SideMarchState &sideState, int side,
    int stationIndex, int previousTransition, const Foil &foil)
{
  const StationInput input = makeStationInput(
      context, side, stationIndex, previousTransition, foil);
  const auto publishEvent = [&](const MarchEvent &event) {
    context.emitMarchInfoLog(event.message);
  };
  BoundaryContext station = prepareMixedModeStation(context, sideState, input);
  const BoundaryContext resolvedStation = finalizeStationResult(
      performMixedModeNewtonIteration(context, sideState, input, station),
      [&](StationMarchResult &result) {
        publishEvent(makeFailureEvent(result.recovery.phase, input.side,
                                      input.stationIndex,
                                      result.station.dmax));
        context.recoverStationAfterFailure(
            input.side, input.stationIndex, result.station, sideState.ami,
            result.recovery.edgeMode, result.recovery.laminarAdvance);
      },
      publishEvent,
      [&](const BoundaryContext &resultStation) {
        context.storeStationStateCommon(input.side, input.stationIndex,
                                        resultStation);
      });

  sideState.sens = sideState.sennew;
  (void)resolvedStation;
  return true;
}

MarcherDu::StationInput MarcherDu::makeStationInput(
    MrchduContext &context, int side, int stationIndex, int previousTransition,
    const Foil &foil) const
{
  return {side,       stationIndex, previousTransition, foil.edge,
          context.readStationModel(side, stationIndex),
          context.isStartOfWake(side, stationIndex)};
}

MarcherDu::StationMarchResult MarcherDu::performMixedModeNewtonIteration(
    MrchduContext &context, SideMarchState &sideState, const StationInput &input,
    MrchduContext::MixedModeStationContext station)
{
  double ueref = 0.0;
  double hkref = 0.0;
  station.flowRegime = context.applyFlowRegimeCandidate(station.flowRegime);
  auto &state = context.mutableState();

  for (int itbl = 1; itbl <= 25; ++itbl)
  {
    blData updatedCurrent =
        context.blprv(state.current(), station.xsi, sideState.ami,
                      station.cti, station.thi, station.dsi, station.dswaki,
                      station.uei);
    state.current() = updatedCurrent;
    context.blkin(state);
    station.flowRegime = context.currentFlowRegime();

    context.checkTransitionIfNeeded(input.side, input.stationIndex,
                                    station.isSimilarity(), 1,
                                    sideState.ami);
    station.flowRegime = context.currentFlowRegime();

    context.updateSystemMatricesForStation(input.edge, input.side,
                                           input.stationIndex, station);

    if (itbl == 1)
    {
      context.initializeFirstIterationState(input.side, input.stationIndex,
                                            input.previousTransition, station,
                                            ueref, hkref);
    }

    if (station.isSimilarity() || input.startOfWake)
    {
      context.configureSimilarityRow(ueref);
    }
    else
    {
      const auto updateMode = du_numerics::computeSensitivityUpdateMode(itbl);
      context.configureViscousRow(hkref, ueref, senswt, updateMode.reset,
                                  updateMode.average, sideState.sens,
                                  sideState.sennew);
    }

    if (context.applyMixedModeNewtonStep(input.side, input.stationIndex,
                                         sideState.ami, station))
    {
      return {station, true, {}, {}};
    }
  }

  return {station,
          false,
          {true, "mrchdu",
           MrchduContext::EdgeVelocityFallbackMode::UsePreviousStation, 2},
          {}};
}
