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
                                   SideMarchState &sideState, int side,
                                   int stationIndex, int previousTransition)
{
  BoundaryContext ctx;
  const auto station = context.readStationModel(side, stationIndex);

  ctx.flowRegime = context.determineRegimeForStation(side, stationIndex);
  ctx.xsi = station.arcLength;
  ctx.uei = station.edgeVelocity;
  ctx.thi = station.momentumThickness;
  ctx.dsi = station.displacementThickness;

  if (stationIndex < previousTransition)
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
      ctx.dsi, ctx.dswaki, ctx.thi, stationIndex, station.trailingEdgeIndex);

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

bool MarcherDu::marchBoundaryLayerSide(MrchduContext &context,
                                       BoundaryLayerState &state, int side,
                                       SideMarchState &sideState,
                                       const Foil &foil,
                                       const StagnationResult &stagnation)
{
  const int previousTransition = context.resetSideState(side, foil, stagnation);
  const int stationCount = context.readSideStationCount(side);

  for (int stationIndex = 0;
       stationIndex < stationCount - 1;
       ++stationIndex)
  {
    if (!processBoundaryLayerStation(context, state, sideState, side,
                                     stationIndex,
                                     previousTransition, foil))
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
  BoundaryContext ctx = prepareMixedModeStation(context, sideState, side,
                                                stationIndex,
                                                previousTransition);
  bool converged = performMixedModeNewtonIteration(
      context, sideState, foil.edge, side, stationIndex, previousTransition,
      ctx);
  if (!converged)
  {
    handleMixedModeNonConvergence(context, sideState, side, stationIndex, ctx);
  }

  sideState.sens = sideState.sennew;
  context.storeStationStateCommon(side, stationIndex, ctx);
  return true;
}

bool MarcherDu::performMixedModeNewtonIteration(
    MrchduContext &context, SideMarchState &sideState, const Edge &edge,
    int side, int ibl, int itrold,
    MrchduContext::MixedModeStationContext &ctx)
{
  bool converged = false;
  double ueref = 0.0;
  double hkref = 0.0;
  ctx.flowRegime = context.applyFlowRegimeCandidate(ctx.flowRegime);
  auto &state = context.mutableState();

  for (int itbl = 1; itbl <= 25; ++itbl)
  {
    blData updatedCurrent =
        context.blprv(state.current(), ctx.xsi, sideState.ami, ctx.cti, ctx.thi,
                      ctx.dsi, ctx.dswaki, ctx.uei);
    state.current() = updatedCurrent;
    context.blkin(state);
    ctx.flowRegime = context.currentFlowRegime();

    context.checkTransitionIfNeeded(side, ibl, ctx.isSimilarity(), 1,
                                    sideState.ami);
    ctx.flowRegime = context.currentFlowRegime();

    const bool startOfWake = context.isStartOfWake(side, ibl);
    context.updateSystemMatricesForStation(edge, side, ibl, ctx);

    if (itbl == 1)
    {
      context.initializeFirstIterationState(side, ibl, itrold, ctx, ueref, hkref);
    }

    if (ctx.isSimilarity() || startOfWake)
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

    if (context.applyMixedModeNewtonStep(side, ibl, sideState.ami, ctx))
    {
      converged = true;
      break;
    }
  }

  return converged;
}

void MarcherDu::handleMixedModeNonConvergence(
    MrchduContext &context, SideMarchState &sideState, int side, int ibl,
    MrchduContext::MixedModeStationContext &ctx)
{
  context.emitMarchFailureLog("mrchdu", side, ibl, ctx.dmax);
  context.recoverStationAfterFailure(
      side, ibl, ctx, sideState.ami,
      MrchduContext::EdgeVelocityFallbackMode::UsePreviousStation, 2);
}
