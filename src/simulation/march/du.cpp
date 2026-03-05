#include "simulation/march/du.hpp"

#include <algorithm>
#include <cmath>

#include "XFoil.h"
#include "core/math_util.hpp"

using BoundaryContext = BoundaryLayerWorkflow::MixedModeStationContext;
using EdgeVelocityFallbackMode =
    BoundaryLayerWorkflow::EdgeVelocityFallbackMode;
using EdgeVelocityDistribution =
    BoundaryLayerWorkflow::EdgeVelocityDistribution;
using QtanResult = BoundaryLayerWorkflow::QtanResult;
using ClContributions = BoundaryLayerWorkflow::ClContributions;
using BoundaryLayerDelta = BoundaryLayerWorkflow::BoundaryLayerDelta;
using BoundaryLayerMetrics = BoundaryLayerWorkflow::BoundaryLayerMetrics;

BoundaryLayerWorkflow::MixedModeStationContext
MarcherDu::prepareMixedModeStation(BoundaryLayerWorkflow &workflow, int side,
                                   int stationIndex, int previousTransition)
{
  BoundaryContext ctx;
  const auto station = workflow.readStationModel(side, stationIndex);

  ctx.flowRegime = workflow.determineRegimeForStation(side, stationIndex);
  ctx.xsi = station.arcLength;
  ctx.uei = station.edgeVelocity;
  ctx.thi = station.momentumThickness;
  ctx.dsi = station.displacementThickness;

  if (stationIndex < previousTransition)
  {
    ami = station.skinFrictionCoeff;
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
  ctx.ami = ami;

  if (ctx.isWake())
  {
    ctx.dswaki = station.wakeGap;
  }
  else
  {
    ctx.dswaki = 0.0;
  }

  double thickness_limit =
      (stationIndex <= station.trailingEdgeIndex) ? 1.02 : 1.00005;
  ctx.dsi =
      std::max(ctx.dsi - ctx.dswaki, thickness_limit * ctx.thi) + ctx.dswaki;

  ctx.flowRegime = workflow.applyFlowRegimeCandidate(ctx.flowRegime);

  return ctx;
}

bool MarcherDu::mrchdu(BoundaryLayerWorkflow &workflow, const Foil &foil,
                       const StagnationResult &stagnation)
{
  return mrchdu(workflow, workflow.state, foil, stagnation);
}

bool MarcherDu::mrchdu(BoundaryLayerWorkflow &workflow,
                       BoundaryLayerState &state, const Foil &foil,
                       const StagnationResult &stagnation)
{

  for (int side = 1; side <= 2; ++side)
  {
    if (!marchBoundaryLayerSide(workflow, state, side, foil, stagnation))
    {
      return false;
    }
  }
  return true;
}

bool MarcherDu::marchBoundaryLayerSide(BoundaryLayerWorkflow &workflow,
                                       BoundaryLayerState &state, int side,
                                       const Foil &foil,
                                       const StagnationResult &stagnation)
{
  const int previousTransition = workflow.resetSideState(side, foil, stagnation);
  const int stationCount = workflow.readSideStationCount(side);

  for (int stationIndex = 0;
       stationIndex < stationCount - 1;
       ++stationIndex)
  {
    if (!processBoundaryLayerStation(workflow, state, side, stationIndex,
                                     previousTransition, foil))
    {
      return false;
    }
  }

  return true;
}

bool MarcherDu::processBoundaryLayerStation(
    BoundaryLayerWorkflow &workflow, BoundaryLayerState &state, int side,
    int stationIndex, int previousTransition, const Foil &foil)
{
  BoundaryContext ctx = prepareMixedModeStation(workflow, side, stationIndex,
                                                previousTransition);
  bool converged = performMixedModeNewtonIteration(
      workflow, foil.edge, side, stationIndex, previousTransition, ctx);
  if (!converged)
  {
    handleMixedModeNonConvergence(workflow, side, stationIndex, ctx);
  }

  sens = sennew;
  workflow.storeStationStateCommon(side, stationIndex, ctx);
  return true;
}

bool MarcherDu::performMixedModeNewtonIteration(
    BoundaryLayerWorkflow &workflow, const Edge &edge, int side, int ibl,
    int itrold, BoundaryLayerWorkflow::MixedModeStationContext &ctx)
{
  bool converged = false;
  double ueref = 0.0;
  double hkref = 0.0;
  ctx.flowRegime = workflow.applyFlowRegimeCandidate(ctx.flowRegime);

  for (int itbl = 1; itbl <= 25; ++itbl)
  {
    blData updatedCurrent =
        workflow.blprv(workflow.state.current(), ctx.xsi, ami, ctx.cti, ctx.thi,
                       ctx.dsi, ctx.dswaki, ctx.uei);
    workflow.state.current() = updatedCurrent;
    workflow.blkin(workflow.state);
    ctx.flowRegime = workflow.currentFlowRegime();

    workflow.checkTransitionIfNeeded(side, ibl, ctx.isSimilarity(), 1, ami);
    ctx.flowRegime = workflow.currentFlowRegime();

    const bool startOfWake = workflow.isStartOfWake(side, ibl);
    workflow.updateSystemMatricesForStation(edge, side, ibl, ctx);

    if (itbl == 1)
    {
      workflow.initializeFirstIterationState(side, ibl, itrold, ctx, ueref,
                                             hkref);
    }

    if (ctx.isSimilarity() || startOfWake)
    {
      workflow.configureSimilarityRow(ueref);
    }
    else
    {
      const bool resetSensitivity = (itbl <= 5);
      const bool averageSensitivity = (itbl > 5 && itbl <= 15);
      workflow.configureViscousRow(hkref, ueref, senswt, resetSensitivity,
                                   averageSensitivity, sens, sennew);
    }

    if (workflow.applyMixedModeNewtonStep(side, ibl, ami, ctx))
    {
      converged = true;
      break;
    }
  }

  return converged;
}

void MarcherDu::handleMixedModeNonConvergence(
    BoundaryLayerWorkflow &workflow, int side, int ibl,
    BoundaryLayerWorkflow::MixedModeStationContext &ctx)
{
  workflow.emitMarchFailureLog("mrchdu", side, ibl, ctx.dmax);
  workflow.recoverStationAfterFailure(
      side, ibl, ctx, ami,
      BoundaryLayerWorkflow::EdgeVelocityFallbackMode::UsePreviousStation, 2);
}
