#include "simulation/march/du.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "XFoil.h"
#include "core/math_util.hpp"
#include "infrastructure/logger.hpp"

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

  ctx.flowRegime = workflow.determineRegimeForStation(side, stationIndex);
  ctx.xsi = workflow.lattice.get(side).arcLengthCoordinates[stationIndex];
  ctx.uei = workflow.lattice.get(side).profiles.edgeVelocity[stationIndex];
  ctx.thi = workflow.lattice.get(side).profiles.momentumThickness[stationIndex];
  ctx.dsi =
      workflow.lattice.get(side).profiles.displacementThickness[stationIndex];

  if (stationIndex < previousTransition)
  {
    ami = workflow.lattice.get(side).profiles.skinFrictionCoeff[stationIndex];
    ctx.cti = 0.03;
  }
  else
  {
    ctx.cti =
        workflow.lattice.get(side).profiles.skinFrictionCoeff[stationIndex];
    if (ctx.cti <= 0.0)
    {
      ctx.cti = 0.03;
    }
  }
  ctx.ami = ami;

  if (ctx.isWake())
  {
    int iw = stationIndex - workflow.lattice.get(side).trailingEdgeIndex;
    ctx.dswaki = workflow.wgap[iw - 1];
  }
  else
  {
    ctx.dswaki = 0.0;
  }

  double thickness_limit =
      (stationIndex <= workflow.lattice.get(side).trailingEdgeIndex) ? 1.02
                                                                     : 1.00005;
  ctx.dsi =
      std::max(ctx.dsi - ctx.dswaki, thickness_limit * ctx.thi) + ctx.dswaki;

  workflow.flowRegime = ctx.flowRegime;

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

  for (int stationIndex = 0;
       stationIndex < workflow.lattice.get(side).stationCount - 1;
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

  for (int itbl = 1; itbl <= 25; ++itbl)
  {
    blData updatedCurrent =
        workflow.blprv(workflow.state.current(), ctx.xsi, ami, ctx.cti, ctx.thi,
                       ctx.dsi, ctx.dswaki, ctx.uei);
    workflow.state.current() = updatedCurrent;
    workflow.blkin(workflow.state);

    workflow.checkTransitionIfNeeded(side, ibl, ctx.isSimilarity(), 1, ami);

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
  std::stringstream ss;
  ss << "     mrchdu: convergence failed at " << ibl << " ,  side " << side
     << ", res=" << std::setw(4) << std::fixed << std::setprecision(3)
     << ctx.dmax << "\n";
  Logger::instance().write(ss.str());

  workflow.resetStationKinematicsAfterFailure(
      side, ibl, ctx,
      BoundaryLayerWorkflow::EdgeVelocityFallbackMode::UsePreviousStation);

  blData updatedCurrent =
      workflow.blprv(workflow.state.current(), ctx.xsi, ami, ctx.cti, ctx.thi,
                     ctx.dsi, ctx.dswaki, ctx.uei);
  workflow.state.current() = updatedCurrent;
  workflow.blkin(workflow.state);

  workflow.checkTransitionIfNeeded(side, ibl, ctx.isSimilarity(), 2, ami);

  workflow.syncStationRegimeStates(side, ibl, ctx.flowRegime);

  ctx.ami = ami;
}
