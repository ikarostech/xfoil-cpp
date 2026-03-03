#pragma once

#include "simulation/march/base.hpp"

class MarcherDu : public Marcher
{
public:
  bool mrchdu(BoundaryLayerWorkflow &workflow, const Foil &foil,
              const StagnationResult &stagnation);
  bool mrchdu(BoundaryLayerWorkflow &workflow, BoundaryLayerState &state,
              const Foil &foil, const StagnationResult &stagnation);

private:
  const double senswt = 1000.0;
  double sens = 0.0;
  double sennew = 0.0;
  double ami;
  bool marchBoundaryLayerSide(BoundaryLayerWorkflow &workflow,
                              BoundaryLayerState &state, int side,
                              const Foil &foil,
                              const StagnationResult &stagnation);
  bool processBoundaryLayerStation(BoundaryLayerWorkflow &workflow,
                                   BoundaryLayerState &state, int side,
                                   int stationIndex, int previousTransition,
                                   const Foil &foil);
  bool performMixedModeNewtonIteration(
      BoundaryLayerWorkflow &workflow, const Edge &edge, int side, int ibl,
      int itrold, BoundaryLayerWorkflow::MixedModeStationContext &ctx);
  void handleMixedModeNonConvergence(
      BoundaryLayerWorkflow &workflow, int side, int ibl,
      BoundaryLayerWorkflow::MixedModeStationContext &ctx);
  BoundaryLayerWorkflow::MixedModeStationContext
  prepareMixedModeStation(BoundaryLayerWorkflow &workflow, int side,
                          int stationIndex, int previousTransition);
};
