#pragma once

#include "simulation/march/base.hpp"

class MarcherDu : public Marcher {
 public:
  bool mrchdu(BoundaryLayerWorkflow& workflow, const Foil& foil,
              const StagnationResult& stagnation);
  bool mrchdu(BoundaryLayerWorkflow& workflow, BoundaryLayerState& state,
              const Foil& foil, const StagnationResult& stagnation);

 private:
  bool marchBoundaryLayerSide(BoundaryLayerWorkflow& workflow,
                              BoundaryLayerState& state, int side,
                              double senswt, double& sens, double& sennew,
                              double& ami, const Foil& foil,
                              const StagnationResult& stagnation);
  bool processBoundaryLayerStation(BoundaryLayerWorkflow& workflow,
                                   BoundaryLayerState& state, int side,
                                   int stationIndex, int previousTransition,
                                   double senswt, double& sens,
                                   double& sennew, double& ami,
                                   const Foil& foil);
  bool performMixedModeNewtonIteration(
      BoundaryLayerWorkflow& workflow, const Edge& edge, int side, int ibl,
      int itrold, BoundaryLayerWorkflow::MixedModeStationContext& ctx,
      double senswt, double& sens, double& sennew, double& ami);
  void handleMixedModeNonConvergence(
      BoundaryLayerWorkflow& workflow, int side, int ibl,
      BoundaryLayerWorkflow::MixedModeStationContext& ctx, double& ami);
  BoundaryLayerWorkflow::MixedModeStationContext prepareMixedModeStation(
      BoundaryLayerWorkflow& workflow, int side, int stationIndex,
      int previousTransition, double& ami);
};
