#pragma once

#include "simulation/march/base.hpp"

class MarcherDu : public Marcher
{
public:
  bool mrchdu(MrchduContext &context, const Foil &foil,
              const StagnationResult &stagnation);
  bool mrchdu(MrchduContext &context, BoundaryLayerState &state,
              const Foil &foil, const StagnationResult &stagnation);

private:
  struct SideMarchState
  {
    double sens = 0.0;
    double sennew = 0.0;
    double ami = 0.0;
  };

  const double senswt = 1000.0;
  bool marchBoundaryLayerSide(MrchduContext &context,
                              BoundaryLayerState &state, int side,
                              SideMarchState &sideState,
                              const Foil &foil,
                              const StagnationResult &stagnation);
  bool processBoundaryLayerStation(MrchduContext &context,
                                   BoundaryLayerState &state,
                                   SideMarchState &sideState, int side,
                                   int stationIndex, int previousTransition,
                                   const Foil &foil);
  bool performMixedModeNewtonIteration(
      MrchduContext &context, SideMarchState &sideState, const Edge &edge,
      int side, int ibl, int itrold,
      MrchduContext::MixedModeStationContext &ctx);
  void handleMixedModeNonConvergence(
      MrchduContext &context, SideMarchState &sideState, int side, int ibl,
      MrchduContext::MixedModeStationContext &ctx);
  MrchduContext::MixedModeStationContext
  prepareMixedModeStation(MrchduContext &context, SideMarchState &sideState,
                          int side,
                          int stationIndex, int previousTransition);
};
