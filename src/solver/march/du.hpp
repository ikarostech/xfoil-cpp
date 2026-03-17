#pragma once

#include "solver/march/base.hpp"

class MarcherDu : public Marcher
{
public:
  bool mrchdu(MrchduContext &context, const Foil &foil,
              const StagnationResult &stagnation);
  bool mrchdu(MrchduContext &context, BoundaryLayerState &state,
              const Foil &foil, const StagnationResult &stagnation);

private:
  struct SideInput
  {
    int side = 0;
    int stationCount = 0;
    int previousTransition = 0;
  };
  struct StationInput
  {
    int side = 0;
    int stationIndex = 0;
    int previousTransition = 0;
    Edge edge;
    MarchContextTypes::StationReadModel stationModel;
    bool startOfWake = false;
  };
  struct SideMarchState
  {
    double sens = 0.0;
    double sennew = 0.0;
    double ami = 0.0;
  };
  using StationMarchResult =
      Marcher::StationMarchResult<MrchduContext::MixedModeStationContext>;

  const double senswt = 1000.0;
  SideInput makeSideInput(MrchduContext &context, int side, const Foil &foil,
                          const StagnationResult &stagnation) const;
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
  StationInput makeStationInput(MrchduContext &context, int side,
                                int stationIndex, int previousTransition,
                                const Foil &foil) const;
  StationMarchResult performMixedModeNewtonIteration(
      MrchduContext &context, SideMarchState &sideState,
      const StationInput &input,
      MrchduContext::MixedModeStationContext station);
  MrchduContext::MixedModeStationContext
  prepareMixedModeStation(MrchduContext &context, SideMarchState &sideState,
                          const StationInput &input);
};
