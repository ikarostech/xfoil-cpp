#pragma once

#include <sstream>

#include "simulation/march/base.hpp"

class MarcherUe : public Marcher
{
private:
  struct SideMarchState
  {
    double thi = 0.0;
    double dsi = 0.0;
    double ami = 0.0;
    double cti = 0.0;
  };

public:
  struct MrchueStationContext
      : public MrchueContext::MixedModeStationContext
  {
    bool direct = true;
    double hmax = 0.0;
    double htarg = 0.0;
  };

  bool mrchue(MrchueContext &context, const Foil &foil,
              const StagnationResult &stagnation);

private:
  bool marchMrchueSide(MrchueContext &context, int side,
                       const Foil &foil, const StagnationResult &stagnation,
                       std::stringstream &ss);
  SideMarchState initializeMrchueSide(MrchueContext &context, int side);
  void prepareMrchueStationContext(MrchueContext &context, int side,
                                   int stationIndex,
                                   const SideMarchState &sideState,
                                   MrchueStationContext &ctx);
  bool performMrchueNewtonLoop(MrchueContext &context, int side,
                               int stationIndex, MrchueStationContext &ctx,
                               const Edge &edge, std::stringstream &ss);
  void updateSideStateFromStation(const MrchueStationContext &ctx,
                                  SideMarchState &sideState);
  void updateSideStateForTrailingEdge(MrchueContext &context, int side,
                                      const Foil &foil, int stationIndex,
                                      SideMarchState &sideState);
  void handleMrchueStationFailure(MrchueContext &context, int side,
                                  int stationIndex, MrchueStationContext &ctx);
  void storeMrchueStationState(MrchueContext &context, int side,
                               int stationIndex,
                               const MrchueStationContext &ctx);
};
