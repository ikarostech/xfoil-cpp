#pragma once

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
  using StationMarchResult = Marcher::StationMarchResult<MrchueStationContext>;

  bool mrchue(MrchueContext &context, const Foil &foil,
              const StagnationResult &stagnation);

private:
  bool marchMrchueSide(MrchueContext &context, int side,
                       const Foil &foil, const StagnationResult &stagnation);
  SideMarchState initializeMrchueSide(MrchueContext &context, int side);
  void prepareMrchueStationContext(MrchueContext &context, int side,
                                   int stationIndex,
                                   const SideMarchState &sideState,
                                   MrchueStationContext &ctx);
  StationMarchResult performMrchueNewtonLoop(MrchueContext &context, int side,
                                             int stationIndex,
                                             MrchueStationContext station,
                                             const Edge &edge);
  void updateSideStateFromStation(const MrchueStationContext &ctx,
                                  SideMarchState &sideState);
  void updateSideStateForTrailingEdge(MrchueContext &context, int side,
                                      const Foil &foil, int stationIndex,
                                      SideMarchState &sideState);
  void storeMrchueStationState(MrchueContext &context, int side,
                               int stationIndex,
                               const MrchueStationContext &ctx);
};
