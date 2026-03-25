#pragma once

#include "solver/march/base.hpp"

class MarcherUe : public Marcher
{
private:
  struct SideInput
  {
    int side = 0;
    int stationCount = 0;
    MarchContextTypes::StationReadModel leadingStationModel;
    double reybl = 0.0;
  };
  struct SideMarchState
  {
    double thi = 0.0;
    double dsi = 0.0;
    double ami = 0.0;
    double cti = 0.0;
  };

public:
  struct StationInput
  {
    int side = 0;
    int stationIndex = 0;
    Edge edge;
    MarchContextTypes::StationReadModel stationModel;
    MarchContextTypes::TrailingEdgeReadModel trailingEdgeModel;
    double hstinv = 0.0;
    double gm1bl = 0.0;
  };
  struct MrchueStationContext
      : public MrchueContext::MixedModeStationContext
  {
    bool direct = true;
    double hmax = 0.0;
    double htarg = 0.0;
  };
  using StationMarchResult = Marcher::StationMarchResult<MrchueStationContext>;

  bool mrchue(MrchueContext &context, const Foil &foil,
              const StagnationFeature &stagnation);

private:
  bool marchMrchueSide(MrchueContext &context, int side, const Foil &foil,
                       const StagnationFeature &stagnation);
  SideInput makeSideInput(MrchueContext &context, int side, const Foil &foil,
                          const StagnationFeature &stagnation) const;
  SideMarchState initializeMrchueSide(const SideInput &input);
  void prepareMrchueStationContext(MrchueContext &context,
                                   const StationInput &input,
                                   const SideMarchState &sideState,
                                   MrchueStationContext &ctx);
  StationInput makeStationInput(MrchueContext &context, int side,
                                int stationIndex, const Foil &foil) const;
  StationMarchResult performMrchueNewtonLoop(MrchueContext &context,
                                             const StationInput &input,
                                             MrchueStationContext station,
                                             const Edge &edge);
  void updateSideStateFromStation(const MrchueStationContext &ctx,
                                  SideMarchState &sideState);
  void updateSideStateForTrailingEdge(const StationInput &input,
                                      const Foil &foil,
                                      SideMarchState &sideState);
  void storeMrchueStationState(MrchueContext &context, int side,
                               int stationIndex,
                               const MrchueStationContext &ctx);
};
