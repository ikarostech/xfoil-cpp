#pragma once

#include <sstream>

#include "simulation/march/base.hpp"

class MarcherUe : public Marcher
{
private:
  double thi = 0.0;
  double dsi = 0.0;
  double ami = 0.0;
  double cti = 0.0;

public:
  struct MrchueStationContext
      : public BoundaryLayerWorkflow::MixedModeStationContext
  {
    bool direct = true;
    double hmax = 0.0;
    double htarg = 0.0;
  };

  bool mrchue(BoundaryLayerWorkflow &workflow, const Foil &foil,
              const StagnationResult &stagnation);

private:
  bool marchMrchueSide(BoundaryLayerWorkflow &workflow, int side,
                       const Foil &foil, const StagnationResult &stagnation,
                       std::stringstream &ss);
  void initializeMrchueSide(BoundaryLayerWorkflow &workflow, int side);
  void prepareMrchueStationContext(BoundaryLayerWorkflow &workflow, int side,
                                   int stationIndex, MrchueStationContext &ctx);
  bool performMrchueNewtonLoop(BoundaryLayerWorkflow &workflow, int side,
                               int stationIndex, MrchueStationContext &ctx,
                               const Edge &edge, std::stringstream &ss);
  void handleMrchueStationFailure(BoundaryLayerWorkflow &workflow, int side,
                                  int stationIndex, MrchueStationContext &ctx);
  void storeMrchueStationState(BoundaryLayerWorkflow &workflow, int side,
                               int stationIndex,
                               const MrchueStationContext &ctx);
};
