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
  {
    bool simi = false;
    bool wake = false;
    bool direct = true;
    double xsi = 0.0;
    double uei = 0.0;
    double thi = 0.0;
    double dsi = 0.0;
    double ami = 0.0;
    double cti = 0.0;
    double dswaki = 0.0;
    double cte = 0.0;
    double dte = 0.0;
    double tte = 0.0;
    double dmax = 0.0;
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
                                  int stationIndex, MrchueStationContext &ctx,
                                  std::stringstream &ss);
  void storeMrchueStationState(BoundaryLayerWorkflow &workflow, int side,
                               int stationIndex,
                               const MrchueStationContext &ctx);
};
