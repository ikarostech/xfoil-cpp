#pragma once

#include <array>

#include "solver/boundary_layer/boundary_layer_builder.hpp"
#include "solver/boundary_layer/initialization/setbl_access.hpp"
#include "solver/boundary_layer/initialization/setbl_jacobian.hpp"
#include "solver/boundary_layer/boundary_layer_geometry.hpp"

namespace setbl_side_sweep {

class SetblSideSweepOps {
 public:
  static void assignXiUle(
      int side, const StagnationResult &stagnation,
      std::array<setbl_jacobian::StationState, 2> &stations) {
    if (side == 1) {
      stations[0].xi_ule = stagnation.sst_go;
      stations[1].xi_ule = -stagnation.sst_gp;
      return;
    }

    stations[0].xi_ule = -stagnation.sst_go;
    stations[1].xi_ule = stagnation.sst_gp;
  }

  static void updateRegimeAfterStation(const BoundaryLayerSetblAccess &access,
                                       int side, int station,
                                       SetblOutputView &output) {
    if (output.flowRegime == FlowRegimeEnum::Transition) {
      output.profiles.get(side).transitionIndex = station;
      access.setTransitionIndex(side, station);
      output.flowRegime = FlowRegimeEnum::Turbulent;
      access.flowRegime() = output.flowRegime;
    }

    if (station == access.trailingEdgeIndex(side)) {
      output.flowRegime = FlowRegimeEnum::Wake;
      access.flowRegime() = output.flowRegime;
      access.solveWakeState();
      access.blmid(FlowRegimeEnum::Wake);
    }
  }

  static void advanceStationWindow(
      const setbl_jacobian::StationState &next,
      setbl_jacobian::StationState &current,
      const BoundaryLayerSetblAccess &access) {
    current.u_m = next.u_m;
    current.d_m = next.d_m;
    current.u_a = next.u_a;
    current.d_a = next.d_a;
    current.due = next.due;
    current.dds = next.dds;
    access.advanceState();
  }
};

} // namespace setbl_side_sweep
