#pragma once

#include <Eigen/Core>

#include "model/boundary_layer/physics.hpp"
#include "model/boundary_layer/state.hpp"
#include "numerics/side_pair.hpp"
#include "solver/boundary_layer/initialization/setbl_access.hpp"

namespace setbl_station_update {

struct StationPrimaryVars {
  double xsi = 0.0;
  double uei = 0.0;
  double thi = 0.0;
  double dsi = 0.0;
  double dswaki = 0.0;
  double ami = 0.0;
  double cti = 0.0;
};

struct StationUpdateResult {
  BoundaryLayerState state;
  Eigen::VectorXd u_m2;
  Eigen::VectorXd d_m2;
  double u_a2 = 0.0;
  double d_a2 = 0.0;
  double due2 = 0.0;
  double dds2 = 0.0;
};

class SetblStationUpdateOps {
 public:
  static StationUpdateResult update(const BoundaryLayerSetblAccess &access,
                                    int side, int station, int iv,
                                    const StationPrimaryVars &vars,
                                    const SidePair<Eigen::VectorXd> &usav,
                                    const SidePair<BoundaryLayerSideProfiles> &profiles,
                                    const BoundaryLayerState &base_state,
                                    int system_size,
                                    const Eigen::MatrixXd &dij) {
    StationUpdateResult result;
    const double d2_m2 = 1.0 / vars.uei;
    const double d2_u2 = -vars.dsi / vars.uei;

    result.u_m2 = Eigen::VectorXd::Zero(system_size);
    result.d_m2 = Eigen::VectorXd::Zero(system_size);
    for (int js = 1; js <= 2; ++js) {
      for (int jbl = 0; jbl < access.stationCount(js) - 1; ++jbl) {
        const int jv = access.stationToSystem(js, jbl);
        result.u_m2[jv] =
            -access.panelInfluenceFactor(side, station) *
            access.panelInfluenceFactor(js, jbl) *
            dij(access.stationToPanel(side, station),
                access.stationToPanel(js, jbl));
        result.d_m2[jv] = d2_u2 * result.u_m2[jv];
      }
    }
    result.d_m2[iv] += d2_m2;

    result.u_a2 = access.inviscidEdgeVelocitySensitivityToAlpha(side, station);
    result.d_a2 = d2_u2 * result.u_a2;
    result.due2 = profiles.get(side).edgeVelocity[station] - usav.get(side)[station];
    result.dds2 = d2_u2 * result.due2;

    result.state = base_state;
    BoundaryLayerPhysics::refreshCurrentStation(
        result.state, access.blCompressibility(), access.blReynolds(),
        vars.xsi, vars.ami, vars.cti, vars.thi, vars.dsi, vars.dswaki,
        vars.uei);
    return result;
  }
};

} // namespace setbl_station_update
