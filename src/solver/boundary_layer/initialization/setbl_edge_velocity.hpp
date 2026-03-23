#pragma once

#include <Eigen/Core>

#include "numerics/side_pair.hpp"
#include "solver/boundary_layer/initialization/setbl_access.hpp"

namespace setbl_edge_velocity {

struct LeTeSensitivities {
    SidePair<Eigen::VectorXd> ule_m;
    SidePair<Eigen::VectorXd> ute_m;
};

struct EdgeVelocitySensitivityResult {
    SidePair<Eigen::VectorXd> edgeVelocity;
    SidePair<Eigen::VectorXd> outputEdgeVelocity;
    SidePair<int> jvte{0, 0};
    SidePair<double> dule{0.0, 0.0};
    SidePair<Eigen::VectorXd> ule_m;
    SidePair<Eigen::VectorXd> ute_m;
    SidePair<double> ule_a{0.0, 0.0};
};

class SetblEdgeVelocityOps {
  public:
    static LeTeSensitivities computeLeTeSensitivities(const BoundaryLayerSetblAccess &access, int nsys,
                                                      const Eigen::MatrixXd &dij) {
        LeTeSensitivities sensitivities;
        sensitivities.ule_m.top    = Eigen::VectorXd::Zero(nsys);
        sensitivities.ule_m.bottom = Eigen::VectorXd::Zero(nsys);
        sensitivities.ute_m.top    = Eigen::VectorXd::Zero(nsys);
        sensitivities.ute_m.bottom = Eigen::VectorXd::Zero(nsys);

        const int ile1 = access.leadingEdgePanelIndex(1);
        const int ile2 = access.leadingEdgePanelIndex(2);
        const int ite1 = access.trailingEdgePanelIndex(1);
        const int ite2 = access.trailingEdgePanelIndex(2);

        for (int js = 1; js <= 2; ++js) {
            for (int jbl = 0; jbl < access.stationCount(js) - 1; ++jbl) {
                const int panel_index     = access.stationToPanel(js, jbl);
                const int system_index    = access.stationToSystem(js, jbl);
                const double panel_factor = access.panelInfluenceFactor(js, jbl);
                sensitivities.ule_m.top[system_index] =
                    -access.leadingEdgePanelInfluenceFactor(1) * panel_factor * dij(ile1, panel_index);
                sensitivities.ule_m.bottom[system_index] =
                    -access.leadingEdgePanelInfluenceFactor(2) * panel_factor * dij(ile2, panel_index);
                sensitivities.ute_m.top[system_index] =
                    -access.trailingEdgePanelInfluenceFactor(1) * panel_factor * dij(ite1, panel_index);
                sensitivities.ute_m.bottom[system_index] =
                    -access.trailingEdgePanelInfluenceFactor(2) * panel_factor * dij(ite2, panel_index);
            }
        }

        return sensitivities;
    }

    static EdgeVelocitySensitivityResult prepare(const BoundaryLayerSetblAccess &access,
                                                 SidePairRef<const BoundaryLayerSideState> profiles,
                                                 const Eigen::MatrixXd &dij, int nsys) {
        EdgeVelocitySensitivityResult result;
        result.edgeVelocity              = access.ueset(dij);
        result.outputEdgeVelocity.top    = profiles.top.edgeVelocity;
        result.outputEdgeVelocity.bottom = profiles.bottom.edgeVelocity;
        result.jvte.top                  = access.trailingEdgeSystemIndex(1);
        result.jvte.bottom               = access.trailingEdgeSystemIndex(2);
        result.dule.top                  = result.outputEdgeVelocity.top[0] - result.edgeVelocity.top[0];
        result.dule.bottom               = result.outputEdgeVelocity.bottom[0] - result.edgeVelocity.bottom[0];

        const auto sensitivities = computeLeTeSensitivities(access, nsys, dij);
        result.ule_m             = sensitivities.ule_m;
        result.ute_m             = sensitivities.ute_m;
        result.ule_a.top         = access.inviscidEdgeVelocitySensitivityToAlpha(1, 0);
        result.ule_a.bottom      = access.inviscidEdgeVelocitySensitivityToAlpha(2, 0);
        return result;
    }
};

} // namespace setbl_edge_velocity
