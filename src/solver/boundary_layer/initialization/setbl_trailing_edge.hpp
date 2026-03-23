#pragma once

#include <Eigen/Core>

#include "model/boundary_layer/state.hpp"
#include "model/foil/edge.hpp"
#include "numerics/side_pair.hpp"
#include "solver/boundary_layer/initialization/setbl_access.hpp"

namespace setbl_te {

struct TeWakeCoefficients {
    double tte      = 0.0;
    double cte      = 0.0;
    double dte      = 0.0;
    double tte_tte1 = 0.0;
    double tte_tte2 = 0.0;
    double dte_mte1 = 0.0;
    double dte_ute1 = 0.0;
    double dte_mte2 = 0.0;
    double dte_ute2 = 0.0;
    double cte_cte1 = 0.0;
    double cte_cte2 = 0.0;
    double cte_tte1 = 0.0;
    double cte_tte2 = 0.0;
};

struct TeWakeUpdateResult {
    bool isStartOfWake = false;
    Eigen::VectorXd d_m1;
    double due1 = 0.0;
    double dds1 = 0.0;
    TeWakeCoefficients coeffs{};
};

struct TeWakeJacobianAdjustments {
    double vz[3][2]                = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
    Eigen::Matrix<double, 3, 2> vb = Eigen::Matrix<double, 3, 2>::Zero();
};

class SetblTrailingEdgeOps {
  public:
    static TeWakeUpdateResult computeWakeUpdate(const BoundaryLayerSetblAccess &access, int side, int station,
                                                const SidePair<Eigen::VectorXd> &usav,
                                                const SidePair<Eigen::VectorXd> &ute_m, const SidePair<int> &jvte,
                                                const Eigen::VectorXd &d_m1_template,
                                                const SidePair<BoundaryLayerSideState> &profiles, const Edge &edge) {
        TeWakeUpdateResult result;
        if (station != access.trailingEdgeIndex(side) + 1) {
            return result;
        }
        result.isStartOfWake = true;

        const int top_te    = access.topTrailingEdgeIndex();
        const int bottom_te = access.bottomTrailingEdgeIndex();

        result.coeffs.tte = profiles.get(1).momentumThickness[top_te] + profiles.get(2).momentumThickness[bottom_te];
        result.coeffs.dte = profiles.get(1).displacementThickness[top_te] +
                            profiles.get(2).displacementThickness[bottom_te] + edge.ante;
        result.coeffs.cte =
            (profiles.get(1).skinFrictionCoeff[top_te] * profiles.get(1).momentumThickness[top_te] +
             profiles.get(2).skinFrictionCoeff[bottom_te] * profiles.get(2).momentumThickness[bottom_te]) /
            result.coeffs.tte;

        result.coeffs.tte_tte1 = 1.0;
        result.coeffs.tte_tte2 = 1.0;
        result.coeffs.dte_mte1 = 1.0 / profiles.top.edgeVelocity[top_te];
        result.coeffs.dte_ute1 = -profiles.get(1).displacementThickness[top_te] / profiles.top.edgeVelocity[top_te];
        result.coeffs.dte_mte2 = 1.0 / profiles.bottom.edgeVelocity[bottom_te];
        result.coeffs.dte_ute2 =
            -profiles.get(2).displacementThickness[bottom_te] / profiles.bottom.edgeVelocity[bottom_te];
        result.coeffs.cte_cte1 = profiles.get(1).momentumThickness[top_te] / result.coeffs.tte;
        result.coeffs.cte_cte2 = profiles.get(2).momentumThickness[bottom_te] / result.coeffs.tte;
        result.coeffs.cte_tte1 = (profiles.get(1).skinFrictionCoeff[top_te] - result.coeffs.cte) / result.coeffs.tte;
        result.coeffs.cte_tte2 = (profiles.get(2).skinFrictionCoeff[bottom_te] - result.coeffs.cte) / result.coeffs.tte;

        result.d_m1 = d_m1_template;
        for (int js = 1; js <= 2; ++js) {
            for (int jbl = 0; jbl < access.stationCount(js) - 1; ++jbl) {
                const int jv    = access.stationToSystem(js, jbl);
                result.d_m1[jv] = result.coeffs.dte_ute1 * ute_m.get(1)[jv] + result.coeffs.dte_ute2 * ute_m.get(2)[jv];
            }
        }
        result.d_m1[jvte.get(1)] += result.coeffs.dte_mte1;
        result.d_m1[jvte.get(2)] += result.coeffs.dte_mte2;

        result.dds1 = result.coeffs.dte_ute1 * (profiles.top.edgeVelocity[top_te] - usav.top[top_te]) +
                      result.coeffs.dte_ute2 * (profiles.bottom.edgeVelocity[bottom_te] - usav.bottom[bottom_te]);

        return result;
    }

    static TeWakeJacobianAdjustments computeJacobianAdjustments(const BoundaryLayerSetblAccess &access,
                                                                const TeWakeCoefficients &coeffs) {
        const auto &blc = access.blcConst();
        TeWakeJacobianAdjustments result;
        result.vz[0][0] = blc.a1(0, 0) * coeffs.cte_cte1;
        result.vz[0][1] = blc.a1(0, 0) * coeffs.cte_tte1 + blc.a1(0, 1) * coeffs.tte_tte1;
        result.vb(0, 0) = blc.a1(0, 0) * coeffs.cte_cte2;
        result.vb(0, 1) = blc.a1(0, 0) * coeffs.cte_tte2 + blc.a1(0, 1) * coeffs.tte_tte2;

        result.vz[1][0] = blc.a1(1, 0) * coeffs.cte_cte1;
        result.vz[1][1] = blc.a1(1, 0) * coeffs.cte_tte1 + blc.a1(1, 1) * coeffs.tte_tte1;
        result.vb(1, 0) = blc.a1(1, 0) * coeffs.cte_cte2;
        result.vb(1, 1) = blc.a1(1, 0) * coeffs.cte_tte2 + blc.a1(1, 1) * coeffs.tte_tte2;

        result.vz[2][0] = blc.a1(2, 0) * coeffs.cte_cte1;
        result.vz[2][1] = blc.a1(2, 0) * coeffs.cte_tte1 + blc.a1(2, 1) * coeffs.tte_tte1;
        result.vb(2, 0) = blc.a1(2, 0) * coeffs.cte_cte2;
        result.vb(2, 1) = blc.a1(2, 0) * coeffs.cte_tte2 + blc.a1(2, 1) * coeffs.tte_tte2;
        return result;
    }
};

} // namespace setbl_te
