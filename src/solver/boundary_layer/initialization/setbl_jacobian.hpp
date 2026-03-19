#pragma once

#include <array>

#include <Eigen/Core>

#include "numerics/side_pair.hpp"
#include "solver/boundary_layer/blsolve.hpp"
#include "solver/boundary_layer/initialization/setbl_access.hpp"

namespace setbl_jacobian {

struct StationState {
  Eigen::VectorXd u_m;
  Eigen::VectorXd d_m;
  double u_a = 0.0;
  double d_a = 0.0;
  double due = 0.0;
  double dds = 0.0;
  double xi_ule = 0.0;
};

struct SideState {
  SidePair<Eigen::VectorXd> ule_m;
  SidePair<double> ule_a{0.0, 0.0};
  SidePair<double> dule{0.0, 0.0};
};

class SetblJacobianOps {
 public:
  static void assembleForStation(
      const BoundaryLayerSetblAccess &access, int iv, int nsys,
      const std::array<StationState, 2> &stations, const SideState &sideState,
      bool controlByAlpha, double re_clmr, double msq_clmr,
      Blsolve::BlNewtonSystem &system) {
    const auto &blc = access.blcConst();
    system.vb[iv] = blc.a1.block(0, 0, 3, 2);
    system.va[iv] = blc.a2.block(0, 0, 3, 2);

    Eigen::Matrix<double, 3, 4> A;
    A.col(0) = blc.a1.col(3).head<3>();
    A.col(1) = blc.a1.col(2).head<3>();
    A.col(2) = blc.a2.col(3).head<3>();
    A.col(3) = blc.a2.col(2).head<3>();

    Eigen::Matrix<double, 4, 2> B;
    B << stations[0].due, stations[0].u_a, stations[0].dds, stations[0].d_a,
        stations[1].due, stations[1].u_a, stations[1].dds, stations[1].d_a;
    const Eigen::Vector3d ax = (blc.a1.col(4) + blc.a2.col(4) + blc.d_xi).head<3>();
    const Eigen::RowVector2d xi =
        stations[0].xi_ule *
            Eigen::RowVector2d(sideState.dule.get(1), sideState.ule_a.get(1)) +
        stations[1].xi_ule *
            Eigen::RowVector2d(sideState.dule.get(2), sideState.ule_a.get(2));
    system.vdel[iv] = A * B + ax * xi;

    for (int jv = 1; jv < nsys; ++jv) {
      const Eigen::Vector4d m(stations[0].u_m(jv), stations[0].d_m(jv),
                              stations[1].u_m(jv), stations[1].d_m(jv));
      const double xi_m =
          stations[0].xi_ule * sideState.ule_m.get(1)(jv) +
          stations[1].xi_ule * sideState.ule_m.get(2)(jv);
      const Eigen::Vector3d vm = A * m + ax * xi_m;
      system.vm.at(0, jv, iv) = vm[0];
      system.vm.at(1, jv, iv) = vm[1];
      system.vm.at(2, jv, iv) = vm[2];
    }

    if (controlByAlpha) {
      system.vdel[iv].col(1).head<3>() =
          blc.d_re.head(3) * re_clmr + blc.d_msq.head(3) * msq_clmr;
    }
  }
};

} // namespace setbl_jacobian
