#include "domain/foil/foil.hpp"

#include <algorithm>
#include <cmath>

#include "simulation/psi.hpp"

namespace {

Eigen::VectorXd buildWakeSpacing(double ds1, double smax, int nn) {
  const int nex = nn - 1;

  const double sigma = smax / ds1;
  const double rnex = static_cast<double>(nex);
  const double rni = 1.0 / rnex;

  const double aaa = rnex * (rnex - 1.0) * (rnex - 2.0) / 6.0;
  const double bbb = rnex * (rnex - 1.0) / 2.0;
  const double ccc = rnex - sigma;

  double disc = std::max(0.0, bbb * bbb - 4.0 * aaa * ccc);
  double ratio = 1.0;
  if (nex == 2) {
    ratio = -ccc / bbb + 1.0;
  } else {
    ratio = (-bbb + std::sqrt(disc)) / (2.0 * aaa) + 1.0;
  }

  for (int iter = 0; iter < 100; iter++) {
    const double sigman = (std::pow(ratio, static_cast<double>(nex)) - 1.0) / (ratio - 1.0);
    const double sigman_rni = std::pow(sigman, rni);
    const double sigma_rni = std::pow(sigma, rni);
    const double res = sigman_rni - sigma_rni;
    const double numerator = rnex * std::pow(ratio, static_cast<double>(nex - 1)) - sigman;
    const double denominator = std::pow(ratio, static_cast<double>(nex)) - 1.0;
    const double dresdr = rni * sigman_rni * numerator / denominator;

    const double dratio = -res / dresdr;
    ratio += dratio;

    if (std::fabs(dratio) < 1.0e-5) {
      break;
    }
  }

  Eigen::VectorXd spline_length(nn);
  spline_length[0] = 0.0;
  double ds = ds1;
  for (int i = 1; i < nn; i++) {
    spline_length[i] = spline_length[i - 1] + ds;
    ds *= ratio;
  }
  return spline_length;
}

Eigen::Vector2d computePsiGradient(const Foil& foil, int node_index,
                                   const Eigen::Vector2d& point,
                                   const Eigen::Matrix2Xd& gamu,
                                   const Eigen::Matrix2Xd& surface_vortex,
                                   double alfa, double qinf,
                                   Eigen::VectorXd& apanel) {
  Eigen::Vector2d psi;
  psi.x() = psilin(foil, node_index, point, {1.0, 0.0}, false,
                   gamu, surface_vortex, alfa, qinf,
                   apanel)
                .psi_ni;
  psi.y() = psilin(foil, node_index, point, {0.0, 1.0}, false,
                   gamu, surface_vortex, alfa, qinf,
                   apanel)
                .psi_ni;
  return psi;
}
}  // namespace

bool Foil::xyWake(int wake_point_count,
                  const Eigen::Matrix2Xd& gamu,
                  const Eigen::Matrix2Xd& surface_vortex, double alfa,
                  double qinf) {
  const int point_count = foil_shape.n;
  if (point_count < 2 || wake_point_count < 2) {
    return false;
  }
  const int total_nodes = point_count + wake_point_count;

  wake_shape.points = Eigen::Matrix2Xd::Zero(2, total_nodes);
  wake_shape.n = wake_point_count;
  wake_shape.normal_vector = Eigen::Matrix2Xd::Zero(2, total_nodes);
  wake_shape.spline_length = Eigen::VectorXd::Zero(total_nodes);
  wake_shape.angle_panel = Eigen::VectorXd::Zero(total_nodes);
  wake_shape.points.block(0, 0, 2, point_count) = foil_shape.points;
  wake_shape.normal_vector.block(0, 0, 2, point_count) = foil_shape.normal_vector;
  wake_shape.spline_length.head(point_count) = foil_shape.spline_length;
  wake_shape.angle_panel.head(point_count) = foil_shape.angle_panel;

  double ds1 = 0.5 * (foil_shape.spline_length[1] - foil_shape.spline_length[0] +
                      foil_shape.spline_length[point_count - 1] -
                      foil_shape.spline_length[point_count - 2]);
  const Eigen::VectorXd wake_spacing = buildWakeSpacing(ds1, edge.chord, wake_point_count);

  Eigen::Vector2d tangent_vector =
      (foil_shape.dpoints_ds.col(point_count - 1) -
       foil_shape.dpoints_ds.col(0))
          .normalized();
  wake_shape.normal_vector.col(point_count) =
      Eigen::Vector2d{tangent_vector.y(), -tangent_vector.x()};

  wake_shape.points.col(point_count) = edge.point_te + 0.0001 * tangent_vector;
  foil_shape.points.col(point_count) = wake_shape.points.col(point_count);
  wake_shape.spline_length[point_count] = wake_shape.spline_length[point_count - 1];

  Eigen::Vector2d psi = computePsiGradient(*this, point_count,
                                           wake_shape.points.col(point_count),
                                           gamu, surface_vortex, alfa, qinf,
                                           wake_shape.angle_panel);
  wake_shape.normal_vector.col(point_count + 1) = -psi.normalized();
  wake_shape.angle_panel[point_count] = std::atan2(psi.y(), psi.x());

  for (int i = point_count + 1; i < total_nodes; i++) {
    const double ds = wake_spacing[i - point_count] -
                      wake_spacing[i - point_count - 1];
    wake_shape.points.col(i).x() =
        wake_shape.points.col(i - 1).x() -
        ds * wake_shape.normal_vector.col(i - 1).y();
    wake_shape.points.col(i).y() =
        wake_shape.points.col(i - 1).y() +
        ds * wake_shape.normal_vector.col(i - 1).x();
    foil_shape.points.col(i) = wake_shape.points.col(i);
    wake_shape.spline_length[i] = wake_shape.spline_length[i - 1] + ds;
    if (i == total_nodes - 1) {
      break;
    }

    Eigen::Vector2d psi_next =
        computePsiGradient(*this, i, wake_shape.points.col(i), gamu,
                           surface_vortex, alfa, qinf, wake_shape.angle_panel);
    wake_shape.normal_vector.col(i + 1) = -psi_next.normalized();
    wake_shape.angle_panel[i] = std::atan2(psi_next.y(), psi_next.x());
  }

  return true;
}
