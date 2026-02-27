#include "XFoil.h"
#include "Eigen/Core"
#include <algorithm>
#include <numbers>

using Eigen::Matrix2Xd;
using Eigen::Vector2d;
using Eigen::VectorXd;

bool XFoil::abcopy(Matrix2Xd copyFrom) {
  constexpr double kPointMergeTolerance = 1.0e-14;
  int point_count = static_cast<int>(copyFrom.cols());

  //---- strip out doubled points
  int r = 1;
  while (r < point_count) {
    const double delta_norm =
        (copyFrom.col(r - 1) - copyFrom.col(r)).norm();
    if (delta_norm <= kPointMergeTolerance) {
      for (int j = r; j < point_count - 1; j++) {
        copyFrom.col(j) = copyFrom.col(j + 1);
      }
      point_count -= 1;
    } else {
      r++;
    }
  }
  //--- number of wake points
  int wake_point_count = point_count / 8 + 2;
  foil.wake_shape.n = wake_point_count;
  Matrix2Xd foil_points = Matrix2Xd::Zero(2, point_count + wake_point_count);
  foil_points.leftCols(point_count) = copyFrom.leftCols(point_count);

  foil.foil_shape.n = point_count;
  initialize();  
  
  foil = Foil(foil_points, point_count);
  foil.wake_shape.n = wake_point_count;
  updateTrailingEdgeState();

  invalidateWakeGeometry();
  invalidatePanelMap();
  setBLInitialized(false);
  invalidateConvergedSolution();

  return true;
}

double XFoil::atanc(double y, double x, double thold) {
  double tpi, thnew, dmomentumThickness, dtcorr;
  tpi = 6.2831853071795864769;
  thnew = atan2(y, x);
  dmomentumThickness = thnew - thold;
  dtcorr = dmomentumThickness - tpi * int((dmomentumThickness + sign(std::numbers::pi, dmomentumThickness)) / tpi);
  return thold + dtcorr;
}

void XFoil::updateTrailingEdgeState() {
  const int node_count = foil.foil_shape.n;

  if (node_count < 2) {
    sigte = 0.0;
    gamte = 0.0;
    return;
  }

  double scs = 0.0;
  double sds = 0.0;
  if (foil.edge.sharp) {
    scs = 1.0;
    sds = 0.0;
  } else if (foil.edge.dste != 0.0) {
    const double inv_dste = 1.0 / foil.edge.dste;
    scs = foil.edge.ante * inv_dste;
    sds = foil.edge.aste * inv_dste;
  }

  if (surface_vortex.rows() > 0 && surface_vortex.cols() >= node_count) {
    const double surface_delta =
        surface_vortex(0, 0) - surface_vortex(0, node_count - 1);
    sigte = 0.5 * surface_delta * scs;
    gamte = -0.5 * surface_delta * sds;
  } else {
    sigte = 0.0;
    gamte = 0.0;
  }
}

bool XFoil::xyWake() {
  if (!foil.xyWake(foil.wake_shape.n, aerodynamicCache.gamu, surface_vortex, analysis_state_.alpha,
                   analysis_state_.qinf)) {
    return false;
  }
  return true;
}

double XFoil::sign(double a, double b) {
  if (b >= 0.0)
    return fabs(a);
  else
    return -fabs(a);
}
