#include "simulation/boundary_layer_aerodynamics.hpp"

#include "XFoil.h"
#include "core/math_util.hpp"

BoundaryLayerEdgeVelocityDistribution
BoundaryLayerAerodynamicsOps::computeNewUeDistribution(
    const SidePair<BoundaryLayerLattice> &lattice, const XFoil &xfoil,
    const BoundaryLayerMatrix3x2dVector &vdel) {
  BoundaryLayerEdgeVelocityDistribution distribution;
  distribution.unew.top = Eigen::VectorXd::Zero(lattice.top.stationCount);
  distribution.unew.bottom = Eigen::VectorXd::Zero(lattice.bottom.stationCount);
  distribution.u_ac.top = Eigen::VectorXd::Zero(lattice.top.stationCount);
  distribution.u_ac.bottom = Eigen::VectorXd::Zero(lattice.bottom.stationCount);

  for (int side = 1; side <= 2; ++side) {
    for (int station = 0; station < lattice.get(side).stationCount - 1; ++station) {
      const int panelIndex = lattice.get(side).stationToPanel[station];
      double dui = 0.0;
      double dui_ac = 0.0;
      for (int otherSide = 1; otherSide <= 2; ++otherSide) {
        for (int otherStation = 0;
             otherStation < lattice.get(otherSide).stationCount - 1;
             ++otherStation) {
          const int otherPanel = lattice.get(otherSide).stationToPanel[otherStation];
          const int systemIndex = lattice.get(otherSide).stationToSystem[otherStation];
          const double influence =
              -lattice.get(side).panelInfluenceFactor[station] *
              lattice.get(otherSide).panelInfluenceFactor[otherStation] *
              xfoil.aerodynamicCache.dij(panelIndex, otherPanel);
          dui += influence * (lattice.get(otherSide).profiles.massFlux[otherStation] +
                              vdel[systemIndex](2, 0));
          dui_ac += influence * (-vdel[systemIndex](2, 1));
        }
      }

      const double inviscidDerivative =
          xfoil.analysis_state_.controlByAlpha
              ? 0.0
              : lattice.get(side).inviscidEdgeVelocityMatrix(1, station);
      distribution.unew.get(side)[station] =
          lattice.get(side).inviscidEdgeVelocityMatrix(0, station) + dui;
      distribution.u_ac.get(side)[station] = inviscidDerivative + dui_ac;
    }
  }
  return distribution;
}

BoundaryLayerClContributions
BoundaryLayerAerodynamicsOps::computeClFromEdgeVelocityDistribution(
    const SidePair<BoundaryLayerLattice> &lattice, const XFoil &xfoil,
    const BoundaryLayerEdgeVelocityDistribution &distribution) {
  BoundaryLayerClContributions contributions;
  const int point_count = xfoil.foil.foil_shape.n;
  if (point_count == 0) {
    return contributions;
  }

  Eigen::VectorXd qnew = Eigen::VectorXd::Zero(point_count);
  Eigen::VectorXd q_ac = Eigen::VectorXd::Zero(point_count);
  for (int side = 1; side <= 2; ++side) {
    const Eigen::VectorXd &unew_vec = distribution.unew.get(side);
    const Eigen::VectorXd &uac_vec = distribution.u_ac.get(side);
    const int limit = lattice.get(side).trailingEdgeIndex;
    for (int station = 0; station < limit; ++station) {
      const int panelIndex = lattice.get(side).stationToPanel[station];
      qnew[panelIndex] = lattice.get(side).panelInfluenceFactor[station] * unew_vec[station];
      q_ac[panelIndex] = lattice.get(side).panelInfluenceFactor[station] * uac_vec[station];
    }
  }

  const auto compressibility = xfoil.buildCompressibilityParams();
  const auto cp_first =
      xfoil.computePressureCoefficient(qnew[0], q_ac[0], compressibility);

  double cpg1 = cp_first.cp;
  double cpg1_ms = cp_first.cp_msq;
  double cpg1_ac = cp_first.cp_velocity_derivative;

  for (int i = 0; i < point_count; i++) {
    const int ip = (i + 1) % point_count;
    const auto cp_next =
        xfoil.computePressureCoefficient(qnew[ip], q_ac[ip], compressibility);

    const double cpg2 = cp_next.cp;
    const double cpg2_ms = cp_next.cp_msq;
    const double cpg2_ac = cp_next.cp_velocity_derivative;

    const Eigen::Vector2d dpoint =
        MathUtil::getRotateMatrix(xfoil.analysis_state_.alpha) *
        (xfoil.foil.foil_shape.points.col(ip) -
         xfoil.foil.foil_shape.points.col(i));

    const double ag = 0.5 * (cpg2 + cpg1);
    const double ag_ms = 0.5 * (cpg2_ms + cpg1_ms);
    const double ag_ac = 0.5 * (cpg2_ac + cpg1_ac);

    contributions.cl += dpoint.x() * ag;
    contributions.cl_a += dpoint.y() * ag;
    contributions.cl_ms += dpoint.x() * ag_ms;
    contributions.cl_ac += dpoint.x() * ag_ac;

    cpg1 = cpg2;
    cpg1_ms = cpg2_ms;
    cpg1_ac = cpg2_ac;
  }

  return contributions;
}

SidePair<Eigen::VectorXd> BoundaryLayerAerodynamicsOps::ueset(
    const SidePair<BoundaryLayerLattice> &lattice, const Eigen::MatrixXd &dij) {
  SidePair<Eigen::VectorXd> edge_velocity;
  edge_velocity.top = lattice.top.profiles.edgeVelocity;
  edge_velocity.bottom = lattice.bottom.profiles.edgeVelocity;
  auto computeInducedVelocity = [&](int side, int station) {
    double dui = 0.0;
    const double side_panel_factor = lattice.get(side).panelInfluenceFactor[station];
    const int side_panel_index = lattice.get(side).stationToPanel[station];

    for (int js = 1; js <= 2; ++js) {
      for (int jbl = 0; jbl < lattice.get(js).stationCount - 1; ++jbl) {
        const double ue_m =
            -side_panel_factor * lattice.get(js).panelInfluenceFactor[jbl] *
            dij(side_panel_index, lattice.get(js).stationToPanel[jbl]);
        dui += ue_m * lattice.get(js).profiles.massFlux[jbl];
      }
    }
    return dui;
  };

  for (int side = 1; side <= 2; side++) {
    for (int station = 0; station < lattice.get(side).stationCount - 1; ++station) {
      const double dui = computeInducedVelocity(side, station);
      edge_velocity.get(side)[station] =
          lattice.get(side).inviscidEdgeVelocityMatrix(0, station) + dui;
    }
  }
  return edge_velocity;
}
