#include "solver/boundary_layer/workflow/relaxation.hpp"

#include <algorithm>

#include "model/boundary_layer/physics.hpp"

namespace {

void applyRelaxationLimit(const Eigen::VectorXd &dn, double dhi, double dlo,
                          double &relaxation) {
  double max_pos = 0.0;
  double min_neg = 0.0;
  for (const double value : dn) {
    if (value > max_pos) {
      max_pos = value;
    }
    if (value < min_neg) {
      min_neg = value;
    }
  }
  if (max_pos > 0.0) {
    relaxation = std::min(relaxation, dhi / max_pos);
  }
  if (min_neg < 0.0) {
    relaxation = std::min(relaxation, dlo / min_neg);
  }
}

} // namespace

BoundaryLayerDelta
BoundaryLayerRelaxationOps::buildBoundaryLayerDelta(
    const BoundaryLayerLattice &lattice_side, const Eigen::VectorXd &unew_side,
    const Eigen::VectorXd &u_ac_side, double dac,
    const BoundaryLayerMatrix3x2dVector &vdel) {
  BoundaryLayerDelta delta;
  const int len = lattice_side.stationCount - 1;
  if (len <= 0) {
    return delta;
  }

  delta.dskinFrictionCoeff = Eigen::VectorXd(len);
  delta.dmomentumThickness = Eigen::VectorXd(len);
  delta.ddisplacementThickness = Eigen::VectorXd(len);
  delta.dedgeVelocity = Eigen::VectorXd(len);

  const auto iv = lattice_side.stationToSystem.segment(0, len);
  Eigen::VectorXd dmass(len);
  for (int j = 0; j < len; ++j) {
    const int idx = iv[j];
    delta.dskinFrictionCoeff[j] = vdel[idx](0, 0) - dac * vdel[idx](0, 1);
    delta.dmomentumThickness[j] = vdel[idx](1, 0) - dac * vdel[idx](1, 1);
    dmass[j] = vdel[idx](2, 0) - dac * vdel[idx](2, 1);
  }

  const Eigen::VectorXd edge_velocity_segment =
      lattice_side.profiles.edgeVelocity.head(len);
  const Eigen::VectorXd displacement_segment =
      lattice_side.profiles.displacementThickness.head(len);

  delta.dedgeVelocity =
      unew_side.head(len) + dac * u_ac_side.head(len) - edge_velocity_segment;
  delta.ddisplacementThickness =
      (dmass - displacement_segment.cwiseProduct(delta.dedgeVelocity))
          .cwiseQuotient(edge_velocity_segment);
  return delta;
}

BoundaryLayerMetrics
BoundaryLayerRelaxationOps::evaluateSegmentRelaxation(
    const BoundaryLayerSideProfiles &profiles,
    const BoundaryLayerDelta &delta, double dhi,
    double dlo, double &relaxation) {
  BoundaryLayerMetrics metrics;
  const int len = delta.dskinFrictionCoeff.size();
  if (len <= 0) {
    return metrics;
  }

  const Eigen::VectorXd skin_friction_segment =
      profiles.skinFrictionCoeff.head(len);
  const Eigen::VectorXd momentum_segment =
      profiles.momentumThickness.head(len);
  const Eigen::VectorXd displacement_segment =
      profiles.displacementThickness.head(len);

  Eigen::VectorXd dn1(len);
  const int transition_index = profiles.transitionIndex;
  for (int idx = 0; idx < len; ++idx) {
    dn1[idx] = idx < transition_index
                   ? delta.dskinFrictionCoeff[idx] / 10.0
                   : delta.dskinFrictionCoeff[idx] / skin_friction_segment[idx];
  }
  const Eigen::VectorXd dn2 =
      delta.dmomentumThickness.cwiseQuotient(momentum_segment);
  const Eigen::VectorXd dn3 =
      delta.ddisplacementThickness.cwiseQuotient(displacement_segment);
  const Eigen::VectorXd dn4 = delta.dedgeVelocity.array().abs() / 0.25;

  applyRelaxationLimit(dn1, dhi, dlo, relaxation);
  applyRelaxationLimit(dn2, dhi, dlo, relaxation);
  applyRelaxationLimit(dn3, dhi, dlo, relaxation);
  applyRelaxationLimit(dn4, dhi, dlo, relaxation);

  metrics.rmsContribution = (dn1.array().square() + dn2.array().square() +
                             dn3.array().square() + dn4.array().square())
                                .sum();
  metrics.maxChange = dn1.cwiseAbs().maxCoeff();
  metrics.maxChange = std::max(metrics.maxChange, dn2.cwiseAbs().maxCoeff());
  metrics.maxChange = std::max(metrics.maxChange, dn3.cwiseAbs().maxCoeff());
  metrics.maxChange = std::max(metrics.maxChange, dn4.cwiseAbs().maxCoeff());
  return metrics;
}

BoundaryLayerSideProfiles BoundaryLayerRelaxationOps::applyBoundaryLayerDelta(
    const BoundaryLayerLattice &lattice_side, const Eigen::VectorXd &wgap,
    const BoundaryLayerDelta &delta, double relaxation,
    double hstinv, double gamm1) {
  BoundaryLayerSideProfiles updated = lattice_side.profiles;
  const int len = delta.dskinFrictionCoeff.size();
  if (len <= 0) {
    return updated;
  }

  updated.skinFrictionCoeff.head(len) += relaxation * delta.dskinFrictionCoeff;
  updated.momentumThickness.head(len) += relaxation * delta.dmomentumThickness;
  updated.displacementThickness.head(len) +=
      relaxation * delta.ddisplacementThickness;
  updated.edgeVelocity.head(len) += relaxation * delta.dedgeVelocity;

  const int transition_index = std::max(0, lattice_side.profiles.transitionIndex);
  for (int idx = transition_index; idx < len; ++idx) {
    updated.skinFrictionCoeff[idx] = std::min(updated.skinFrictionCoeff[idx], 0.25);
  }

  for (int ibl = 0; ibl < len; ++ibl) {
    double dswaki = 0.0;
    if (ibl > lattice_side.trailingEdgeIndex) {
      const int wake_index = ibl - (lattice_side.trailingEdgeIndex + 1);
      dswaki = wgap[wake_index];
    }

    const double hklim = ibl <= lattice_side.trailingEdgeIndex ? 1.02 : 1.00005;
    const double edge_velocity = updated.edgeVelocity[ibl];
    const double edge_velocity_sq = edge_velocity * edge_velocity;
    const double denom = 1.0 - 0.5 * edge_velocity_sq * hstinv;
    const double msq = edge_velocity_sq * hstinv / (gamm1 * denom);
    double dsw = updated.displacementThickness[ibl] - dswaki;
    dsw = BoundaryLayerPhysics::adjustDisplacementForHkLimit(
        dsw, updated.momentumThickness[ibl], msq, hklim);
    updated.displacementThickness[ibl] = dsw + dswaki;
    updated.massFlux[ibl] = updated.displacementThickness[ibl] * edge_velocity;
  }

  return updated;
}
