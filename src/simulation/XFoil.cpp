/****************************************************************************

        XFoil Class
        Copyright (C) 2000 Mark Drela
        Copyright (C) 2003 Andre Deperrois techwinder@gmail.com

        This program is free software; you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation; either version 2 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program; if not, write to the Free Software
        Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

*****************************************************************************/

#include "XFoil.h"
#include "domain/boundary_layer/boundary_layer_diff_solver.hpp"
#include "Eigen/Core"
#include <cmath>
#include <limits>
using namespace Eigen;

void ClearInitState(const XFoil& xfoil);

XFoil::CompressibilityParams XFoil::buildCompressibilityParams() const {
  const double current_mach = analysis_state_.currentMach;
  const double beta = std::sqrt(1.0 - current_mach * current_mach);
  const double beta_msq = -0.5 / beta;
  const double prandtlGlauertFactor =
      0.5 * current_mach * current_mach / (1.0 + beta);
  const double prandtlGlauertFactor_msq =
      0.5 / (1.0 + beta) -
      prandtlGlauertFactor / (1.0 + beta) * beta_msq;
  const double karmanTsienFactor =
      (current_mach / (1.0 + beta)) * (current_mach / (1.0 + beta));
  const double karmanTsienFactor_msq =
      1.0 / ((1.0 + beta) * (1.0 + beta)) -
      2.0 * karmanTsienFactor / (1.0 + beta) * beta_msq;
  return {beta, beta_msq, karmanTsienFactor, karmanTsienFactor_msq,
          prandtlGlauertFactor, prandtlGlauertFactor_msq};
}

XFoil::PressureCoefficientResult XFoil::computePressureCoefficient(
    double tangential_velocity, double velocity_derivative,
    const CompressibilityParams &params) const {
  const double velocity_ratio = tangential_velocity / analysis_state_.qinf;
  const double cginc = 1.0 - velocity_ratio * velocity_ratio;
  const double denom = params.beta + params.prandtlGlauertFactor * cginc;
  const double pressure_coefficient = cginc / denom;
  const double cp_msq =
      -pressure_coefficient / denom *
      (params.beta_msq + params.prandtlGlauertFactor_msq * cginc);

  double cp_velocity_derivative = 0.0;
  if (velocity_derivative != 0.0) {
    const double freestream_speed = analysis_state_.qinf;
    const double cpi =
        -2.0 * tangential_velocity / (freestream_speed * freestream_speed);
    const double cpc_cpi =
        (1.0 - params.prandtlGlauertFactor * pressure_coefficient) / denom;
    cp_velocity_derivative = cpc_cpi * cpi * velocity_derivative;
  }

  return {pressure_coefficient, cp_msq, cp_velocity_derivative};
}

XFoil::XFoil() : analysis_state_() {
  // fortran seems to initializes variables to 0
  mvisc = 0.0;

  // initialize transition parameters until user changes them
  acrit = 9.0;
  auto& boundary_layer_lattice = boundaryLayerWorkflow.lattice;
  boundary_layer_lattice.top.transitionLocation = 1.0;
  boundary_layer_lattice.bottom.transitionLocation = 1.0;

  //---- initialize freestream mach number to zero
  analysis_state_.machType = MachType::CONSTANT;
  analysis_state_.referenceMach = 0.0;

  //---- drop tolerance for bl system solver
  setVAccel(0.01);
  //---- default viscous parameters
  analysis_state_.reynoldsType = ReynoldsType::CONSTANT;
  analysis_state_.referenceRe = 0.0;
}

XFoil::~XFoil() {
  ClearInitState(*this);
}

bool XFoil::isBLInitialized() const {
  if (!hasPanelMap()) {
    return false;
  }
  const auto& top = boundaryLayerWorkflow.lattice.top;
  const auto& bottom = boundaryLayerWorkflow.lattice.bottom;
  if (top.stationCount <= 1 || bottom.stationCount <= 1) {
    return false;
  }
  if (top.profiles.edgeVelocity.size() < top.stationCount ||
      bottom.profiles.edgeVelocity.size() < bottom.stationCount ||
      top.profiles.skinFrictionCoeff.size() < top.stationCount ||
      bottom.profiles.skinFrictionCoeff.size() < bottom.stationCount ||
      top.profiles.momentumThickness.size() < top.stationCount ||
      bottom.profiles.momentumThickness.size() < bottom.stationCount ||
      top.profiles.displacementThickness.size() < top.stationCount ||
      bottom.profiles.displacementThickness.size() < bottom.stationCount) {
    return false;
  }

  const int top_count = top.stationCount - 1;
  const int bottom_count = bottom.stationCount - 1;
  const auto top_theta = top.profiles.momentumThickness.head(top_count);
  const auto bottom_theta = bottom.profiles.momentumThickness.head(bottom_count);
  const auto top_delta = top.profiles.displacementThickness.head(top_count);
  const auto bottom_delta = bottom.profiles.displacementThickness.head(bottom_count);

  if (!top_theta.allFinite() || !bottom_theta.allFinite() ||
      !top_delta.allFinite() || !bottom_delta.allFinite()) {
    return false;
  }

  return top_theta.maxCoeff() > 0.0 && bottom_theta.maxCoeff() > 0.0 &&
         top_delta.maxCoeff() > 0.0 && bottom_delta.maxCoeff() > 0.0;
}

void XFoil::setBLInitialized(bool bInitialized) {
  if (bInitialized) {
    return;
  }
  auto invalidateSide = [](BoundaryLayerLattice& lattice) {
    if (lattice.profiles.edgeVelocity.size() > 0)
      lattice.profiles.edgeVelocity.setZero();
    if (lattice.profiles.skinFrictionCoeff.size() > 0)
      lattice.profiles.skinFrictionCoeff.setZero();
    if (lattice.profiles.momentumThickness.size() > 0)
      lattice.profiles.momentumThickness.setZero();
    if (lattice.profiles.displacementThickness.size() > 0)
      lattice.profiles.displacementThickness.setZero();
    if (lattice.profiles.massFlux.size() > 0)
      lattice.profiles.massFlux.setZero();
  };
  invalidateSide(boundaryLayerWorkflow.lattice.top);
  invalidateSide(boundaryLayerWorkflow.lattice.bottom);
  invalidateConvergedSolution();
}

void XFoil::invalidateConvergedSolution() {
  avisc = std::numeric_limits<double>::quiet_NaN();
  mvisc = std::numeric_limits<double>::quiet_NaN();
}

void XFoil::invalidateWakeGeometry() {
  const int point_count = foil.foil_shape.n;
  const int wake_point_count = foil.wake_shape.n;
  const int total_nodes = point_count + wake_point_count;
  if (wake_point_count > 0 && aerodynamicCache.dij.rows() >= total_nodes &&
      aerodynamicCache.dij.cols() >= total_nodes) {
    aerodynamicCache.dij.block(point_count, point_count,
                               wake_point_count, wake_point_count)
        .setConstant(std::numeric_limits<double>::quiet_NaN());
  }
  foil.wake_shape.points.resize(0, 0);
  foil.wake_shape.normal_vector.resize(0, 0);
  foil.wake_shape.spline_length.resize(0);
  foil.wake_shape.angle_panel.resize(0);
}

void XFoil::invalidatePanelMap() {
  boundaryLayerWorkflow.lattice.top.stationCount = 0;
  boundaryLayerWorkflow.lattice.bottom.stationCount = 0;
}

bool XFoil::hasWakeGeometry() const {
  const int point_count = foil.foil_shape.n;
  const int wake_point_count = foil.wake_shape.n;
  const int total_nodes = point_count + wake_point_count;
  if (wake_point_count < 2) {
    return false;
  }
  if (foil.wake_shape.points.cols() < total_nodes ||
      foil.wake_shape.normal_vector.cols() < total_nodes ||
      foil.wake_shape.spline_length.size() < total_nodes ||
      foil.wake_shape.angle_panel.size() < total_nodes) {
    return false;
  }

  return foil.wake_shape.points.block(0, point_count, 2, wake_point_count).allFinite();
}

bool XFoil::hasPanelMap() const {
  const auto& top = boundaryLayerWorkflow.lattice.top;
  const auto& bottom = boundaryLayerWorkflow.lattice.bottom;
  const int point_count = foil.foil_shape.n;
  const int total_nodes = point_count + foil.wake_shape.n;

  auto isValidSide = [total_nodes](const BoundaryLayerLattice& lattice) {
    if (lattice.stationCount <= 1 ||
        lattice.stationToPanel.size() < lattice.stationCount ||
        lattice.panelInfluenceFactor.size() < lattice.stationCount) {
      return false;
    }
    for (int i = 0; i < lattice.stationCount - 1; ++i) {
      const int panel = lattice.stationToPanel[i];
      if (panel < 0 || panel >= total_nodes) {
        return false;
      }
    }
    return true;
  };

  if (!isValidSide(top) || !isValidSide(bottom)) {
    return false;
  }

  return top.trailingEdgeIndex >= 0 &&
         top.trailingEdgeIndex < top.stationCount &&
         bottom.trailingEdgeIndex >= 0 &&
         bottom.trailingEdgeIndex < bottom.stationCount;
}

bool XFoil::hasAirfoilInfluenceMatrix() const {
  const int point_count = foil.foil_shape.n;
  const int total_nodes = point_count + foil.wake_shape.n;
  if (aerodynamicCache.dij.rows() < total_nodes ||
      aerodynamicCache.dij.cols() < total_nodes) {
    return false;
  }
  const auto block = aerodynamicCache.dij.block(0, 0, point_count, point_count);
  return block.allFinite() && block.cwiseAbs().maxCoeff() > 0.0;
}

bool XFoil::hasWakeInfluenceMatrix() const {
  const int point_count = foil.foil_shape.n;
  const int wake_point_count = foil.wake_shape.n;
  const int total_nodes = point_count + wake_point_count;
  if (aerodynamicCache.dij.rows() < total_nodes ||
      aerodynamicCache.dij.cols() < total_nodes) {
    return false;
  }
  const auto block = aerodynamicCache.dij.block(point_count, point_count,
                                                wake_point_count, wake_point_count);
  return block.allFinite() && block.cwiseAbs().maxCoeff() > 0.0;
}

bool XFoil::hasConvergedSolution() const {
  if (!std::isfinite(avisc) || !std::isfinite(mvisc)) {
    return false;
  }
  const double alpha_tol = 1.0e-12;
  const double mach_tol = 1.0e-12;
  if (std::fabs(analysis_state_.alpha - avisc) > alpha_tol ||
      std::fabs(analysis_state_.currentMach - mvisc) > mach_tol) {
    return false;
  }

  const int total_nodes_with_wake = foil.foil_shape.n + foil.wake_shape.n;
  return qvis.size() >= total_nodes_with_wake &&
         qvis.head(total_nodes_with_wake).allFinite();
}

double XFoil::VAccel() {
  return vaccel_;
}

void XFoil::setVAccel(double accel) {
  vaccel_ = accel;
}
