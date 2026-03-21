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

#include "application/xfoil/XFoil.h"
#include "Eigen/Core"
#include <cmath>
#include <limits>
using namespace Eigen;

XFoil::CompressibilityParams XFoil::buildCompressibilityParams() const {
  const double current_mach = analysis_state_.currentMach;
  const double beta = std::sqrt(1.0 - current_mach * current_mach);
  const double beta_msq = -0.5 / beta;
  const double prandtlGlauertFactor =
      0.5 * current_mach * current_mach / (1.0 + beta);
  const double prandtlGlauertFactor_msq =
      0.5 / (1.0 + beta) - prandtlGlauertFactor / (1.0 + beta) * beta_msq;
  const double karmanTsienFactor =
      (current_mach / (1.0 + beta)) * (current_mach / (1.0 + beta));
  const double karmanTsienFactor_msq =
      1.0 / ((1.0 + beta) * (1.0 + beta)) -
      2.0 * karmanTsienFactor / (1.0 + beta) * beta_msq;
  return {beta,
          beta_msq,
          karmanTsienFactor,
          karmanTsienFactor_msq,
          prandtlGlauertFactor,
          prandtlGlauertFactor_msq};
}

XFoil::PressureCoefficientResult
XFoil::computePressureCoefficient(double tangential_velocity,
                                  double velocity_derivative,
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
  state_.viscous.convergedMach = 0.0;

  // initialize transition parameters until user changes them
  acrit = 9.0;
  boundaryLayer.setTransitionLocations(1.0, 1.0);

  //---- initialize freestream mach number to zero
  analysis_state_.machType = MachType::CONSTANT;
  analysis_state_.referenceMach = 0.0;

  //---- drop tolerance for bl system solver
  setVAccel(0.01);
  //---- default viscous parameters
  analysis_state_.reynoldsType = ReynoldsType::CONSTANT;
  analysis_state_.referenceRe = 0.0;
}

XFoil::~XFoil() = default;

bool XFoil::isBLInitialized() const {
  if (!hasPanelMap()) {
    return false;
  }
  const auto top = boundaryLayer.readSideModel(1);
  const auto bottom = boundaryLayer.readSideModel(2);
  if (!top.hasStations || !bottom.hasStations) {
    return false;
  }
  return top.hasFiniteThickness && bottom.hasFiniteThickness;
}

void XFoil::setBLInitialized(bool bInitialized) {
  if (bInitialized) {
    return;
  }
  boundaryLayer.zeroProfiles();
  invalidateConvergedSolution();
}

void XFoil::invalidateConvergedSolution() {
  state_.viscous.convergedAlpha = std::numeric_limits<double>::quiet_NaN();
  state_.viscous.convergedMach = std::numeric_limits<double>::quiet_NaN();
}

void XFoil::invalidateWakeGeometry() {
  const int point_count = foil.foil_shape.n;
  const int wake_point_count = foil.wake_shape.n;
  if (const int total_nodes = point_count + wake_point_count;
      wake_point_count > 0 &&
          state_.inviscid.cache.dij.rows() >= total_nodes &&
      state_.inviscid.cache.dij.cols() >= total_nodes) {
    state_.inviscid.cache.dij
        .block(point_count, point_count, wake_point_count, wake_point_count)
        .setConstant(std::numeric_limits<double>::quiet_NaN());
  }
  foil.wake_shape.points.resize(0, 0);
  foil.wake_shape.normal_vector.resize(0, 0);
  foil.wake_shape.spline_length.resize(0);
  foil.wake_shape.angle_panel.resize(0);
}

void XFoil::invalidatePanelMap() {
  boundaryLayer.clearPanelMap();
}

bool XFoil::hasPanelMap() const {
  const int point_count = foil.foil_shape.n;
  const int total_nodes = point_count + foil.wake_shape.n;
  return boundaryLayer.hasValidPanelMap(total_nodes);
}

bool XFoil::hasAirfoilInfluenceMatrix() const {
  const int point_count = foil.foil_shape.n;
  if (const int total_nodes = point_count + foil.wake_shape.n;
      state_.inviscid.cache.dij.rows() < total_nodes ||
      state_.inviscid.cache.dij.cols() < total_nodes) {
    return false;
  }
  const auto block =
      state_.inviscid.cache.dij.block(0, 0, point_count, point_count);
  return block.allFinite() && block.cwiseAbs().maxCoeff() > 0.0;
}

bool XFoil::hasWakeInfluenceMatrix() const {
  const int point_count = foil.foil_shape.n;
  const int wake_point_count = foil.wake_shape.n;
  if (const int total_nodes = point_count + wake_point_count;
      state_.inviscid.cache.dij.rows() < total_nodes ||
      state_.inviscid.cache.dij.cols() < total_nodes) {
    return false;
  }
  const auto block = state_.inviscid.cache.dij.block(
      point_count, point_count, wake_point_count, wake_point_count);
  return block.allFinite() && block.cwiseAbs().maxCoeff() > 0.0;
}

bool XFoil::hasConvergedSolution() const {
  if (!std::isfinite(state_.viscous.convergedAlpha) ||
      !std::isfinite(state_.viscous.convergedMach)) {
    return false;
  }
  const double alpha_tol = 1.0e-12;
  if (const double mach_tol = 1.0e-12;
      std::fabs(analysis_state_.alpha - state_.viscous.convergedAlpha) >
          alpha_tol ||
      std::fabs(analysis_state_.currentMach - state_.viscous.convergedMach) >
          mach_tol) {
    return false;
  }

  const int total_nodes_with_wake = foil.foil_shape.n + foil.wake_shape.n;
  return state_.viscous.qvis.size() >= total_nodes_with_wake &&
         state_.viscous.qvis.head(total_nodes_with_wake).allFinite();
}

double XFoil::VAccel() { return vaccel_; }

void XFoil::setVAccel(double accel) { vaccel_ = accel; }
