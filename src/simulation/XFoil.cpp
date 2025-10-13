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
#include "Eigen/Core"
#include <cmath>
using namespace Eigen;

void ClearAerodynamicsState(const XFoil& xfoil);
void ClearInitState(const XFoil& xfoil);

namespace {
bool g_cancel_flag = false;
double g_vaccel = 0.01;
}  // namespace

XFoil::CompressibilityParams XFoil::buildCompressibilityParams() const {
  const double beta = std::sqrt(1.0 - minf * minf);
  const double beta_msq = -0.5 / beta;
  const double bfac = 0.5 * minf * minf / (1.0 + beta);
  const double bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;
  return {beta, beta_msq, bfac, bfac_msq};
}

XFoil::PressureCoefficientResult XFoil::computePressureCoefficient(
    double tangential_velocity, double velocity_derivative,
    const CompressibilityParams &params) const {
  const double velocity_ratio = tangential_velocity / qinf;
  const double cginc = 1.0 - velocity_ratio * velocity_ratio;
  const double denom = params.beta + params.bfac * cginc;
  const double pressure_coefficient = cginc / denom;
  const double cp_msq =
      -pressure_coefficient / denom *
      (params.beta_msq + params.bfac_msq * cginc);

  double cp_velocity_derivative = 0.0;
  if (velocity_derivative != 0.0) {
    const double cpi = -2.0 * tangential_velocity / (qinf * qinf);
    const double cpc_cpi = (1.0 - params.bfac * pressure_coefficient) / denom;
    cp_velocity_derivative = cpc_cpi * cpi * velocity_derivative;
  }

  return {pressure_coefficient, cp_msq, cp_velocity_derivative};
}

Matrix2d XFoil::buildBodyToFreestreamRotation() const {
  const double ca = std::cos(alfa);
  const double sa = std::sin(alfa);
  Matrix2d rotation;
  rotation << ca, sa, -sa, ca;
  return rotation;
}

XFoil::XFoil()
    : blData1(boundaryLayerState.station1),
      blData2(boundaryLayerState.station2),
      ipan(boundaryLayerLattice.stationToPanel),
      isys(boundaryLayerLattice.stationToSystem),
      iblte(boundaryLayerLattice.trailingEdgeIndex),
      nbl(boundaryLayerLattice.stationCount),
      xstrip(boundaryLayerLattice.transitionLocation),
      thet(boundaryLayerLattice.thet),
      ctau(boundaryLayerLattice.ctau),
      dstr(boundaryLayerLattice.dstr),
      uedg(boundaryLayerLattice.uedg),
      ctq(boundaryLayerLattice.ctq),
      xssi(boundaryLayerLattice.xssi),
      uinv(boundaryLayerLattice.uinv),
      uinv_a(boundaryLayerLattice.uinv_a),
      mass(boundaryLayerLattice.mass),
      vti(boundaryLayerLattice.vti),
      itran(boundaryLayerLattice.transitionIndex) {
  m_pOutStream = nullptr;

  // fortran seems to initializes variables to 0
  mvisc = 0.0;

  // initialize transition parameters until user changes them
  acrit = 9.0;
  xstrip.top = 1.0;
  xstrip.bottom = 1.0;

  //---- initialize freestream mach number to zero
  mach_type = MachType::CONSTANT;
  minf1 = 0.0;

  //---- drop tolerance for bl system solver
  setVAccel(0.01);
  //---- default viscous parameters
  reynolds_type = ReynoldsType::CONSTANT;
  reinf1 = 0.0;
}

XFoil::~XFoil() {
  ClearAerodynamicsState(*this);
  ClearInitState(*this);
}

double XFoil::VAccel() {
  return g_vaccel;
}

void XFoil::setVAccel(double accel) {
  g_vaccel = accel;
}

bool XFoil::isCancelled() {
  return g_cancel_flag;
}

void XFoil::setCancel(bool bCancel) {
  g_cancel_flag = bCancel;
}
