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
using namespace Eigen;

// determinant
double cross2(const Eigen::Vector2d &a, const Eigen::Vector2d &b) {
  return a[0] * b[1] - a[1] * b[0];
}

bool XFoil::s_bCancel = false;
double XFoil::vaccel = 0.01;

XFoil::XFoil() {
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
  vaccel = 0.01;
  //---- default viscous parameters
  reynolds_type = ReynoldsType::CONSTANT;
  reinf1 = 0.0;
}

XFoil::~XFoil() {}
