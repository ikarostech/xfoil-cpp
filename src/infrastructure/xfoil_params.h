/****************************************************************************

        XFoil Parameters

        Copyright (C) 2008-2018 Andre Deperrois

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

#pragma once

// XFoil Direct Parameters - refer to XFoil documentation
#define IQX 302  /**< 300 = number of surface panel nodes + 6 */
#define IWX 50   /**< number of wake panel nodes */
#define IBX 604  /**< 600 number of buffer airfoil nodes = 2*IQX */
#define IZX 350  /**< 350 = number of panel nodes [airfoil + wake] */
#define IVX \
  302 /**< 300 = number of nodes along bl on one side of airfoil and wake. */
