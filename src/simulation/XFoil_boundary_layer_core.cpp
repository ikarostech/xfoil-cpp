#include "XFoil.h"
#include "domain/coefficient/skin_friction.hpp"
#include "domain/coefficient/dissipation.hpp"
#include <algorithm>
#include <cmath>
using namespace Eigen;

bool XFoil::blkin(BoundaryLayerState& state) {
  //----------------------------------------------------------
  //     calculates turbulence-independent secondary "2"
  //     variables from the primary "2" variables.
  //----------------------------------------------------------
  blData& current = state.current();
  //---- set edge mach number ** 2
  current.param.mz =
      current.param.uz * current.param.uz * hstinv /
      (gm1bl * (1.0 - 0.5 * current.param.uz * current.param.uz * hstinv));
  double tr2 = 1.0 + 0.5 * gm1bl * current.param.mz;
  current.param.mz_uz = 2.0 * current.param.mz * tr2 / current.param.uz;
  current.param.mz_ms =
      current.param.uz * current.param.uz * tr2 /
      (gm1bl * (1.0 - 0.5 * current.param.uz * current.param.uz * hstinv)) *
      hstinv_ms;

  //---- set edge density (isentropic relation)
  current.param.rz = rstbl * pow(tr2, (-1.0 / gm1bl));
  current.param.rz_uz = -current.param.rz / tr2 * 0.5 * current.param.mz_uz;
  current.param.rz_ms = -current.param.rz / tr2 * 0.5 * current.param.mz_ms +
                        rstbl_ms * pow(tr2, (-1.0 / gm1bl));

  //---- set shape parameter
  current.param.hz = current.param.dz / current.param.tz;
  current.param.hz_dz = 1.0 / current.param.tz;
  current.param.hz_tz = -current.param.hz / current.param.tz;

  //---- set edge static/stagnation enthalpy
  double herat = 1.0 - 0.5 * current.param.uz * current.param.uz * hstinv;
  double he_u2 = -current.param.uz * hstinv;
  double he_ms = -0.5 * current.param.uz * current.param.uz * hstinv_ms;
  //---- set molecular viscosity
  double v2_he = (1.5 / herat - 1.0 / (herat + hvrat));

  //---- set kinematic shape parameter
  boundary_layer::KineticShapeParameterResult hkin_result =
      boundary_layer::hkin(current.param.hz, current.param.mz);
  current.hkz.scalar = hkin_result.hk;

  current.hkz.u() = hkin_result.hk_msq * current.param.mz_uz;
  current.hkz.t() = hkin_result.hk_h * current.param.hz_tz;
  current.hkz.d() = hkin_result.hk_h * current.param.hz_dz;
  current.hkz.ms() = hkin_result.hk_msq * current.param.mz_ms;

  //---- set momentum thickness reynolds number
  current.rtz.scalar =
      current.param.rz * current.param.uz * current.param.tz /
      (sqrt(herat * herat * herat) * (1.0 + hvrat) / (herat + hvrat) / reybl);
  current.rtz.u() = current.rtz.scalar *
                    (1.0 / current.param.uz +
                     current.param.rz_uz / current.param.rz - v2_he * he_u2);
  current.rtz.t() = current.rtz.scalar / current.param.tz;
  current.rtz.ms() =
      current.rtz.scalar * (current.param.rz_ms / current.param.rz +
                            (1 / reybl * reybl_ms - v2_he * he_ms));
  current.rtz.re() = current.rtz.scalar * (reybl_re / reybl);

  return true;
}


SkinFrictionCoefficients XFoil::blmid(BoundaryLayerState& state,
                                             FlowRegimeEnum flowRegimeType) {

  blData& previous = state.previous();
  blData& current = state.current();

  //---- set similarity variables if not defined
  if (simi) {
    previous.hkz = current.hkz;
    previous.rtz = current.rtz;
    previous.param.mz = current.param.mz;
    previous.param.mz_uz = current.param.mz_uz;
    previous.param.mz_ms = current.param.mz_ms;
  }

  //---- define stuff for midpoint cf
  double hka = 0.5 * (previous.hkz.scalar + current.hkz.scalar);
  double rta = 0.5 * (previous.rtz.scalar + current.rtz.scalar);
  double ma = 0.5 * (previous.param.mz + current.param.mz);

  //---- compute midpoint skin friction coefficient
  skin_friction::C_f cf_res = skin_friction::getSkinFriction(hka, rta, ma, flowRegimeType);

  SkinFrictionCoefficients coeffs;
  coeffs.cfm = cf_res.cf;
  double cfm_hka = cf_res.hk;
  double cfm_rta = cf_res.rt;
  double cfm_ma = cf_res.msq;

  coeffs.cfm_u1 = 0.5 * (cfm_hka * previous.hkz.u() + cfm_ma * previous.param.mz_uz +
                         cfm_rta * previous.rtz.u());
  coeffs.cfm_t1 = 0.5 * (cfm_hka * previous.hkz.t() + cfm_rta * previous.rtz.t());
  coeffs.cfm_d1 = 0.5 * (cfm_hka * previous.hkz.d());

  coeffs.cfm_u2 = 0.5 * (cfm_hka * current.hkz.u() + cfm_ma * current.param.mz_uz +
                         cfm_rta * current.rtz.u());
  coeffs.cfm_t2 = 0.5 * (cfm_hka * current.hkz.t() + cfm_rta * current.rtz.t());
  coeffs.cfm_d2 = 0.5 * (cfm_hka * current.hkz.d());

  coeffs.cfm_ms =
      0.5 * (cfm_hka * previous.hkz.ms() + cfm_ma * previous.param.mz_ms +
             cfm_rta * previous.rtz.ms() + cfm_hka * current.hkz.ms() +
             cfm_ma * current.param.mz_ms + cfm_rta * current.rtz.ms());
  coeffs.cfm_re = 0.5 * (cfm_rta * previous.rtz.re() + cfm_rta * current.rtz.re());

  return coeffs;
}

SkinFrictionCoefficients XFoil::blmid(FlowRegimeEnum flowRegimeType) {
  return blmid(boundaryLayerWorkflow.state, flowRegimeType);
}


/** ----------------------------------------------------------
 *     set bl primary "2" variables from parameter list
 *  ---------------------------------------------------------- */
blData XFoil::blprv(blData data, double xsi, double ami, double cti,
                    double thi, double dsi, double dswaki, double uei) const {
  data.param.xz = xsi;
  data.param.amplz = ami;
  data.param.sz = cti;
  data.param.tz = thi;
  data.param.dz = dsi - dswaki;
  data.param.dwz = dswaki;

  data.param.uz =
      uei * (1.0 - tkbl) / (1.0 - tkbl * (uei / qinfbl) * (uei / qinfbl));
  data.param.uz_uei =
      (1.0 + tkbl * (2.0 * data.param.uz * uei / qinfbl / qinfbl - 1.0)) /
      (1.0 - tkbl * (uei / qinfbl) * (uei / qinfbl));
  data.param.uz_ms =
      (data.param.uz * (uei / qinfbl) * (uei / qinfbl) - uei) * tkbl_ms /
      (1.0 - tkbl * (uei / qinfbl) * (uei / qinfbl));
  return data;
}


/** ------------------------------------------------------------------
 *
 *      sets up the bl newton system governing the current interval:
 *
 *      |       ||da1|     |       ||da2|       |     |
 *      |  vs1  ||dt1|  +  |  vs2  ||dt2|   =   |vsrez|
 *      |       ||dd1|     |       ||dd2|       |     |
 *               |du1|              |du2|
 *               |dx1|              |dx2|
 *
 *         3x5    5x1         3x5    5x1          3x1
 *
 *      the system as shown corresponds to a laminar station
 *      if tran, then  ds2  replaces  da2
 *      if turb, then  ds1, ds2  replace  da1, da2
 *
 * ------------------------------------------------------------------ */
