/**
 * BL Newton system assembly split from XFoil.cpp
 */

#include "XFoil.h"

using Eigen::RowVector;
using Eigen::Vector3d;

struct LogarithmicDifferences {
    double xlog;
    double ulog;
    double tlog;
    double hlog;
    double ddlog;
  };
LogarithmicDifferences getLogarithmicDifferences(int flowRegimeType, BoundaryLayerState boundaryLayerState) {
  LogarithmicDifferences logDiffs;
  if (flowRegimeType == 0) {
        //----- similarity logarithmic differences  (prescribed)
    logDiffs.xlog = 1.0;
    logDiffs.ulog = 1.0;
    logDiffs.tlog = 0.0;
    logDiffs.hlog = 0.0;
    logDiffs.ddlog = 0.0;
    logDiffs.ddlog = 0.0;
  } else {
    logDiffs.xlog = log(boundaryLayerState.station2.param.xz / boundaryLayerState.station1.param.xz);
    logDiffs.ulog = log(boundaryLayerState.station2.param.uz / boundaryLayerState.station1.param.uz);
    logDiffs.tlog = log(boundaryLayerState.station2.param.tz / boundaryLayerState.station1.param.tz);
    logDiffs.hlog = log(boundaryLayerState.station2.hsz.scalar / boundaryLayerState.station1.hsz.scalar);
    logDiffs.ddlog = 1.0;
  }
  return logDiffs;
}
XFoil::BlSystemCoeffs XFoil::bldif(int flowRegimeType) const {
  LogarithmicDifferences logDiffs = getLogarithmicDifferences(flowRegimeType, boundaryLayerState);

  BlSystemCoeffs coeffs;
  coeffs.clear();

  //---- set triggering constant for local upwinding
  //---- use less upwinding in the wake
  double hdcon;
  if (flowRegimeType == 3) {
    hdcon = 1.0 / blData2.hkz.scalar / blData2.hkz.scalar;
  } else {
    hdcon = 5.0 / MathUtil::pow(blData2.hkz.scalar, 2);
  }
  double hd_hk2 = -hdcon * 2.0 / blData2.hkz.scalar;

  //---- local upwinding is based on local change in  log(hk-1)
  //-    (mainly kicks in at transition)
  double hl = log(fabs((blData2.hkz.scalar - 1.0) / (blData1.hkz.scalar - 1.0)));
  double hl_hk1 = -1.0 / (blData1.hkz.scalar - 1.0);
  double hl_hk2 = 1.0 / (blData2.hkz.scalar - 1.0);

  double hlsq = std::min(hl * hl, 15.0);
  double ehh = exp(-hlsq * hdcon);

  double upw = 1.0 - 0.5 * ehh;
  double upw_hl = ehh * hl * hdcon;
  double upw_hd = ehh * hlsq / 2;

  double upw_hk1 = upw_hl * hl_hk1;
  double upw_hk2 = upw_hl * hl_hk2 + upw_hd * hd_hk2;

  Vector3d upw1 = upw_hk1 * blData1.hkz.pos_vector();
  Vector3d upw2 = upw_hk2 * blData2.hkz.pos_vector();
  double upw_ms = upw_hk1 * blData1.hkz.ms() + upw_hk2 * blData2.hkz.ms();

  if (flowRegimeType == 0) {
    //***** le point -->  set zero amplification factor
    coeffs.a2(0, 0) = 1.0;
    coeffs.d_re[0] = 0.0;
    coeffs.rhs[0] = -blData2.param.amplz;
  } else if (flowRegimeType == 1) {
    //----- build laminar amplification equation
    bldifLaminar(coeffs);
  } else {
    //----- build turbulent or wake shear lag equation
    bldifTurbulent(static_cast<FlowRegimeEnum>(flowRegimeType), upw, upw1, upw2,
                   upw_ms, logDiffs.ulog, coeffs);
  }

  //**** set up momentum equation
  bldifMomentum(logDiffs.xlog, logDiffs.ulog, logDiffs.tlog, logDiffs.ddlog, coeffs);

  //**** set up shape parameter equation
  bldifShape(upw, logDiffs.xlog, logDiffs.ulog, logDiffs.hlog, logDiffs.ddlog, upw1, upw2, upw_ms, coeffs);

  return coeffs;
}

void XFoil::bldifLaminar(BlSystemCoeffs& coeffs) const {
  AxResult ax_result =
      axset(blData1.hkz.scalar, blData1.param.tz, blData1.rtz.scalar,
            blData1.param.amplz, blData2.hkz.scalar, blData2.param.tz,
            blData2.rtz.scalar, blData2.param.amplz, amcrit);

  double rezc = blData2.param.amplz - blData1.param.amplz -
                ax_result.ax * (blData2.param.xz - blData1.param.xz);
  double z_ax = -(blData2.param.xz - blData1.param.xz);

  coeffs.a1(0, 0) = z_ax * ax_result.ax_a1 - 1.0;
  coeffs.a1(0, 1) = z_ax * (ax_result.ax_hk1 * blData1.hkz.t() + ax_result.ax_t1 +
                      ax_result.ax_rt1 * blData1.rtz.t());
  coeffs.a1(0, 2) = z_ax * (ax_result.ax_hk1 * blData1.hkz.d());
  coeffs.a1(0, 3) = z_ax * (ax_result.ax_hk1 * blData1.hkz.u() +
                      ax_result.ax_rt1 * blData1.rtz.u());
  coeffs.a1(0, 4) = ax_result.ax;
  coeffs.a2(0, 0) = z_ax * ax_result.ax_a2 + 1.0;
  coeffs.a2(0, 1) = z_ax * (ax_result.ax_hk2 * blData2.hkz.t() + ax_result.ax_t2 +
                      ax_result.ax_rt2 * blData2.rtz.t());
  coeffs.a2(0, 2) = z_ax * (ax_result.ax_hk2 * blData2.hkz.d());
  coeffs.a2(0, 3) = z_ax * (ax_result.ax_hk2 * blData2.hkz.u() +
                      ax_result.ax_rt2 * blData2.rtz.u());
  coeffs.a2(0, 4) = -ax_result.ax;
  coeffs.d_msq[0] = z_ax * (ax_result.ax_hk1 * blData1.hkz.ms() +
                   ax_result.ax_rt1 * blData1.rtz.ms() +
                   ax_result.ax_hk2 * blData2.hkz.ms() +
                   ax_result.ax_rt2 * blData2.rtz.ms());
  coeffs.d_re[0] = z_ax * (ax_result.ax_rt1 * blData1.rtz.re() +
                   ax_result.ax_rt2 * blData2.rtz.re());
  coeffs.d_xi[0] = 0.0;
  coeffs.rhs[0] = -rezc;
}

void XFoil::bldifTurbulent(FlowRegimeEnum flowRegimeType, double upw,
                           const Vector3d &upw1, const Vector3d &upw2,
                           double upw_ms, double ulog, BlSystemCoeffs& coeffs) const {
  double sa = (1.0 - upw) * blData1.param.sz + upw * blData2.param.sz;
  double cqa = (1.0 - upw) * blData1.cqz.scalar + upw * blData2.cqz.scalar;
  double cfa = (1.0 - upw) * blData1.cfz.scalar + upw * blData2.cfz.scalar;
  double hka = (1.0 - upw) * blData1.hkz.scalar + upw * blData2.hkz.scalar;

  double usa = 0.5 * (blData1.usz.scalar + blData2.usz.scalar);
  double rta = 0.5 * (blData1.rtz.scalar + blData2.rtz.scalar);
  double dea = 0.5 * (blData1.dez.scalar + blData2.dez.scalar);
  double da = 0.5 * (blData1.param.dz + blData2.param.dz);

  double ald = (flowRegimeType == FlowRegimeEnum::Wake) ? dlcon : 1.0;

  double gcc, hkc, hkc_hka;
  if (flowRegimeType == FlowRegimeEnum::Turbulent) {
    gcc = gccon;
    hkc = hka - 1.0 - gcc / rta;
    hkc_hka = 1.0;
    if (hkc < 0.01) {
      hkc = 0.01;
      hkc_hka = 0.0;
    }
  } else {
    hkc = hka - 1.0;
    hkc_hka = 1.0;
  }

  double hr = hkc / (gacon * ald * hka);
  double hr_hka = hkc_hka / (gacon * ald * hka) - hr / hka;

  double uq = (0.5 * cfa - hr * hr) / (gbcon * da);
  double uq_hka = -2.0 * hr * hr_hka / (gbcon * da);
  double uq_cfa = 0.5 / (gbcon * da);
  double uq_da = -uq / da;

  double scc = sccon * 1.333 / (1.0 + usa);
  double scc_usa = -scc / (1.0 + usa);

  double slog = log(blData2.param.sz / blData1.param.sz);
  double dxi = blData2.param.xz - blData1.param.xz;

  double rezc = scc * (cqa - sa * ald) * dxi - dea * 2.0 * slog +
                dea * 2.0 * (uq * dxi - ulog);

  double z_cfa = dea * 2.0 * uq_cfa * dxi;
  double z_hka = dea * 2.0 * uq_hka * dxi;
  double z_da = dea * 2.0 * uq_da * dxi;
  double z_sl = -dea * 2.0;
  double z_ul = -dea * 2.0;
  double z_dxi = scc * (cqa - sa * ald) + dea * 2.0 * uq;
  double z_usa = scc_usa * (cqa - sa * ald) * dxi;
  double z_cqa = scc * dxi;
  double z_sa = -scc * dxi * ald;
  double z_dea = 2.0 * (uq * dxi - ulog - slog);

  double z_upw = z_cqa * (blData2.cqz.scalar - blData1.cqz.scalar) +
                 z_sa * (blData2.param.sz - blData1.param.sz) +
                 z_cfa * (blData2.cfz.scalar - blData1.cfz.scalar) +
                 z_hka * (blData2.hkz.scalar - blData1.hkz.scalar);
  double z_de = 0.5 * z_dea;
  double z_us = 0.5 * z_usa;
  double z_d = 0.5 * z_da;
  double z_u1 = -z_ul / blData1.param.uz;
  double z_u2 = z_ul / blData2.param.uz;
  double z_x1 = -z_dxi;
  double z_x2 = z_dxi;
  double z_s1 = (1.0 - upw) * z_sa - z_sl / blData1.param.sz;
  double z_s2 = upw * z_sa + z_sl / blData2.param.sz;
  double z_cq1 = (1.0 - upw) * z_cqa;
  double z_cq2 = upw * z_cqa;
  double z_cf1 = (1.0 - upw) * z_cfa;
  double z_cf2 = upw * z_cfa;
  double z_hk1 = (1.0 - upw) * z_hka;
  double z_hk2 = upw * z_hka;

  RowVector<double, 5> vs1_row;
  RowVector<double, 5> vs2_row;

  vs1_row << z_s1, 0.0, z_d, z_u1, z_x1;
  vs2_row << z_s2, 0.0, z_d, z_u2, z_x2;

  Vector3d vs1_vec =
      z_upw * upw1 + z_de * blData1.dez.pos_vector() +
      z_us * blData1.usz.pos_vector() + z_cq1 * blData1.cqz.pos_vector() +
      z_cf1 * blData1.cfz.pos_vector() + z_hk1 * blData1.hkz.pos_vector();

  Vector3d vs2_vec =
      z_upw * upw2 + z_de * blData2.dez.pos_vector() +
      z_us * blData2.usz.pos_vector() + z_cq2 * blData2.cqz.pos_vector() +
      z_cf2 * blData2.cfz.pos_vector() + z_hk2 * blData2.hkz.pos_vector();

  vs1_row.segment<3>(1) += vs1_vec;
  vs2_row.segment<3>(1) += vs2_vec;

  coeffs.a1.row(0) = vs1_row;
  coeffs.a2.row(0) = vs2_row;

  coeffs.d_msq[0] = z_upw * upw_ms + z_de * blData1.dez.ms() + z_us * blData1.usz.ms() +
           z_de * blData2.dez.ms() + z_us * blData2.usz.ms() +
           z_cq1 * blData1.cqz.ms() + z_cf1 * blData1.cfz.ms() +
           z_hk1 * blData1.hkz.ms() + z_cq2 * blData2.cqz.ms() +
           z_cf2 * blData2.cfz.ms() + z_hk2 * blData2.hkz.ms();

  coeffs.d_re[0] = z_cq1 * blData1.cqz.re() + z_cf1 * blData1.cfz.re() +
           z_cq2 * blData2.cqz.re() + z_cf2 * blData2.cfz.re();
  coeffs.d_xi[0] = 0.0;
  coeffs.rhs[0] = -rezc;
}

void XFoil::bldifMomentum(double xlog, double ulog, double tlog, double ddlog,
                          BlSystemCoeffs& coeffs) const {
  double ha = 0.5 * (blData1.param.hz + blData2.param.hz);
  double ma = 0.5 * (blData1.param.mz + blData2.param.mz);
  double xa = 0.5 * (blData1.param.xz + blData2.param.xz);
  double ta = 0.5 * (blData1.param.tz + blData2.param.tz);
  double hwa = 0.5 * (blData1.param.dwz / blData1.param.tz +
                      blData2.param.dwz / blData2.param.tz);

  double cfx =
      0.50 * cfm * xa / ta +
      0.25 * (blData1.cfz.scalar * blData1.param.xz / blData1.param.tz +
              blData2.cfz.scalar * blData2.param.xz / blData2.param.tz);
  double cfx_xa = 0.50 * cfm / ta;
  double cfx_ta = -.50 * cfm * xa / ta / ta;

  double cfx_x1 = 0.25 * blData1.cfz.scalar / blData1.param.tz + cfx_xa * 0.5;
  double cfx_x2 = 0.25 * blData2.cfz.scalar / blData2.param.tz + cfx_xa * 0.5;
  double cfx_t1 = -.25 * blData1.cfz.scalar * blData1.param.xz /
                      blData1.param.tz / blData1.param.tz +
                  cfx_ta * 0.5;
  double cfx_t2 = -.25 * blData2.cfz.scalar * blData2.param.xz /
                      blData2.param.tz / blData2.param.tz +
                  cfx_ta * 0.5;
  double cfx_cf1 = 0.25 * blData1.param.xz / blData1.param.tz;
  double cfx_cf2 = 0.25 * blData2.param.xz / blData2.param.tz;
  double cfx_cfm = 0.50 * xa / ta;

  double btmp = ha + 2.0 - ma + hwa;

  double rezt = tlog + btmp * ulog - xlog * 0.5 * cfx;
  double z_cfx = -xlog * 0.5;
  double z_ha = ulog;
  double z_hwa = ulog;
  double z_ma = -ulog;
  double z_xl = -ddlog * 0.5 * cfx;
  double z_ul = ddlog * btmp;
  double z_tl = ddlog;

  double z_cfm = z_cfx * cfx_cfm;
  double z_cf1 = z_cfx * cfx_cf1;
  double z_cf2 = z_cfx * cfx_cf2;

  double z_t1 =
      -z_tl / blData1.param.tz + z_cfx * cfx_t1 +
      z_hwa * 0.5 * (-blData1.param.dwz / blData1.param.tz / blData1.param.tz);
  double z_t2 =
      z_tl / blData2.param.tz + z_cfx * cfx_t2 +
      z_hwa * 0.5 * (-blData2.param.dwz / blData2.param.tz / blData2.param.tz);
  double z_x1 = -z_xl / blData1.param.xz + z_cfx * cfx_x1;
  double z_x2 = z_xl / blData2.param.xz + z_cfx * cfx_x2;
  double z_u1 = -z_ul / blData1.param.uz;
  double z_u2 = z_ul / blData2.param.uz;

  // Row k=1: build with vector operations (t,d,u columns as a segment)
  {
    RowVector<double, 5> row1_a1 = RowVector<double, 5>::Zero();
    RowVector<double, 5> row1_a2 = RowVector<double, 5>::Zero();

    Vector3d hterm1(blData1.param.hz_tz, blData1.param.hz_dz, 0.0);
    Vector3d hterm2(blData2.param.hz_tz, blData2.param.hz_dz, 0.0);
    Vector3d cfm1v(cfm_t1, cfm_d1, cfm_u1);
    Vector3d cfm2v(cfm_t2, cfm_d2, cfm_u2);
    Vector3d cfz1v = blData1.cfz.pos_vector();
    Vector3d cfz2v = blData2.cfz.pos_vector();
    Vector3d mz1(0.0, 0.0, 0.5 * z_ma * blData1.param.mz_uz);
    Vector3d mz2(0.0, 0.0, 0.5 * z_ma * blData2.param.mz_uz);

    Vector3d seg1 = 0.5 * z_ha * hterm1 + z_cfm * cfm1v + z_cf1 * cfz1v + mz1 + Vector3d(z_t1, 0.0, z_u1);
    Vector3d seg2 = 0.5 * z_ha * hterm2 + z_cfm * cfm2v + z_cf2 * cfz2v + mz2 + Vector3d(z_t2, 0.0, z_u2);

    row1_a1.segment<3>(1) = seg1;
    row1_a2.segment<3>(1) = seg2;
    row1_a1(4) = z_x1;
    row1_a2(4) = z_x2;

    coeffs.a1.row(1) = row1_a1;
    coeffs.a2.row(1) = row1_a2;
  }

  coeffs.d_msq[1] = 0.5 * z_ma * blData1.param.mz_ms + z_cfm * cfm_ms +
           z_cf1 * blData1.cfz.ms() + 0.5 * z_ma * blData2.param.mz_ms +
           z_cf2 * blData2.cfz.ms();
  coeffs.d_re[1] = z_cfm * cfm_re + z_cf1 * blData1.cfz.re() + z_cf2 * blData2.cfz.re();
  coeffs.d_xi[1] = 0.0;
  coeffs.rhs[1] = -rezt;
}

void XFoil::bldifShape(double upw, double xlog, double ulog, double hlog,
                       double ddlog, const Vector3d &upw1, const Vector3d &upw2,
                       double upw_ms, BlSystemCoeffs& coeffs) const {
  double xot1 = blData1.param.xz / blData1.param.tz;
  double xot2 = blData2.param.xz / blData2.param.tz;

  double ha = 0.5 * (blData1.param.hz + blData2.param.hz);
  double hsa = 0.5 * (blData1.hsz.scalar + blData2.hsz.scalar);
  double hca = 0.5 * (blData1.hcz.scalar + blData2.hcz.scalar);
  double hwa = 0.5 * (blData1.param.dwz / blData1.param.tz +
                      blData2.param.dwz / blData2.param.tz);

  double dix =
      (1.0 - upw) * blData1.diz.scalar * xot1 + upw * blData2.diz.scalar * xot2;
  double cfx =
      (1.0 - upw) * blData1.cfz.scalar * xot1 + upw * blData2.cfz.scalar * xot2;
  double dix_upw = blData2.diz.scalar * xot2 - blData1.diz.scalar * xot1;
  double cfx_upw = blData2.cfz.scalar * xot2 - blData1.cfz.scalar * xot1;

  double btmp = 2.0 * hca / hsa + 1.0 - ha - hwa;

  double rezh = hlog + btmp * ulog + xlog * (0.5 * cfx - dix);
  double z_cfx = xlog * 0.5;
  double z_dix = -xlog;
  double z_hca = 2.0 * ulog / hsa;
  double z_ha = -ulog;
  double z_hwa = -ulog;
  double z_xl = ddlog * (0.5 * cfx - dix);
  double z_ul = ddlog * btmp;
  double z_hl = ddlog;

  double z_upw = z_cfx * cfx_upw + z_dix * dix_upw;

  double z_hs1 = -hca * ulog / hsa / hsa - z_hl / blData1.hsz.scalar;
  double z_hs2 = -hca * ulog / hsa / hsa + z_hl / blData2.hsz.scalar;

  double z_cf1 = (1.0 - upw) * z_cfx * xot1;
  double z_cf2 = upw * z_cfx * xot2;
  double z_di1 = (1.0 - upw) * z_dix * xot1;
  double z_di2 = upw * z_dix * xot2;

  double z_t1 = (1.0 - upw) *
                (z_cfx * blData1.cfz.scalar + z_dix * blData1.diz.scalar) *
                (-xot1 / blData1.param.tz);
  double z_t2 = upw *
                (z_cfx * blData2.cfz.scalar + z_dix * blData2.diz.scalar) *
                (-xot2 / blData2.param.tz);
  double z_x1 = (1.0 - upw) *
                    (z_cfx * blData1.cfz.scalar + z_dix * blData1.diz.scalar) /
                    blData1.param.tz -
                z_xl / blData1.param.xz;
  double z_x2 = upw *
                    (z_cfx * blData2.cfz.scalar + z_dix * blData2.diz.scalar) /
                    blData2.param.tz +
                z_xl / blData2.param.xz;
  double z_u1 = -z_ul / blData1.param.uz;
  double z_u2 = z_ul / blData2.param.uz;

  z_t1 +=
      z_hwa * 0.5 * (-blData1.param.dwz / blData1.param.tz / blData1.param.tz);
  z_t2 +=
      z_hwa * 0.5 * (-blData2.param.dwz / blData2.param.tz / blData2.param.tz);

  // Row k=2: vector assembly
  {
    RowVector<double, 5> row2_a1 = RowVector<double, 5>::Zero();
    RowVector<double, 5> row2_a2 = RowVector<double, 5>::Zero();

    Vector3d base1 = z_hs1 * blData1.hsz.pos_vector() +
                     z_cf1 * blData1.cfz.pos_vector() +
                     z_di1 * blData1.diz.pos_vector() +
                     Vector3d(z_t1, 0.0, z_u1);
    Vector3d base2 = z_hs2 * blData2.hsz.pos_vector() +
                     z_cf2 * blData2.cfz.pos_vector() +
                     z_di2 * blData2.diz.pos_vector() +
                     Vector3d(z_t2, 0.0, z_u2);

    row2_a1(0) = z_di1 * blData1.diz.s();
    row2_a2(0) = z_di2 * blData2.diz.s();
    row2_a1.segment<3>(1) = base1;
    row2_a2.segment<3>(1) = base2;
    row2_a1(4) = z_x1;
    row2_a2(4) = z_x2;

    coeffs.a1.row(2) = row2_a1;
    coeffs.a2.row(2) = row2_a2;
  }
  coeffs.d_msq[2] = z_hs1 * blData1.hsz.ms() + z_cf1 * blData1.cfz.ms() +
           z_di1 * blData1.diz.ms() + z_hs2 * blData2.hsz.ms() +
           z_cf2 * blData2.cfz.ms() + z_di2 * blData2.diz.ms();
  coeffs.d_re[2] = z_hs1 * blData1.hsz.re() + z_cf1 * blData1.cfz.re() +
           z_di1 * blData1.diz.re() + z_hs2 * blData2.hsz.re() +
           z_cf2 * blData2.cfz.re() + z_di2 * blData2.diz.re();

  // Column t,d,u increments as a vector
  {
    Vector3d inc1 = 0.5 * (z_hca * blData1.hcz.pos_vector() +
                           z_ha * Vector3d(blData1.param.hz_tz, blData1.param.hz_dz, 0.0)) +
                    z_upw * upw1;
    Vector3d inc2 = 0.5 * (z_hca * blData2.hcz.pos_vector() +
                           z_ha * Vector3d(blData2.param.hz_tz, blData2.param.hz_dz, 0.0)) +
                    z_upw * upw2;
    coeffs.a1.row(2).segment<3>(1) += inc1;
    coeffs.a2.row(2).segment<3>(1) += inc2;
  }

  coeffs.d_msq[2] = 0.5 * (z_hca * blData1.hcz.ms()) + z_upw * upw_ms +
           0.5 * (z_hca * blData2.hcz.ms());

  coeffs.d_xi[2] = 0.0;
  coeffs.rhs[2] = -rezh;
}
