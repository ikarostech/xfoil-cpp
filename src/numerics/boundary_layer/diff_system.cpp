#include "numerics/boundary_layer/diff_system.hpp"

#include <algorithm>
#include <cmath>

#include <Eigen/Core>

#include "numerics/boundary_layer_util.hpp"
#include "numerics/math_util.hpp"

using Eigen::RowVector;
using Eigen::Vector3d;
using std::exp;
using std::fabs;
using std::log;

constexpr double kSccon = 5.6;
constexpr double kGacon = 6.70;
constexpr double kGbcon = 0.75;
constexpr double kGccon = 18.0;
constexpr double kDlcon = 0.9;

struct LogarithmicDifferences {
    double xlog;
    double ulog;
    double tlog;
    double hlog;
    double ddlog;
};

static LogarithmicDifferences getLogarithmicDifferences(FlowRegimeEnum flowRegimeType,
                                                        BoundaryLayerStationWindow BoundaryLayerStationWindow) {
    LogarithmicDifferences logDiffs;
    if (flowRegimeType == FlowRegimeEnum::Similarity) {
        logDiffs.xlog  = 1.0;
        logDiffs.ulog  = 1.0;
        logDiffs.tlog  = 0.0;
        logDiffs.hlog  = 0.0;
        logDiffs.ddlog = 0.0;
    } else {
        logDiffs.xlog =
            log(BoundaryLayerStationWindow.station2.param.xz / BoundaryLayerStationWindow.station1.param.xz);
        logDiffs.ulog =
            log(BoundaryLayerStationWindow.station2.param.uz / BoundaryLayerStationWindow.station1.param.uz);
        logDiffs.tlog =
            log(BoundaryLayerStationWindow.station2.param.tz / BoundaryLayerStationWindow.station1.param.tz);
        logDiffs.hlog =
            log(BoundaryLayerStationWindow.station2.hsz.scalar / BoundaryLayerStationWindow.station1.hsz.scalar);
        logDiffs.ddlog = 1.0;
    }
    return logDiffs;
}

BlSystemCoeffs BlDiffSolver::solve(FlowRegimeEnum flowRegimeType, BoundaryLayerStationWindow BoundaryLayerStationWindow,
                                   const SkinFrictionCoefficients &skinFriction, double amcrit) {
    const LogarithmicDifferences logDiffs = getLogarithmicDifferences(flowRegimeType, BoundaryLayerStationWindow);

    blData &station1 = BoundaryLayerStationWindow.station1;
    blData &station2 = BoundaryLayerStationWindow.station2;

    BlSystemCoeffs coeffs;
    coeffs.clear();

    double hdcon;
    if (flowRegimeType == FlowRegimeEnum::Wake) {
        hdcon = 1.0 / station2.hkz.scalar / station2.hkz.scalar;
    } else {
        hdcon = 5.0 / MathUtil::pow(station2.hkz.scalar, 2);
    }
    const double hd_hk2 = -hdcon * 2.0 / station2.hkz.scalar;

    const double hl     = log(fabs((station2.hkz.scalar - 1.0) / (station1.hkz.scalar - 1.0)));
    const double hl_hk1 = -1.0 / (station1.hkz.scalar - 1.0);
    const double hl_hk2 = 1.0 / (station2.hkz.scalar - 1.0);

    const double hlsq = std::min(hl * hl, 15.0);
    const double ehh  = exp(-hlsq * hdcon);

    const double upw    = 1.0 - 0.5 * ehh;
    const double upw_hl = ehh * hl * hdcon;
    const double upw_hd = ehh * hlsq / 2;

    const double upw_hk1 = upw_hl * hl_hk1;
    const double upw_hk2 = upw_hl * hl_hk2 + upw_hd * hd_hk2;

    const Vector3d upw1 = upw_hk1 * station1.hkz.pos_vector();
    const Vector3d upw2 = upw_hk2 * station2.hkz.pos_vector();
    const double upw_ms = upw_hk1 * station1.hkz.ms() + upw_hk2 * station2.hkz.ms();

    if (flowRegimeType == FlowRegimeEnum::Similarity) {
        coeffs.a2(0, 0) = 1.0;
        coeffs.d_re[0]  = 0.0;
        coeffs.rhs[0]   = -station2.param.amplz;
    } else if (flowRegimeType == FlowRegimeEnum::Laminar) {
        bldifLaminar(BoundaryLayerStationWindow, amcrit, coeffs);
    } else {
        bldifTurbulent(BoundaryLayerStationWindow, flowRegimeType, upw, upw1, upw2, upw_ms, logDiffs.ulog, coeffs);
    }

    bldifMomentum(BoundaryLayerStationWindow, logDiffs.xlog, logDiffs.ulog, logDiffs.tlog, logDiffs.ddlog, skinFriction,
                  coeffs);
    bldifShape(BoundaryLayerStationWindow, upw, logDiffs.xlog, logDiffs.ulog, logDiffs.hlog, logDiffs.ddlog, upw1, upw2,
               upw_ms, coeffs);

    return coeffs;
}

void BlDiffSolver::bldifLaminar(BoundaryLayerStationWindow &BoundaryLayerStationWindow, double amcrit,
                                BlSystemCoeffs &coeffs) {
    blData &station1 = BoundaryLayerStationWindow.station1;
    blData &station2 = BoundaryLayerStationWindow.station2;

    const auto ax_result = BoundaryLayerUtil::axset(station1.hkz.scalar, station1.param.tz, station1.rtz.scalar,
                                                    station1.param.amplz, station2.hkz.scalar, station2.param.tz,
                                                    station2.rtz.scalar, station2.param.amplz, amcrit);

    const double rezc =
        station2.param.amplz - station1.param.amplz - ax_result.ax * (station2.param.xz - station1.param.xz);
    const double z_ax = -(station2.param.xz - station1.param.xz);

    coeffs.a1(0, 0) = z_ax * ax_result.ax_a1 - 1.0;
    coeffs.a1(0, 1) =
        z_ax * (ax_result.ax_hk1 * station1.hkz.t() + ax_result.ax_t1 + ax_result.ax_rt1 * station1.rtz.t());
    coeffs.a1(0, 2) = z_ax * (ax_result.ax_hk1 * station1.hkz.d());
    coeffs.a1(0, 3) = z_ax * (ax_result.ax_hk1 * station1.hkz.u() + ax_result.ax_rt1 * station1.rtz.u());
    coeffs.a1(0, 4) = ax_result.ax;
    coeffs.a2(0, 0) = z_ax * ax_result.ax_a2 + 1.0;
    coeffs.a2(0, 1) =
        z_ax * (ax_result.ax_hk2 * station2.hkz.t() + ax_result.ax_t2 + ax_result.ax_rt2 * station2.rtz.t());
    coeffs.a2(0, 2) = z_ax * (ax_result.ax_hk2 * station2.hkz.d());
    coeffs.a2(0, 3) = z_ax * (ax_result.ax_hk2 * station2.hkz.u() + ax_result.ax_rt2 * station2.rtz.u());
    coeffs.a2(0, 4) = -ax_result.ax;
    coeffs.d_msq[0] = z_ax * (ax_result.ax_hk1 * station1.hkz.ms() + ax_result.ax_rt1 * station1.rtz.ms() +
                              ax_result.ax_hk2 * station2.hkz.ms() + ax_result.ax_rt2 * station2.rtz.ms());
    coeffs.d_re[0]  = z_ax * (ax_result.ax_rt1 * station1.rtz.re() + ax_result.ax_rt2 * station2.rtz.re());
    coeffs.d_xi[0]  = 0.0;
    coeffs.rhs[0]   = -rezc;
}

void BlDiffSolver::bldifTurbulent(BoundaryLayerStationWindow &BoundaryLayerStationWindow, FlowRegimeEnum flowRegimeType,
                                  double upw, const Vector3d &upw1, const Vector3d &upw2, double upw_ms, double ulog,
                                  BlSystemCoeffs &coeffs) {
    blData &station1 = BoundaryLayerStationWindow.station1;
    blData &station2 = BoundaryLayerStationWindow.station2;

    const double sa  = (1.0 - upw) * station1.param.sz + upw * station2.param.sz;
    const double cqa = (1.0 - upw) * station1.cqz.scalar + upw * station2.cqz.scalar;
    const double cfa = (1.0 - upw) * station1.cfz.scalar + upw * station2.cfz.scalar;
    const double hka = (1.0 - upw) * station1.hkz.scalar + upw * station2.hkz.scalar;

    const double usa = 0.5 * (station1.usz.scalar + station2.usz.scalar);
    const double rta = 0.5 * (station1.rtz.scalar + station2.rtz.scalar);
    const double dea = 0.5 * (station1.dez.scalar + station2.dez.scalar);
    const double da  = 0.5 * (station1.param.dz + station2.param.dz);

    const double ald = (flowRegimeType == FlowRegimeEnum::Wake) ? kDlcon : 1.0;

    double hkc;
    double hkc_hka;
    if (flowRegimeType == FlowRegimeEnum::Turbulent) {
        hkc     = hka - 1.0 - kGccon / rta;
        hkc_hka = 1.0;
        if (hkc < 0.01) {
            hkc     = 0.01;
            hkc_hka = 0.0;
        }
    } else {
        hkc     = hka - 1.0;
        hkc_hka = 1.0;
    }

    const double hr     = hkc / (kGacon * ald * hka);
    const double hr_hka = hkc_hka / (kGacon * ald * hka) - hr / hka;

    const double uq     = (0.5 * cfa - hr * hr) / (kGbcon * da);
    const double uq_hka = -2.0 * hr * hr_hka / (kGbcon * da);
    const double uq_cfa = 0.5 / (kGbcon * da);
    const double uq_da  = -uq / da;

    const double scc     = kSccon * 1.333 / (1.0 + usa);
    const double scc_usa = -scc / (1.0 + usa);

    const double slog = log(station2.param.sz / station1.param.sz);
    const double dxi  = station2.param.xz - station1.param.xz;

    const double rezc = scc * (cqa - sa * ald) * dxi - dea * 2.0 * slog + dea * 2.0 * (uq * dxi - ulog);

    const double z_cfa = dea * 2.0 * uq_cfa * dxi;
    const double z_hka = dea * 2.0 * uq_hka * dxi;
    const double z_da  = dea * 2.0 * uq_da * dxi;
    const double z_sl  = -dea * 2.0;
    const double z_ul  = -dea * 2.0;
    const double z_dxi = scc * (cqa - sa * ald) + dea * 2.0 * uq;
    const double z_usa = scc_usa * (cqa - sa * ald) * dxi;
    const double z_cqa = scc * dxi;
    const double z_sa  = -scc * dxi * ald;
    const double z_dea = 2.0 * (uq * dxi - ulog - slog);

    const double z_upw =
        z_cqa * (station2.cqz.scalar - station1.cqz.scalar) + z_sa * (station2.param.sz - station1.param.sz) +
        z_cfa * (station2.cfz.scalar - station1.cfz.scalar) + z_hka * (station2.hkz.scalar - station1.hkz.scalar);
    const double z_de  = 0.5 * z_dea;
    const double z_us  = 0.5 * z_usa;
    const double z_d   = 0.5 * z_da;
    const double z_u1  = -z_ul / station1.param.uz;
    const double z_u2  = z_ul / station2.param.uz;
    const double z_x1  = -z_dxi;
    const double z_x2  = z_dxi;
    const double z_s1  = (1.0 - upw) * z_sa - z_sl / station1.param.sz;
    const double z_s2  = upw * z_sa + z_sl / station2.param.sz;
    const double z_cq1 = (1.0 - upw) * z_cqa;
    const double z_cq2 = upw * z_cqa;
    const double z_cf1 = (1.0 - upw) * z_cfa;
    const double z_cf2 = upw * z_cfa;
    const double z_hk1 = (1.0 - upw) * z_hka;
    const double z_hk2 = upw * z_hka;

    RowVector<double, 5> vs1_row;
    RowVector<double, 5> vs2_row;

    vs1_row << z_s1, 0.0, z_d, z_u1, z_x1;
    vs2_row << z_s2, 0.0, z_d, z_u2, z_x2;

    const Vector3d vs1_vec = z_upw * upw1 + z_de * station1.dez.pos_vector() + z_us * station1.usz.pos_vector() +
                             z_cq1 * station1.cqz.pos_vector() + z_cf1 * station1.cfz.pos_vector() +
                             z_hk1 * station1.hkz.pos_vector();

    const Vector3d vs2_vec = z_upw * upw2 + z_de * station2.dez.pos_vector() + z_us * station2.usz.pos_vector() +
                             z_cq2 * station2.cqz.pos_vector() + z_cf2 * station2.cfz.pos_vector() +
                             z_hk2 * station2.hkz.pos_vector();

    vs1_row.segment<3>(1) += vs1_vec;
    vs2_row.segment<3>(1) += vs2_vec;

    coeffs.a1.row(0) = vs1_row;
    coeffs.a2.row(0) = vs2_row;

    coeffs.d_msq[0] = z_upw * upw_ms + z_de * station1.dez.ms() + z_us * station1.usz.ms() + z_de * station2.dez.ms() +
                      z_us * station2.usz.ms() + z_cq1 * station1.cqz.ms() + z_cf1 * station1.cfz.ms() +
                      z_hk1 * station1.hkz.ms() + z_cq2 * station2.cqz.ms() + z_cf2 * station2.cfz.ms() +
                      z_hk2 * station2.hkz.ms();

    coeffs.d_re[0] =
        z_cq1 * station1.cqz.re() + z_cf1 * station1.cfz.re() + z_cq2 * station2.cqz.re() + z_cf2 * station2.cfz.re();
    coeffs.d_xi[0] = 0.0;
    coeffs.rhs[0]  = -rezc;
}

void BlDiffSolver::bldifMomentum(BoundaryLayerStationWindow &BoundaryLayerStationWindow, double xlog, double ulog,
                                 double tlog, double ddlog, const SkinFrictionCoefficients &skinFriction,
                                 BlSystemCoeffs &coeffs) {
    blData &station1 = BoundaryLayerStationWindow.station1;
    blData &station2 = BoundaryLayerStationWindow.station2;

    const double ha  = 0.5 * (station1.param.hz + station2.param.hz);
    const double ma  = 0.5 * (station1.param.mz + station2.param.mz);
    const double xa  = 0.5 * (station1.param.xz + station2.param.xz);
    const double ta  = 0.5 * (station1.param.tz + station2.param.tz);
    const double hwa = 0.5 * (station1.param.dwz / station1.param.tz + station2.param.dwz / station2.param.tz);

    const double cfx =
        0.50 * skinFriction.cfm * xa / ta + 0.25 * (station1.cfz.scalar * station1.param.xz / station1.param.tz +
                                                    station2.cfz.scalar * station2.param.xz / station2.param.tz);
    const double cfx_xa = 0.50 * skinFriction.cfm / ta;
    const double cfx_ta = -.50 * skinFriction.cfm * xa / ta / ta;

    const double cfx_x1 = 0.25 * station1.cfz.scalar / station1.param.tz + cfx_xa * 0.5;
    const double cfx_x2 = 0.25 * station2.cfz.scalar / station2.param.tz + cfx_xa * 0.5;
    const double cfx_t1 =
        -.25 * station1.cfz.scalar * station1.param.xz / station1.param.tz / station1.param.tz + cfx_ta * 0.5;
    const double cfx_t2 =
        -.25 * station2.cfz.scalar * station2.param.xz / station2.param.tz / station2.param.tz + cfx_ta * 0.5;
    const double cfx_cf1 = 0.25 * station1.param.xz / station1.param.tz;
    const double cfx_cf2 = 0.25 * station2.param.xz / station2.param.tz;
    const double cfx_cfm = 0.50 * xa / ta;

    const double btmp = ha + 2.0 - ma + hwa;

    const double rezt  = tlog + btmp * ulog - xlog * 0.5 * cfx;
    const double z_cfx = -xlog * 0.5;
    const double z_ha  = ulog;
    const double z_hwa = ulog;
    const double z_ma  = -ulog;
    const double z_xl  = -ddlog * 0.5 * cfx;
    const double z_ul  = ddlog * btmp;
    const double z_tl  = ddlog;

    const double z_cfm = z_cfx * cfx_cfm;
    const double z_cf1 = z_cfx * cfx_cf1;
    const double z_cf2 = z_cfx * cfx_cf2;

    const double z_t1 = -z_tl / station1.param.tz + z_cfx * cfx_t1 +
                        z_hwa * 0.5 * (-station1.param.dwz / station1.param.tz / station1.param.tz);
    const double z_t2 = z_tl / station2.param.tz + z_cfx * cfx_t2 +
                        z_hwa * 0.5 * (-station2.param.dwz / station2.param.tz / station2.param.tz);
    const double z_x1 = -z_xl / station1.param.xz + z_cfx * cfx_x1;
    const double z_x2 = z_xl / station2.param.xz + z_cfx * cfx_x2;
    const double z_u1 = -z_ul / station1.param.uz;
    const double z_u2 = z_ul / station2.param.uz;

    {
        RowVector<double, 5> row1_a1 = RowVector<double, 5>::Zero();
        RowVector<double, 5> row1_a2 = RowVector<double, 5>::Zero();

        const Vector3d hterm1(station1.param.hz_tz, station1.param.hz_dz, 0.0);
        const Vector3d hterm2(station2.param.hz_tz, station2.param.hz_dz, 0.0);
        const Vector3d cfm1v(skinFriction.cfm_t1, skinFriction.cfm_d1, skinFriction.cfm_u1);
        const Vector3d cfm2v(skinFriction.cfm_t2, skinFriction.cfm_d2, skinFriction.cfm_u2);
        const Vector3d cfz1v = station1.cfz.pos_vector();
        const Vector3d cfz2v = station2.cfz.pos_vector();
        const Vector3d mz1(0.0, 0.0, 0.5 * z_ma * station1.param.mz_uz);
        const Vector3d mz2(0.0, 0.0, 0.5 * z_ma * station2.param.mz_uz);

        const Vector3d seg1 = 0.5 * z_ha * hterm1 + z_cfm * cfm1v + z_cf1 * cfz1v + mz1 + Vector3d(z_t1, 0.0, z_u1);
        const Vector3d seg2 = 0.5 * z_ha * hterm2 + z_cfm * cfm2v + z_cf2 * cfz2v + mz2 + Vector3d(z_t2, 0.0, z_u2);

        row1_a1.segment<3>(1) = seg1;
        row1_a2.segment<3>(1) = seg2;
        row1_a1(4)            = z_x1;
        row1_a2(4)            = z_x2;

        coeffs.a1.row(1) = row1_a1;
        coeffs.a2.row(1) = row1_a2;
    }

    coeffs.d_msq[1] = 0.5 * z_ma * station1.param.mz_ms + z_cfm * skinFriction.cfm_ms + z_cf1 * station1.cfz.ms() +
                      0.5 * z_ma * station2.param.mz_ms + z_cf2 * station2.cfz.ms();
    coeffs.d_re[1] = z_cfm * skinFriction.cfm_re + z_cf1 * station1.cfz.re() + z_cf2 * station2.cfz.re();
    coeffs.d_xi[1] = 0.0;
    coeffs.rhs[1]  = -rezt;
}

void BlDiffSolver::bldifShape(BoundaryLayerStationWindow &BoundaryLayerStationWindow, double upw, double xlog,
                              double ulog, double hlog, double ddlog, const Vector3d &upw1, const Vector3d &upw2,
                              double upw_ms, BlSystemCoeffs &coeffs) {
    blData &station1 = BoundaryLayerStationWindow.station1;
    blData &station2 = BoundaryLayerStationWindow.station2;

    const double xot1 = station1.param.xz / station1.param.tz;
    const double xot2 = station2.param.xz / station2.param.tz;

    const double ha  = 0.5 * (station1.param.hz + station2.param.hz);
    const double hsa = 0.5 * (station1.hsz.scalar + station2.hsz.scalar);
    const double hca = 0.5 * (station1.hcz.scalar + station2.hcz.scalar);
    const double hwa = 0.5 * (station1.param.dwz / station1.param.tz + station2.param.dwz / station2.param.tz);

    const double dix     = (1.0 - upw) * station1.diz.scalar * xot1 + upw * station2.diz.scalar * xot2;
    const double cfx     = (1.0 - upw) * station1.cfz.scalar * xot1 + upw * station2.cfz.scalar * xot2;
    const double dix_upw = station2.diz.scalar * xot2 - station1.diz.scalar * xot1;
    const double cfx_upw = station2.cfz.scalar * xot2 - station1.cfz.scalar * xot1;

    const double btmp = 2.0 * hca / hsa + 1.0 - ha - hwa;

    const double rezh  = hlog + btmp * ulog + xlog * (0.5 * cfx - dix);
    const double z_cfx = xlog * 0.5;
    const double z_dix = -xlog;
    const double z_hca = 2.0 * ulog / hsa;
    const double z_ha  = -ulog;
    const double z_hwa = -ulog;
    const double z_xl  = ddlog * (0.5 * cfx - dix);
    const double z_ul  = ddlog * btmp;
    const double z_hl  = ddlog;

    const double z_upw = z_cfx * cfx_upw + z_dix * dix_upw;

    const double z_hs1 = -hca * ulog / hsa / hsa - z_hl / station1.hsz.scalar;
    const double z_hs2 = -hca * ulog / hsa / hsa + z_hl / station2.hsz.scalar;

    const double z_cf1 = (1.0 - upw) * z_cfx * xot1;
    const double z_cf2 = upw * z_cfx * xot2;
    const double z_di1 = (1.0 - upw) * z_dix * xot1;
    const double z_di2 = upw * z_dix * xot2;

    double z_t1 =
        (1.0 - upw) * (z_cfx * station1.cfz.scalar + z_dix * station1.diz.scalar) * (-xot1 / station1.param.tz);
    double z_t2       = upw * (z_cfx * station2.cfz.scalar + z_dix * station2.diz.scalar) * (-xot2 / station2.param.tz);
    const double z_x1 = (1.0 - upw) * (z_cfx * station1.cfz.scalar + z_dix * station1.diz.scalar) / station1.param.tz -
                        z_xl / station1.param.xz;
    const double z_x2 = upw * (z_cfx * station2.cfz.scalar + z_dix * station2.diz.scalar) / station2.param.tz +
                        z_xl / station2.param.xz;
    const double z_u1 = -z_ul / station1.param.uz;
    const double z_u2 = z_ul / station2.param.uz;

    z_t1 += z_hwa * 0.5 * (-station1.param.dwz / station1.param.tz / station1.param.tz);
    z_t2 += z_hwa * 0.5 * (-station2.param.dwz / station2.param.tz / station2.param.tz);

    {
        RowVector<double, 5> row2_a1 = RowVector<double, 5>::Zero();
        RowVector<double, 5> row2_a2 = RowVector<double, 5>::Zero();

        const Vector3d base1 = z_hs1 * station1.hsz.pos_vector() + z_cf1 * station1.cfz.pos_vector() +
                               z_di1 * station1.diz.pos_vector() + Vector3d(z_t1, 0.0, z_u1);
        const Vector3d base2 = z_hs2 * station2.hsz.pos_vector() + z_cf2 * station2.cfz.pos_vector() +
                               z_di2 * station2.diz.pos_vector() + Vector3d(z_t2, 0.0, z_u2);

        row2_a1(0)            = z_di1 * station1.diz.s();
        row2_a2(0)            = z_di2 * station2.diz.s();
        row2_a1.segment<3>(1) = base1;
        row2_a2.segment<3>(1) = base2;
        row2_a1(4)            = z_x1;
        row2_a2(4)            = z_x2;

        coeffs.a1.row(2) = row2_a1;
        coeffs.a2.row(2) = row2_a2;
    }
    coeffs.d_msq[2] = z_hs1 * station1.hsz.ms() + z_cf1 * station1.cfz.ms() + z_di1 * station1.diz.ms() +
                      z_hs2 * station2.hsz.ms() + z_cf2 * station2.cfz.ms() + z_di2 * station2.diz.ms();
    coeffs.d_re[2] = z_hs1 * station1.hsz.re() + z_cf1 * station1.cfz.re() + z_di1 * station1.diz.re() +
                     z_hs2 * station2.hsz.re() + z_cf2 * station2.cfz.re() + z_di2 * station2.diz.re();

    {
        const Vector3d inc1 = 0.5 * (z_hca * station1.hcz.pos_vector() +
                                     z_ha * Vector3d(station1.param.hz_tz, station1.param.hz_dz, 0.0)) +
                              z_upw * upw1;
        const Vector3d inc2 = 0.5 * (z_hca * station2.hcz.pos_vector() +
                                     z_ha * Vector3d(station2.param.hz_tz, station2.param.hz_dz, 0.0)) +
                              z_upw * upw2;
        coeffs.a1.row(2).segment<3>(1) += inc1;
        coeffs.a2.row(2).segment<3>(1) += inc2;
    }

    coeffs.d_msq[2] = 0.5 * (z_hca * station1.hcz.ms()) + z_upw * upw_ms + 0.5 * (z_hca * station2.hcz.ms());

    coeffs.d_xi[2] = 0.0;
    coeffs.rhs[2]  = -rezh;
}
