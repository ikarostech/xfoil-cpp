/**
 * Boundary-layer scalar computations split from XFoil.cpp
 */

#include "XFoil.h"
// dependencies are included by XFoil.h except coefficient headers
#include "model/coefficient/skin_friction.hpp"
#include "model/coefficient/dissipation.hpp"

using Eigen::Vector3d;

void XFoil::computeShapeParameters(blData &ref, FlowRegimeEnum flowRegimeType) {
  if (flowRegimeType == FlowRegimeEnum::Wake)
    ref.hkz.scalar = std::max(ref.hkz.scalar, 1.00005);
  if (flowRegimeType != FlowRegimeEnum::Wake)
    ref.hkz.scalar = std::max(ref.hkz.scalar, 1.05000);

  auto hct_result = boundary_layer::hct(ref.hkz.scalar, ref.param.mz);
  ref.hcz.scalar = hct_result.hc;
  ref.hcz.pos_vector() = hct_result.hc_hk * ref.hkz.pos_vector();
  ref.hcz.u() += hct_result.hc_msq * ref.param.mz_uz;
  ref.hcz.ms() += hct_result.hc_msq * ref.param.mz_ms;

  boundary_layer::ThicknessShapeParameterResult hs_result;
  if (flowRegimeType == FlowRegimeEnum::Laminar) {
    hs_result = boundary_layer::hsl(ref.hkz.scalar);
    ref.hsz.scalar = hs_result.hs;
  } else {
    hs_result =
        boundary_layer::hst(ref.hkz.scalar, ref.rtz.scalar, ref.param.mz);
    ref.hsz.scalar = hs_result.hs;
  }

  ref.hsz.vector =
      hs_result.hs_hk * ref.hkz.vector + hs_result.hs_rt * ref.rtz.vector;
  ref.hsz.u() += hs_result.hs_msq * ref.param.mz_uz;
  ref.hsz.ms() += hs_result.hs_msq * ref.param.mz_ms;

  double us2_hs2 =
      0.5 * (1.0 - (ref.hkz.scalar - 1.0) / (gbcon * ref.param.hz));
  double us2_hk2 = 0.5 * ref.hsz.scalar * (-1.0 / (gbcon * ref.param.hz));
  double us2_h2 = 0.5 * ref.hsz.scalar * (ref.hkz.scalar - 1.0) /
                  (gbcon * ref.param.hz * ref.param.hz);

  ref.usz.scalar = 0.5 * ref.hsz.scalar *
                   (1.0 - (ref.hkz.scalar - 1.0) / (gbcon * ref.param.hz));
  ref.usz.vector = us2_hs2 * ref.hsz.vector + us2_hk2 * ref.hkz.vector;
  ref.usz.t() += us2_h2 * ref.param.hz_tz;
  ref.usz.d() += us2_h2 * ref.param.hz_dz;

  if (flowRegimeType == FlowRegimeEnum::Laminar ||
      flowRegimeType == FlowRegimeEnum::Turbulent) {
    if (ref.usz.scalar > 0.95) {
      ref.usz.scalar = 0.98;
      ref.usz.vector = Vector<double, 6>::Zero();
    }
  }

  if (flowRegimeType == FlowRegimeEnum::Wake && ref.usz.scalar > 0.99995) {
    ref.usz.scalar = 0.99995;
    ref.usz.vector = Vector<double, 6>::Zero();
  }
}

void XFoil::computeCoefficients(blData &ref, FlowRegimeEnum flowRegimeType) {
  double hkc = ref.hkz.scalar - 1.0;
  double hkc_hk2 = 1.0;
  double hkc_rt2 = 0.0;
  if (flowRegimeType == FlowRegimeEnum::Turbulent) {
    const double gcc = gccon;
    hkc = ref.hkz.scalar - 1.0 - gcc / ref.rtz.scalar;
    hkc_hk2 = 1.0;
    hkc_rt2 = gcc / ref.rtz.scalar / ref.rtz.scalar;
    if (hkc < 0.01) {
      hkc = 0.01;
      hkc_hk2 = 0.0;
      hkc_rt2 = 0.0;
    }
  }

  double hkb = ref.hkz.scalar - 1.0;
  double usb = 1.0 - ref.usz.scalar;
  ref.cqz.scalar = hkc / ref.hkz.scalar *
                   sqrt(ctcon * ref.hsz.scalar * hkb / (usb * ref.param.hz));
  double cq2_hs2 = ref.cqz.scalar / ref.hsz.scalar / 2;
  double cq2_us2 = ref.cqz.scalar * (0.5 / usb);
  double cq2_hk2 =
      ref.cqz.scalar * (0.5 / hkb - 1.0 / ref.hkz.scalar + hkc_hk2 / hkc);
  double cq2_rt2 = ref.cqz.scalar * (hkc_rt2 / hkc);
  double cq2_h2 = ref.cqz.scalar * (-0.5 / ref.param.hz);

  ref.cqz.vector = cq2_hs2 * ref.hsz.vector + cq2_us2 * ref.usz.vector +
                   cq2_hk2 * ref.hkz.vector + cq2_rt2 * ref.rtz.vector;
  ref.cqz.t() += cq2_h2 * ref.param.hz_tz;
  ref.cqz.d() += cq2_h2 * ref.param.hz_dz;

  skin_friction::C_f c_f = skin_friction::getSkinFriction(ref.hkz.scalar, ref.rtz.scalar, ref.param.mz, flowRegimeType);
  ref.cfz.scalar = c_f.cf;
  double cf2_hk2 = c_f.hk;
  double cf2_rt2 = c_f.rt;
  double cf2_m2 = c_f.msq;
  
  ref.cfz.vector = cf2_hk2 * ref.hkz.vector + cf2_rt2 * ref.rtz.vector;
  ref.cfz.u() += cf2_m2 * ref.param.mz_uz;
  ref.cfz.ms() += cf2_m2 * ref.param.mz_ms;
}

void XFoil::computeDissipationAndThickness(blData &ref,
                                           FlowRegimeEnum flowRegimeType) {
  double di2l;
  if (flowRegimeType == FlowRegimeEnum::Laminar) {
    auto dissipation_result = dissipation::getDissipation(
        ref.hkz.scalar, ref.rtz.scalar, flowRegimeType);
    ref.diz.scalar = dissipation_result.di;
    ref.diz.vector = dissipation_result.di_hk * ref.hkz.vector +
                     dissipation_result.di_rt * ref.rtz.vector;
  } else {
    if (flowRegimeType == FlowRegimeEnum::Turbulent) {
      auto c_ft = skin_friction::getSkinFriction(ref.hkz.scalar, ref.rtz.scalar, ref.param.mz, flowRegimeType);
      double cf2t = c_ft.cf;
      double cf2t_hk2 = c_ft.hk;
      double cf2t_rt2 = c_ft.rt;
      double cf2t_m2 = c_ft.msq;
      Vector<double, 6> cf2t_vector =
          cf2t_hk2 * ref.hkz.vector + cf2t_rt2 * ref.rtz.vector +
          Vector<double, 6>{cf2t_m2 * ref.param.mz_uz, 0, 0, 0,
                            cf2t_m2 * ref.param.mz_ms, 0};

      ref.diz.scalar = (0.5 * cf2t * ref.usz.scalar) * 2.0 / ref.hsz.scalar;
      double di2_hs2 = -(0.5 * cf2t * ref.usz.scalar) * 2.0 / ref.hsz.scalar /
                       ref.hsz.scalar;
      double di2_us2 = (0.5 * cf2t) * 2.0 / ref.hsz.scalar;
      double di2_cf2t = (0.5 * ref.usz.scalar) * 2.0 / ref.hsz.scalar;

      ref.diz.vector = di2_hs2 * ref.hsz.vector + di2_us2 * ref.usz.vector +
                       di2_cf2t * cf2t_vector;

      double grt = log(ref.rtz.scalar);
      double hmin = 1.0 + 2.1 / grt;
      double hm_rt2 = -(2.1 / grt / grt) / ref.rtz.scalar;

      double fl = (ref.hkz.scalar - 1.0) / (hmin - 1.0);
      double fl_hk2 = 1.0 / (hmin - 1.0);
      double fl_rt2 = (-fl / (hmin - 1.0)) * hm_rt2;

      double tfl = tanh(fl);
      double dfac = 0.5 + 0.5 * tfl;
      double df_fl = 0.5 * (1.0 - tfl * tfl);

      double df_hk2 = df_fl * fl_hk2;
      double df_rt2 = df_fl * fl_rt2;

      ref.diz.vector =
          dfac * ref.diz.vector +
          ref.diz.scalar * (df_hk2 * ref.hkz.vector + df_rt2 * ref.rtz.vector);
      ref.diz.re() -= df_rt2 * ref.rtz.re();

      ref.diz.scalar = ref.diz.scalar * dfac;
    } else {
      ref.diz.scalar = 0.0;
      ref.diz.vector = Vector<double, 6>::Zero();
    }
  }

  if (flowRegimeType != FlowRegimeEnum::Laminar) {
    double dd = ref.param.sz * ref.param.sz * (0.995 - ref.usz.scalar) * 2.0 /
                ref.hsz.scalar;
    double dd_hs2 = -ref.param.sz * ref.param.sz * (0.995 - ref.usz.scalar) *
                    2.0 / ref.hsz.scalar / ref.hsz.scalar;
    double dd_us2 = -ref.param.sz * ref.param.sz * 2.0 / ref.hsz.scalar;
    const double dd_s2 =
        ref.param.sz * 2.0 * (0.995 - ref.usz.scalar) * 2.0 / ref.hsz.scalar;

    ref.diz.scalar = ref.diz.scalar + dd;
    ref.diz.s() = dd_s2;
    ref.diz.vector =
        ref.diz.vector + dd_hs2 * ref.hsz.vector + dd_us2 * ref.usz.vector;

    dd = 0.15 * (0.995 - ref.usz.scalar) * (0.995 - ref.usz.scalar) /
         ref.rtz.scalar * 2.0 / ref.hsz.scalar;
    dd_us2 = -0.15 * (0.995 - ref.usz.scalar) * 2.0 / ref.rtz.scalar * 2.0 /
             ref.hsz.scalar;
    dd_hs2 = -dd / ref.hsz.scalar;
    const double dd_rt2 = -dd / ref.rtz.scalar;

    ref.diz.scalar = ref.diz.scalar + dd;
    ref.diz.vector = ref.diz.vector + dd_hs2 * ref.hsz.vector +
                     dd_us2 * ref.usz.vector + dd_rt2 * ref.rtz.vector;
  }

  if (flowRegimeType == FlowRegimeEnum::Turbulent) {
    auto dissipation_result = dissipation::getDissipation(
        ref.hkz.scalar, ref.rtz.scalar, flowRegimeType);
    if (dissipation_result.di > ref.diz.scalar) {
      ref.diz.scalar = dissipation_result.di;
      ref.diz.vector = dissipation_result.di_hk * ref.hkz.vector +
                       dissipation_result.di_rt * ref.rtz.vector;
    }
  }

  if (flowRegimeType == FlowRegimeEnum::Wake) {
    auto dissipation_result = dissipation::getDissipation(
        ref.hkz.scalar, ref.rtz.scalar, flowRegimeType);
    di2l = dissipation_result.di;
    if (di2l > ref.diz.scalar) {
      ref.diz.scalar = dissipation_result.di;
      ref.diz.vector = dissipation_result.di_hk * ref.hkz.vector +
                       dissipation_result.di_rt * ref.rtz.vector;
    }

    ref.diz.scalar = ref.diz.scalar * 2.0;
    ref.diz.vector *= 2;
  }

  ref.dez.scalar =
      (3.15 + 1.72 / (ref.hkz.scalar - 1.0)) * ref.param.tz + ref.param.dz;
  double de2_hk2 =
      (-1.72 / (ref.hkz.scalar - 1.0) / (ref.hkz.scalar - 1.0)) * ref.param.tz;
  ref.dez.vector = de2_hk2 * ref.hkz.vector;
  ref.dez.t() += (3.15 + 1.72 / (ref.hkz.scalar - 1.0));
  ref.dez.d() += 1.0;

  double hdmax = 12.0;
  if (ref.dez.scalar > hdmax * ref.param.tz) {
    ref.dez.scalar = hdmax * ref.param.tz;
    ref.dez.vector = Vector<double, 6>::Zero();
    ref.dez.t() = hdmax;
  }
}

bool XFoil::dslim(double &dstr, double thet, double msq, double hklim) {
  const double h = (dstr) / thet;

  boundary_layer::KineticShapeParameterResult hkin_result =
      boundary_layer::hkin(h, msq);

  const double dh = std::max(0.0, hklim - hkin_result.hk) / hkin_result.hk_h;
  dstr = (dstr) + dh * thet;

  return true;
}
