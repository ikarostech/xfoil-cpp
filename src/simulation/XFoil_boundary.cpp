/**
 * Boundary-layer scalar computations split from XFoil.cpp
 */

#include "XFoil.h"
// dependencies are included by XFoil.h except coefficient headers
#include "domain/coefficient/skin_friction.hpp"
#include "domain/coefficient/dissipation.hpp"

using Eigen::Vector3d;

namespace {

skin_friction::C_f evaluateSkinFriction(const blData &data,
                                        FlowRegimeEnum flowRegimeType) {
  return skin_friction::getSkinFriction(data.hkz.scalar, data.rtz.scalar,
                                        data.param.mz, flowRegimeType);
}

dissipation::DissipationResult evaluateDissipation(const blData &data,
                                                   FlowRegimeEnum flowRegimeType) {
  return dissipation::getDissipation(data.hkz.scalar, data.rtz.scalar,
                                     flowRegimeType);
}

}  // namespace

blData XFoil::computeShapeParameters(const blData &ref,
                                     FlowRegimeEnum flowRegimeType) const {
  blData result = ref;

  if (flowRegimeType == FlowRegimeEnum::Wake)
    result.hkz.scalar = std::max(result.hkz.scalar, 1.00005);
  if (flowRegimeType != FlowRegimeEnum::Wake)
    result.hkz.scalar = std::max(result.hkz.scalar, 1.05000);

  const auto hct_result = boundary_layer::hct(result.hkz.scalar, result.param.mz);
  result.hcz.scalar = hct_result.hc;
  result.hcz.vector.segment<3>(0) =
      hct_result.hc_hk * result.hkz.vector.segment<3>(0);
  result.hcz.u() += hct_result.hc_msq * result.param.mz_uz;
  result.hcz.ms() += hct_result.hc_msq * result.param.mz_ms;

  boundary_layer::ThicknessShapeParameterResult hs_result;
  if (flowRegimeType == FlowRegimeEnum::Laminar) {
    hs_result = boundary_layer::hsl(result.hkz.scalar);
    result.hsz.scalar = hs_result.hs;
  } else {
    hs_result =
        boundary_layer::hst(result.hkz.scalar, result.rtz.scalar, result.param.mz);
    result.hsz.scalar = hs_result.hs;
  }

  result.hsz.vector =
      hs_result.hs_hk * result.hkz.vector + hs_result.hs_rt * result.rtz.vector;
  result.hsz.u() += hs_result.hs_msq * result.param.mz_uz;
  result.hsz.ms() += hs_result.hs_msq * result.param.mz_ms;

  const double us2_hs2 = 0.5 *
      (1.0 - (result.hkz.scalar - 1.0) / (gbcon * result.param.hz));
  const double us2_hk2 =
      0.5 * result.hsz.scalar * (-1.0 / (gbcon * result.param.hz));
  const double us2_h2 = 0.5 * result.hsz.scalar * (result.hkz.scalar - 1.0) /
                        (gbcon * result.param.hz * result.param.hz);

  result.usz.scalar = 0.5 * result.hsz.scalar *
                      (1.0 - (result.hkz.scalar - 1.0) /
                                   (gbcon * result.param.hz));
  result.usz.vector =
      us2_hs2 * result.hsz.vector + us2_hk2 * result.hkz.vector;
  result.usz.t() += us2_h2 * result.param.hz_tz;
  result.usz.d() += us2_h2 * result.param.hz_dz;

  if (flowRegimeType == FlowRegimeEnum::Laminar ||
      flowRegimeType == FlowRegimeEnum::Turbulent) {
    if (result.usz.scalar > 0.95) {
      result.usz.scalar = 0.98;
      result.usz.vector = Vector<double, 6>::Zero();
    }
  }

  if (flowRegimeType == FlowRegimeEnum::Wake && result.usz.scalar > 0.99995) {
    result.usz.scalar = 0.99995;
    result.usz.vector = Vector<double, 6>::Zero();
  }

  return result;
}

blData XFoil::computeShearCoefficients(const blData &ref,
                                       FlowRegimeEnum flowRegimeType) const {
  blData result = ref;

  double hkc = result.hkz.scalar - 1.0;
  double hkc_hk2 = 1.0;
  double hkc_rt2 = 0.0;
  if (flowRegimeType == FlowRegimeEnum::Turbulent) {
    const double gcc = gccon;
    hkc = result.hkz.scalar - 1.0 - gcc / result.rtz.scalar;
    hkc_hk2 = 1.0;
    hkc_rt2 = gcc / result.rtz.scalar / result.rtz.scalar;
    if (hkc < 0.01) {
      hkc = 0.01;
      hkc_hk2 = 0.0;
      hkc_rt2 = 0.0;
    }
  }

  double hkb = result.hkz.scalar - 1.0;
  double usb = 1.0 - result.usz.scalar;
  result.cqz.scalar =
      hkc / result.hkz.scalar *
      sqrt(ctcon * result.hsz.scalar * hkb / (usb * result.param.hz));
  double cq2_hs2 = result.cqz.scalar / result.hsz.scalar / 2;
  double cq2_us2 = result.cqz.scalar * (0.5 / usb);
  double cq2_hk2 =
      result.cqz.scalar *
      (0.5 / hkb - 1.0 / result.hkz.scalar + hkc_hk2 / hkc);
  double cq2_rt2 = result.cqz.scalar * (hkc_rt2 / hkc);
  double cq2_h2 = result.cqz.scalar * (-0.5 / result.param.hz);

  result.cqz.vector = cq2_hs2 * result.hsz.vector +
                      cq2_us2 * result.usz.vector +
                      cq2_hk2 * result.hkz.vector +
                      cq2_rt2 * result.rtz.vector;
  result.cqz.t() += cq2_h2 * result.param.hz_tz;
  result.cqz.d() += cq2_h2 * result.param.hz_dz;

  return result;
}

blData XFoil::computeSkinFrictionCoefficients(
    const blData &ref, FlowRegimeEnum flowRegimeType) const {
  blData result = ref;

  skin_friction::C_f c_f = evaluateSkinFriction(result, flowRegimeType);
  result.cfz.scalar = c_f.cf;
  double cf2_hk2 = c_f.hk;
  double cf2_rt2 = c_f.rt;
  double cf2_m2 = c_f.msq;

  result.cfz.vector =
      cf2_hk2 * result.hkz.vector + cf2_rt2 * result.rtz.vector;
  result.cfz.u() += cf2_m2 * result.param.mz_uz;
  result.cfz.ms() += cf2_m2 * result.param.mz_ms;

  return result;
}

blData XFoil::computeDissipation(const blData &ref,
                                 FlowRegimeEnum flowRegimeType) const {
  blData result = ref;
  if (flowRegimeType == FlowRegimeEnum::Laminar) {
    auto dissipation_result = evaluateDissipation(result, flowRegimeType);
    result.diz.scalar = dissipation_result.di;
    result.diz.vector = dissipation_result.di_hk * result.hkz.vector +
                        dissipation_result.di_rt * result.rtz.vector;
  } else {
    if (flowRegimeType == FlowRegimeEnum::Turbulent) {
      auto c_ft = evaluateSkinFriction(result, flowRegimeType);
      double cf2t = c_ft.cf;
      double cf2t_hk2 = c_ft.hk;
      double cf2t_rt2 = c_ft.rt;
      double cf2t_m2 = c_ft.msq;
      Vector<double, 6> cf2t_vector =
          cf2t_hk2 * result.hkz.vector + cf2t_rt2 * result.rtz.vector +
          Vector<double, 6>{cf2t_m2 * result.param.mz_uz, 0, 0, 0,
                            cf2t_m2 * result.param.mz_ms, 0};

      result.diz.scalar =
          (0.5 * cf2t * result.usz.scalar) * 2.0 / result.hsz.scalar;
      double di2_hs2 =
          -(0.5 * cf2t * result.usz.scalar) * 2.0 / result.hsz.scalar /
          result.hsz.scalar;
      double di2_us2 = (0.5 * cf2t) * 2.0 / result.hsz.scalar;
      double di2_cf2t = (0.5 * result.usz.scalar) * 2.0 / result.hsz.scalar;

      result.diz.vector = di2_hs2 * result.hsz.vector +
                          di2_us2 * result.usz.vector +
                          di2_cf2t * cf2t_vector;

      double grt = log(result.rtz.scalar);
      double hmin = 1.0 + 2.1 / grt;
      double hm_rt2 = -(2.1 / grt / grt) / result.rtz.scalar;

      double fl = (result.hkz.scalar - 1.0) / (hmin - 1.0);
      double fl_hk2 = 1.0 / (hmin - 1.0);
      double fl_rt2 = (-fl / (hmin - 1.0)) * hm_rt2;

      double tfl = tanh(fl);
      double dfac = 0.5 + 0.5 * tfl;
      double df_fl = 0.5 * (1.0 - tfl * tfl);

      double df_hk2 = df_fl * fl_hk2;
      double df_rt2 = df_fl * fl_rt2;

      result.diz.vector =
          dfac * result.diz.vector +
          result.diz.scalar *
              (df_hk2 * result.hkz.vector + df_rt2 * result.rtz.vector);
      result.diz.re() -= df_rt2 * result.rtz.re();

      result.diz.scalar = result.diz.scalar * dfac;
    } else {
      result.diz.scalar = 0.0;
      result.diz.vector = Vector<double, 6>::Zero();
    }
  }

  if (flowRegimeType != FlowRegimeEnum::Laminar) {
    double dd = result.param.sz * result.param.sz *
                (0.995 - result.usz.scalar) * 2.0 / result.hsz.scalar;
    double dd_hs2 = -result.param.sz * result.param.sz *
                    (0.995 - result.usz.scalar) * 2.0 /
                    result.hsz.scalar / result.hsz.scalar;
    double dd_us2 =
        -result.param.sz * result.param.sz * 2.0 / result.hsz.scalar;
    const double dd_s2 =
        result.param.sz * 2.0 * (0.995 - result.usz.scalar) * 2.0 /
        result.hsz.scalar;

    result.diz.scalar = result.diz.scalar + dd;
    result.diz.s() = dd_s2;
    result.diz.vector =
        result.diz.vector + dd_hs2 * result.hsz.vector +
        dd_us2 * result.usz.vector;

    dd = 0.15 * (0.995 - result.usz.scalar) * (0.995 - result.usz.scalar) /
         result.rtz.scalar * 2.0 / result.hsz.scalar;
    dd_us2 = -0.15 * (0.995 - result.usz.scalar) * 2.0 /
             result.rtz.scalar * 2.0 / result.hsz.scalar;
    dd_hs2 = -dd / result.hsz.scalar;
    const double dd_rt2 = -dd / result.rtz.scalar;

    result.diz.scalar = result.diz.scalar + dd;
    result.diz.vector = result.diz.vector + dd_hs2 * result.hsz.vector +
                        dd_us2 * result.usz.vector +
                        dd_rt2 * result.rtz.vector;
  }

  if (flowRegimeType == FlowRegimeEnum::Turbulent) {
    auto dissipation_result = evaluateDissipation(result, flowRegimeType);
    if (dissipation_result.di > result.diz.scalar) {
      result.diz.scalar = dissipation_result.di;
      result.diz.vector = dissipation_result.di_hk * result.hkz.vector +
                          dissipation_result.di_rt * result.rtz.vector;
    }
  }

  if (flowRegimeType == FlowRegimeEnum::Wake) {
    auto dissipation_result = evaluateDissipation(result, flowRegimeType);
    const double di2l = dissipation_result.di;
    if (di2l > result.diz.scalar) {
      result.diz.scalar = dissipation_result.di;
      result.diz.vector = dissipation_result.di_hk * result.hkz.vector +
                          dissipation_result.di_rt * result.rtz.vector;
    }

    result.diz.scalar = result.diz.scalar * 2.0;
    result.diz.vector *= 2;
  }

  return result;
}

blData XFoil::computeThickness(const blData &ref,
                               [[maybe_unused]] FlowRegimeEnum flowRegimeType) const {
  blData result = ref;

  result.dez.scalar = (3.15 + 1.72 / (result.hkz.scalar - 1.0)) *
                      result.param.tz + result.param.dz;
  double de2_hk2 = (-1.72 / (result.hkz.scalar - 1.0) /
                    (result.hkz.scalar - 1.0)) * result.param.tz;
  result.dez.vector = de2_hk2 * result.hkz.vector;
  result.dez.t() += (3.15 + 1.72 / (result.hkz.scalar - 1.0));
  result.dez.d() += 1.0;

  double hdmax = 12.0;
  if (result.dez.scalar > hdmax * result.param.tz) {
    result.dez.scalar = hdmax * result.param.tz;
    result.dez.vector = Vector<double, 6>::Zero();
    result.dez.t() = hdmax;
  }

  return result;
}

bool XFoil::dslim(double &dstr, double thet, double msq, double hklim) {
  const double h = (dstr) / thet;

  boundary_layer::KineticShapeParameterResult hkin_result =
      boundary_layer::hkin(h, msq);

  const double dh = std::max(0.0, hklim - hkin_result.hk) / hkin_result.hk_h;
  dstr = (dstr) + dh * thet;

  return true;
}
