#pragma once

#include <Eigen/Core>

/**
 * @brief Internal station-level boundary-layer state owned by the model.
 *
 * This is not a generic solver scratch buffer. It is the model's local state
 * representation for one marching station, including both primary values and
 * the derived sensitivities that are repeatedly reused by boundary-layer
 * computations.
 */
class BoundaryLayerStationState {
public:
  class StateVector {
  public:
    double scalar = 0.0;
    Eigen::Vector<double, 6> vector = Eigen::Vector<double, 6>::Zero();
    Eigen::Vector<double, 3> pos_vector() { return vector.segment(0, 3); }
    Eigen::Vector<double, 3> pos_vector() const { return vector.segment(0, 3); }
    double &t() { return vector[0]; }
    double &d() { return vector[1]; }
    double &u() { return vector[2]; }
    double &s() { return vector[3]; }
    double &ms() { return vector[4]; }
    double &re() { return vector[5]; }
  };

  class PrimaryState {
  public:
    double xz = 0.0;
    double uz = 0.0;
    double tz = 0.0;
    double dz = 0.0;
    double sz = 0.0;
    double amplz = 0.0;
    double uz_uei = 0.0;
    double uz_ms = 0.0;
    double dwz = 0.0;
    double hz = 0.0;
    double hz_tz = 0.0;
    double hz_dz = 0.0;
    double mz = 0.0;
    double mz_uz = 0.0;
    double mz_ms = 0.0;
    double rz = 0.0;
    double rz_uz = 0.0;
    double rz_ms = 0.0;
  };

  using blVector = StateVector;
  using blParam = PrimaryState;

  StateVector hkz;
  StateVector hsz;
  StateVector hcz;
  StateVector rtz;
  StateVector cfz;
  StateVector diz;
  StateVector usz;
  StateVector cqz;
  StateVector dez;
  PrimaryState param;

  void assignPrimaryStationData(double xsi, double amplification,
                                double shearCoefficient,
                                double momentumThickness,
                                double displacementThickness,
                                double wakeGapThickness) {
    param.xz = xsi;
    param.amplz = amplification;
    param.sz = shearCoefficient;
    param.tz = momentumThickness;
    param.dz = displacementThickness - wakeGapThickness;
    param.dwz = wakeGapThickness;
  }

  void assignCompressibleEdgeVelocity(double edgeVelocity, double qinf,
                                      double karmanTsienFactor,
                                      double karmanTsienFactor_msq) {
    const double velocity_ratio = edgeVelocity / qinf;
    const double denom =
        1.0 - karmanTsienFactor * velocity_ratio * velocity_ratio;
    param.uz = edgeVelocity * (1.0 - karmanTsienFactor) / denom;
    param.uz_uei =
        (1.0 + karmanTsienFactor *
                   (2.0 * param.uz * velocity_ratio / qinf - 1.0)) /
        denom;
    param.uz_ms =
        (param.uz * velocity_ratio * velocity_ratio - edgeVelocity) *
        karmanTsienFactor_msq / denom;
  }
};

/**
 * @brief Internal sensitivity state owned by the boundary-layer model.
 *
 * Transition and marching operations use this as a compact representation of
 * local sensitivities. Keeping it in the model clarifies that callers should
 * depend on the model behavior, not on the raw storage layout.
 */
class BoundaryLayerSensitivityState {
public:
  Eigen::Vector<double, 13> vector = Eigen::Vector<double, 13>::Zero();
  double scalar = 0.0;
  double &a1() { return vector[0]; }
  double &t1() { return vector[1]; }
  double &d1() { return vector[2]; }
  double &u1() { return vector[3]; }
  double &x1() { return vector[4]; }
  double &t2() { return vector[5]; }
  double &d2() { return vector[6]; }
  double &u2() { return vector[7]; }
  double &x2() { return vector[8]; }
  double &a2() { return vector[9]; }
  double &ms() { return vector[10]; }
  double &re() { return vector[11]; }
  double &xf() { return vector[12]; }
};

using blData = BoundaryLayerStationState;
using blDiff = BoundaryLayerSensitivityState;

class boundary_layer {
public:
  // TODO 副作用の除去
  class DensityShapeParameterResult {
  public:
    double hc;
    double hc_hk;
    double hc_msq;
  };
  static DensityShapeParameterResult hct(double hk, double msq);

  class KineticShapeParameterResult {
  public:
    double hk;
    double hk_h;
    double hk_msq;
  };
  static KineticShapeParameterResult hkin(double h, double msq);

  class ThicknessShapeParameterResult {
  public:
    double hs;
    double hs_hk;
    double hs_rt;
    double hs_msq;
  };
  static ThicknessShapeParameterResult hsl(double hk);
  static ThicknessShapeParameterResult hst(double hk, double rt, double msq);
};
