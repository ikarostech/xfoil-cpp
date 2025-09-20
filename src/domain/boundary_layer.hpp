#pragma once
#include <Eigen/Core>

class blData {
    public:
  class blVector {
  public:
    double scalar;
    Eigen::Vector<double, 6> vector = Eigen::Vector<double, 6>::Zero();
    Eigen::Vector<double, 3> pos_vector() {
      return vector.segment(0, 3);
    }
    double& t() {
      return vector[0];
    }
    double& d() {
      return vector[1];
    }
    double& u() {
      return vector[2];
    }
    double& s() {
      return vector[3];
    }
    double& ms() {
      return vector[4];
    }
    double& re() {
      return vector[5];
    }
  };
  blVector hkz, hsz, hcz, rtz, cfz, diz, usz, cqz, dez;
  class blParam {
    public:
    double xz, uz, tz, dz, sz, amplz, uz_uei, uz_ms, dwz, 
      hz, hz_tz, hz_dz, 
      mz, mz_uz, mz_ms, 
      rz, rz_uz, rz_ms;
  };
  blParam param; /*
  double xz, uz, tz, dz, sz, amplz, uz_uei, uz_ms, dwz, 
      hz, hz_tz, hz_dz, 
      mz, mz_uz, mz_ms, 
      rz, rz_uz, rz_ms;
      */
};
class boundary_layer {
    public:

    //TODO 副作用の除去
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