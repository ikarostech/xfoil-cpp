#include "../Eigen/Core"

class blData {
    public:
  class blVector {
  public:
    double scalar;
    Eigen::Vector<double, 6> vector = Eigen::Vector<double, 6>::Zero();
    Eigen::Vector<double, 3> pos_vector() {
      return vector.segment(0, 3);
    }
    double& u() {
      return vector[0];
    }
    double& t() {
      return vector[1];
    }
    double& d() {
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
  double xz, uz, tz, dz, sz, amplz, uz_uei, uz_ms, dwz, 
      hz, hz_tz, hz_dz, 
      mz, mz_uz, mz_ms, 
      rz, rz_uz, rz_ms;
};
class boundary_layer {
    public:
    static bool hct(double hk, double msq, double &hc, double &hc_hk, double &hc_msq);
    static bool hkin(double h, double msq, double &hk, double &hk_h, double &hk_msq);
    static bool hsl(double hk, double &hs, double &hs_hk, double &hs_rt, double &hs_msq);
    static bool hst(double hk, double rt, double msq, double &hs, double &hs_hk, double &hs_rt, double &hs_msq);
};