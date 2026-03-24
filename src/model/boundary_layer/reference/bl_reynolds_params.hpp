#pragma once
#include "model/boundary_layer/reference/bl_compressibility_params.hpp"

class BlReynoldsParams {
  public:
    BlReynoldsParams() = default;
    BlReynoldsParams(double currentRe, const BlCompressibilityParams &compressibility);

    double reybl() const {
        return reybl_;
    }
    double reybl_re() const {
        return reybl_re_;
    }
    double reybl_ms() const {
        return reybl_ms_;
    }

  private:
    double reybl_    = 0.0;
    double reybl_re_ = 0.0;
    double reybl_ms_ = 0.0;
};
