#include "model/boundary_layer/reference/bl_compressibility_params.hpp"

#include <cmath>

#include "numerics/math_util.hpp"

BlCompressibilityParams::BlCompressibilityParams(double qinf,
                                                 double currentMach) {
    const double beta                = std::sqrt(1.0 - currentMach * currentMach);
    const double karman_tsien_factor = std::pow(currentMach / (1.0 + beta), 2);

    gm1bl_     = 0.4;
    qinfbl_    = qinf;
    tkbl_      = karman_tsien_factor;
    tkbl_ms_   = 1.0 / ((1.0 + beta) * (1.0 + beta)) - 2.0 * karman_tsien_factor / (1.0 + beta) * (-0.5 / beta);
    rstbl_     = std::pow(1.0 + 0.5 * gm1bl_ * currentMach * currentMach, 1.0 / gm1bl_);
    rstbl_ms_  = 0.5 * rstbl_ / (1.0 + 0.5 * gm1bl_ * currentMach * currentMach);
    hstinv_    = gm1bl_ * MathUtil::pow(currentMach / qinf, 2) / (1.0 + 0.5 * gm1bl_ * currentMach * currentMach);
    hstinv_ms_ = gm1bl_ * MathUtil::pow(1.0 / qinf, 2) / (1.0 + 0.5 * gm1bl_ * currentMach * currentMach) -
                 0.5 * gm1bl_ * hstinv_ / (1.0 + 0.5 * gm1bl_ * currentMach * currentMach);
}
