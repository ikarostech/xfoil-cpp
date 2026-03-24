#include "model/boundary_layer/reference/bl_reynolds_params.hpp"

#include <cmath>

namespace {
constexpr double kHvrat = 0.35;
}

BlReynoldsParams::BlReynoldsParams(double currentRe, const BlCompressibilityParams &compressibility) {
    const double herat    = 1.0 - 0.5 * compressibility.qinfbl() * compressibility.qinfbl() * compressibility.hstinv();
    const double herat_ms = -0.5 * compressibility.qinfbl() * compressibility.qinfbl() * compressibility.hstinv_ms();
    const double factor   = std::sqrt(herat * herat * herat) * (1.0 + kHvrat) / (herat + kHvrat);
    reybl_                = currentRe * factor;
    reybl_re_             = factor;
    reybl_ms_             = reybl_ * (1.5 / herat - 1.0 / (herat + kHvrat)) * herat_ms;
}
