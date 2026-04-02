#include "model/analysis/stagnation_feature.hpp"

// TODO StagnationFeatureValuesを除却
namespace {
struct StagnationFeatureValues {
    int index     = 0;
    double sst    = 0.0;
    double sst_go = 0.0;
    double sst_gp = 0.0;
    bool found    = false;
};

StagnationFeatureValues computeStagnationFeatureValues(const Eigen::Matrix2Xd &surface_vortex,
                                                       const Eigen::VectorXd &spline_length) {
    StagnationFeatureValues values;
    const int point_count = static_cast<int>(surface_vortex.cols());
    for (int i = 0; i < point_count - 1; ++i) {
        if (surface_vortex(0, i) >= 0.0 && surface_vortex(0, i + 1) < 0.0) {
            values.index = i;
            values.found = true;
            break;
        }
    }

    if (!values.found) {
        values.index = point_count / 2;
    }

    const double dgam = surface_vortex(0, values.index + 1) - surface_vortex(0, values.index);
    const double ds   = spline_length[values.index + 1] - spline_length[values.index];

    if (surface_vortex(0, values.index) < -surface_vortex(0, values.index + 1)) {
        values.sst = spline_length[values.index] - ds * (surface_vortex(0, values.index) / dgam);
    } else {
        values.sst = spline_length[values.index + 1] - ds * (surface_vortex(0, values.index + 1) / dgam);
    }

    if (values.sst <= spline_length[values.index]) {
        values.sst = spline_length[values.index] + 0.0000001;
    }
    if (values.sst >= spline_length[values.index + 1]) {
        values.sst = spline_length[values.index + 1] - 0.0000001;
    }

    values.sst_go = (values.sst - spline_length[values.index + 1]) / dgam;
    values.sst_gp = (spline_length[values.index] - values.sst) / dgam;
    return values;
}
} // namespace

StagnationFeature::StagnationFeature(const Eigen::Matrix2Xd &surface_vortex, const Eigen::VectorXd &spline_length) {
    const auto values = computeStagnationFeatureValues(surface_vortex, spline_length);
    index_            = values.index;
    sst_              = values.sst;
    sst_go_           = values.sst_go;
    sst_gp_           = values.sst_gp;
    found_            = values.found;
}

StagnationFeature::StagnationFeature(int index, double sst, double sst_go, double sst_gp, bool found)
    : index_(index), sst_(sst), sst_go_(sst_go), sst_gp_(sst_gp), found_(found) {}
