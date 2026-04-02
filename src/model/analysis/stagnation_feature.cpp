#include "model/analysis/stagnation_feature.hpp"

StagnationFeature::StagnationFeature(int index, double sst, double sst_go, double sst_gp, bool found)
    : index_(index), sst_(sst), sst_go_(sst_go), sst_gp_(sst_gp), found_(found) {}
