#include <gtest/gtest.h>

#include "dissipation.hpp"

TEST(DissipationTest, DilBelow4) {
    double hk = 3.0;
    double rt = 1e5;
    auto res = dissipation::dil(hk, rt);
    EXPECT_NEAR(2.0904999999999997e-06, res.di, 1e-12);
    EXPECT_NEAR(-1.1275e-07, res.di_hk, 1e-12);
    EXPECT_NEAR(-2.0904999999999998e-11, res.di_rt, 1e-17);
}

TEST(DissipationTest, DilAbove4) {
    double hk = 5.0;
    double rt = 1e5;
    auto res = dissipation::dil(hk, rt);
    EXPECT_NEAR(2.0543137254901963e-06, res.di, 1e-12);
    EXPECT_NEAR(-3.0757400999615534e-08, res.di_hk, 1e-12);
    EXPECT_NEAR(-2.0543137254901964e-11, res.di_rt, 1e-17);
}

TEST(DissipationTest, DilwWake) {
    double hk = 5.0;
    double rt = 1e5;
    auto res = dissipation::dilw(hk, rt);
    EXPECT_NEAR(1.841404463247928e-06, res.di, 1e-12);
    EXPECT_NEAR(-5.568121217349049e-07, res.di_hk, 1e-12);
    EXPECT_NEAR(-1.8414044632479282e-11, res.di_rt, 1e-17);
}

TEST(DissipationTest, GetDissipationFlowRegime) {
    double hk = 5.0;
    double rt = 1e5;
    auto lam = dissipation::getDissipation(hk, rt, XFoil::FlowRegimeEnum::Laminar);
    auto lamExp = dissipation::dil(hk, rt);
    EXPECT_DOUBLE_EQ(lamExp.di, lam.di);
    EXPECT_DOUBLE_EQ(lamExp.di_hk, lam.di_hk);
    EXPECT_DOUBLE_EQ(lamExp.di_rt, lam.di_rt);

    auto wake = dissipation::getDissipation(hk, rt, XFoil::FlowRegimeEnum::Wake);
    auto wakeExp = dissipation::dilw(hk, rt);
    EXPECT_DOUBLE_EQ(wakeExp.di, wake.di);
    EXPECT_DOUBLE_EQ(wakeExp.di_hk, wake.di_hk);
    EXPECT_DOUBLE_EQ(wakeExp.di_rt, wake.di_rt);
}

int main() {
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}

