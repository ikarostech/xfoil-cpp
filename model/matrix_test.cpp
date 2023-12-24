#include <gtest/gtest.h>
#include <vector>

#include "matrix.hpp"

TEST(trisol, test_value) {
    //given
    std::vector<double> a{0.0, 2, 7, -3};
    std::vector<double> b{0.0, -1, 2};
    std::vector<double> c{0.0, 1, 4};
    std::vector<double> d{0.0, 4, 25, -5};

    //when
    matrix::trisol(a.data(), b.data(), c.data(), d.data(), 3);

    //then
    ASSERT_EQ(0, d[0]);
    ASSERT_EQ(0.80555555555555558, d[1]);
    ASSERT_EQ(2.3888888888888888, d[2]);
    ASSERT_EQ(1.6666666666666667, d[3]);

    ASSERT_EQ(0, a[0]);
    ASSERT_EQ(2, a[1]);
    ASSERT_EQ(6, a[2]);
    ASSERT_EQ(-3, a[3]);

    ASSERT_EQ(0, b[0]);
    ASSERT_EQ(-1, b[1]);
    ASSERT_EQ(2, b[2]);
}

int main() {
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}