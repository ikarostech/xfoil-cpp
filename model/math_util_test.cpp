#include <gtest/gtest.h>
#include <vector>

#include "math_util.hpp"

TEST(tridiagonalSolve, test_value) {
    //given
    Eigen::MatrixXd A(3, 3);
    A <<  3, 1, 0,
          1, 4, 2,
          0, 2, 5;
    Eigen::VectorXd d(3);
    d << 5, 15, 19;

    Eigen::VectorXd expectedX(3);
    expectedX << 1, 2, 3;

    //when
    ThomasAlgorithmResult actual = MathUtil::tridiagonalSolve(A, d);
    std::cout<<actual.A<<std::endl;

    //then
    ASSERT_FLOAT_EQ(actual.x(0), expectedX(0));
    ASSERT_FLOAT_EQ(actual.x(1), expectedX(1));
    ASSERT_FLOAT_EQ(actual.x(2), expectedX(2));
}
TEST(pow, test_value) {
    //when
    double actual = MathUtil::pow(3, 2);

    //then
    ASSERT_FLOAT_EQ(actual, 9);
}
int main() {
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}