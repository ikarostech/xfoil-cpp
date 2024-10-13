#include "math_util.hpp"

#include "../Eigen/Core"
#include <iostream>
const int INDEX_START_WITH = 1;

ThomasAlgorithmResult MathUtil::tridiagonalSolve(Eigen::MatrixXd A, Eigen::VectorXd d) {
    int n = A.rows();
    Eigen::VectorXd x(n);

    // Create temporary vectors for the algorithm
    Eigen::VectorXd a = A.diagonal(-1);
    Eigen::VectorXd b = A.diagonal();
    Eigen::VectorXd c = A.diagonal(1);
    
    // Forward elimination
    for (int i = 0; i < n - 1; i++) {
        c(i) /= b(i); 
        d(i) /= b(i);

        b(i + 1) -= a(i) * c(i);
        d(i + 1) -= a(i) * d(i);

    }

    // Backward substitution
    d(n - 1) = d(n - 1) / b(n - 1);

    for (int i = n - 2; i >= 0; i--) {
        d(i) -= c(i) * d(i + 1);
    }
    
    ThomasAlgorithmResult result = ThomasAlgorithmResult();
    A.diagonal(-1) = a;
    A.diagonal() = b;
    A.diagonal(1) = c;
    result.A = A;
    result.x = d;
    return result;
}

double MathUtil::pow(double a, int b) {
    if (b < 1) {
        return 1;
    }
    else {
        return a * pow(a, b - 1);
    }
}