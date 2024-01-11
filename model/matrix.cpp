#include "matrix.hpp"

#include "../Eigen/Core"
#include <iostream>
const int INDEX_START_WITH = 1;
//TODO trisolの除却
bool matrix::trisol(double a[], double b[], double c[], double d[], int kk) {
    //-----------------------------------------
    //     solves kk long, tri-diagonal system |
    //                                         |
    //             a c          d              |
    //             b a c        d              |
    //               b a .      .              |
    //                 . . c    .              |
    //                   b a    d              |
    //                                         |
    //     the righthand side d is replaced by |
    //     the solution.  a, c are destroyed.  |
    //-----------------------------------------

    for (int k = INDEX_START_WITH; k < kk; k++) {
        c[k] = c[k] / a[k];
        d[k] = d[k] / a[k];
        a[k + 1] = a[k + 1] - b[k + 1] * c[k];
        d[k + 1] = d[k + 1] - b[k + 1] * d[k];
    }

    d[kk] = d[kk] / a[kk];

    for (int k = kk - 1; k >= INDEX_START_WITH; k--) {
        d[k] = d[k] - c[k] * d[k + 1];
    }
    return true;
}

ThomasAlgorithmResult matrix::tridiagonalSolve(Eigen::MatrixXd A, Eigen::VectorXd d) {
    int n = A.rows();
    Eigen::VectorXd x(n);

    // Create temporary vectors for the algorithm
    Eigen::VectorXd a = A.diagonal(-1);
    Eigen::VectorXd b = A.diagonal();
    Eigen::VectorXd c = A.diagonal(1);
    
    // Forward elimination
    for (int i = 1; i < n; i++) {
        std::cout<<i<<std::endl;
        double factor = a(i - 1) / b(i - 1);
        b(i) -= factor * c(i - 1);
        d(i) -= factor * d(i - 1);
    }
    std::cout<<"forward"<<std::endl;
    // Backward substitution
    x(n - 1) = d(n - 1) / b(n - 1);

    for (int i = n - 2; i >= 0; i--) {
        x(i) = (d(i) - c(i) * x(i + 1)) / b(i);
    }
    
    ThomasAlgorithmResult result = ThomasAlgorithmResult();
    A.diagonal(-1) = a;
    A.diagonal() = b;
    A.diagonal(1) = c;
    result.A = A;
    result.d = d;
    result.x = x;
    return result;
}