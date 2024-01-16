#include "matrix.hpp"

#include "../Eigen/Core"
#include <iostream>
const int INDEX_START_WITH = 1;
//TODO trisolの除却
bool matrix::trisol(double a[], double b[], double c[], double d[], int kk) {
    //-----------------------------------------
    //     solves kk long, tri-diagonal system |
    //                                         |
    //             a1 c1          d1           |
    //             b2 a2 c2       d2           |
    //               b3 a3 .       .           |
    //                 . . cn-1    .           |
    //                   bn an    dn           |
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