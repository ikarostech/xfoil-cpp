#include "matrix.hpp"

#include <iostream>
const int INDEX_START_WITH = 1;
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