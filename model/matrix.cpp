#include "matrix.hpp"

#include <iostream>
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

  int k;
  for (k = 2; k <= kk; k++) {
    int km = k - 1;
    c[km] = c[km] / a[km];
    d[km] = d[km] / a[km];
    a[k] = a[k] - b[k] * c[km];
    d[k] = d[k] - b[k] * d[km];
  }

  d[kk] = d[kk] / a[kk];

  for (k = kk - 1; k >= 1; k--) {
    d[k] = d[k] - c[k] * d[k + 1];
  }
  return true;
}