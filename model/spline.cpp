#include "spline.hpp"
#include "matrix.hpp"
#include <iostream>
Eigen::VectorXd spline::scalc(Eigen::MatrixX2d points, int n, const int s_size) {
    Eigen::VectorXd s = Eigen::VectorXd::Zero(s_size);
    
    for (int i = 1; i < n; i++) {
        s[i] = s[i - 1] + (points.row(i) - points.row(i - 1)).norm();
    }

    return s;
}
/**	  Calculates x(ss)
 *	   xs array must have been calculated by spline */
double spline::seval(double ss, Eigen::VectorXd x, Eigen::VectorXd xs, Eigen::VectorXd s, int n) {

    int index = std::distance(s.begin(), std::lower_bound(s.begin(), s.begin() + n, ss));

    double ds = s[index] - s[index - 1];
    double t = (ss - s[index - 1]) / ds;
    double cx_former = ds * xs[index - 1] - x[index] + x[index - 1];
    double cx_later = ds * xs[index] - x[index] + x[index - 1];
    return t * x[index] + (1.0 - t) * x[index - 1] +
            (t - t * t) * ((1.0 - t) * cx_former - t * cx_later);
}

/** --------------------------------------------------
*	   calculates dx/ds(ss)                         |
*	   xs array must have been calculated by spline |
* -------------------------------------------------- */
double spline::deval(double ss, Eigen::VectorXd x, Eigen::VectorXd xs, Eigen::VectorXd s, int n) {

    int i = std::distance(s.begin(), std::lower_bound(s.begin(), s.begin() + n, ss));
    
    const double ds = s[i] - s[i - 1];
    const double dx = x[i] - x[i - 1];
    const double t = (ss - s[i - 1]) / ds;
    const double cx1 = ds * xs[i - 1] - dx;
    const double cx2 = ds * xs[i] - dx;
    return (dx + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 +
            t * (3.0 * t - 2.0) * cx2) / ds;
}

/** --------------------------------------------------
 *      calculates d2x/ds2(ss)                          /
 *      xs array must have been calculated by spline    /
 * --------------------------------------------------- */
double spline::d2val(double ss, Eigen::VectorXd x, Eigen::VectorXd xs, Eigen::VectorXd s, int n) {

    int i = std::distance(s.begin(), std::lower_bound(s.begin(), s.begin() + n, ss));
    
    const double ds = s[i] - s[i - 1];
    const double t = (ss - s[i - 1]) / ds;
    const double cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    const double cx2 = ds * xs[i] - x[i] + x[i - 1];
    return ((6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2) / pow(ds, 2.0);
}

/** -------------------------------------------------------
 *      Calculates spline coefficients for x(s).          |
 *      Specified 1st derivative and/or usual zero 2nd    |
 *      derivative end conditions are used.               |
 *                                                        |
 *      To evaluate the spline at some value of s,        |
 *      use spline::seval and/or deval.                   |
 *                                                        |
 *      s        independent variable array (input)       |
 *      x        dependent variable array   (input)       |
 *      xs       dx/ds array                (calculated)  |
 *      n        number of points           (input)       |
 *      xs1,xs2  endpoint derivatives       (input)       |
 *               if = 999.0, then usual zero second       |
 *               derivative end condition(s) are used     |
 *               if = -999.0, then zero third             |
 *               derivative end condition(s) are used     |
 *                                                        |
 * ------------------------------------------------------- */
Eigen::VectorXd spline::splind(Eigen::VectorXd x, Eigen::VectorXd s, int n) {
  Eigen::VectorXd xs = Eigen::VectorXd::Zero(x.size());
  Eigen::MatrixXd matrixA = Eigen::MatrixXd(n, n);
  Eigen::VectorXd vectorD = Eigen::VectorXd(n);

  for (int i = 1; i < n - 1; i++) {
    const double dsm = s[i + INDEX_START_WITH] - s[i - 1 + INDEX_START_WITH];
    const double dsp = s[i + 1 + INDEX_START_WITH] - s[i + INDEX_START_WITH];
    matrixA(i, i - 1) = dsp;
    matrixA(i, i) = 2.0 * (dsm + dsp);
    matrixA(i, i + 1) = dsm;
    vectorD(i) =
        3.0 * ((x[i + 1 + INDEX_START_WITH] - x[i + INDEX_START_WITH]) * dsm / dsp + (x[i + INDEX_START_WITH] - x[i - 1 + INDEX_START_WITH]) * dsp / dsm);
  }

  //----- set zero third derivative end condition
  matrixA(0, 0) = 1.0;
  matrixA(0, 1) = 1.0;
  vectorD(0) = 2.0 * (x[1 + INDEX_START_WITH] - x[0 + INDEX_START_WITH]) / (s[1 + INDEX_START_WITH] - s[0 + INDEX_START_WITH]);

  matrixA(n - 1, n - 2) = 1.0;
  matrixA(n - 1, n - 1) = 1.0;
  vectorD(n - 1) = 2.0 * (x[n - 1 + INDEX_START_WITH] - x[n - 2 + INDEX_START_WITH]) / (s[n - 1  + INDEX_START_WITH] - s[n - 2 + INDEX_START_WITH]);

  //---- solve for derivative array xs
  Eigen::VectorXd vectorXs = matrix::tridiagonalSolve(matrixA, vectorD).x;
  //FIXME xsの0とn移行が0でないと結果がおかしくなるバグが存在
  xs[0] = 0;
  xs[n] = 0;
  for (int i=0; i<n; i++) {
    xs[i + 1] = vectorXs(i);
  }
  return xs;
}