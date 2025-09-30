#include "spline.hpp"
#include "math_util.hpp"
#include <iostream>

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

/**
 * 	   Calculates the "inverse" spline function s(x).
 * 	   Since s(x) can be multi-valued or not defined,
 * 	   this is not a "black-box" routine.  The calling
 * 	   program must pass via si a sufficiently good
 * 	   initial guess for s(xi).
 *
 * 	   xi	   specified x value	   (input)
 * 	   si	   calculated s(xi) value  (input,output)
 * 	   x,xs,s  usual spline arrays	   (input)
 */
double spline::sinvrt(double si, double xi, Eigen::VectorXd x, Eigen::VectorXd xs, Eigen::VectorXd spline_length, int n) {
  int iter;
  double sisav;
  sisav = si;

  for (iter = 1; iter <= 10; iter++) {
    const double res = spline::seval(si, x, xs, spline_length, n) - xi;
    const double resp = spline::deval(si, x, xs, spline_length, n);
    const double ds = -res / resp;
    si = si + ds;
    if (fabs(ds / (spline_length[n] - spline_length[1])) < 1.0e-5) return si;
  }

  si = sisav;

  return si;
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
 *                                                        |
 * ------------------------------------------------------- */
Eigen::VectorXd spline::splind(Eigen::VectorXd x, Eigen::VectorXd s) {
  Eigen::VectorXd xs = Eigen::VectorXd::Zero(x.size());
  int n = s.size();
  Eigen::MatrixXd matrixA = Eigen::MatrixXd(n, n);
  Eigen::VectorXd vectorD = Eigen::VectorXd(n);

  for (int i = 1; i < n - 1; i++) {
    const double dsm = s[i] - s[i - 1];
    const double dsp = s[i + 1] - s[i];
    matrixA(i, i - 1) = dsp;
    matrixA(i, i) = 2.0 * (dsm + dsp);
    matrixA(i, i + 1) = dsm;
    vectorD(i) =
        3.0 * ((x[i + 1] - x[i]) * dsm / dsp + (x[i] - x[i - 1]) * dsp / dsm);
  }

  //----- set zero third derivative end condition
  matrixA(0, 0) = 1.0;
  matrixA(0, 1) = 1.0;
  vectorD(0) = 2.0 * (x[1] - x[0]) / (s[1] - s[0]);

  matrixA(n - 1, n - 2) = 1.0;
  matrixA(n - 1, n - 1) = 1.0;
  vectorD(n - 1) = 2.0 * (x[n - 1] - x[n - 2]) / (s[n - 1] - s[n - 2]);

  //---- solve for derivative array xs
  Eigen::VectorXd vectorXs = MathUtil::tridiagonalSolve(matrixA, vectorD).x;
  return vectorXs;
}