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

    //TODO: INDEX_START_WITH=1
    int l = 1;
    int r = n;

    while (r - l > 1) {
        int imid = (int)((l + r) / 2);
        if (ss < s[imid]) {
            r = imid;
        }
        else {
            l = imid;
        }
    }
    int index = r;

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
 *       A simple averaging of adjacent segment slopes    |
 *      is used to achieve non-oscillatory curve.         |
 *      End conditions are set by end segment slope.      |
 *      To evaluate the spline at some value of s,        |
 *      use spline::seval and/or deval.                   |
 *                                                        |
 *      s        independent variable array (input)       |
 *      x        dependent variable array   (input)       |
 *      xs       dx/ds array                (calculated)  |
 *      n        number of points           (input)       |
 *                                                        |
 * -------------------------------------------------------*/
Eigen::VectorXd spline::splina(Eigen::VectorXd x, Eigen::VectorXd s, int n, int xs_size) {
    Eigen::VectorXd xs = Eigen::VectorXd::Zero(xs_size);
    double xs_former, xs_later;
    xs_former = xs_later = 0.0;

    bool is_duplicate = true;
    for (int i = spline::INDEX_START_WITH; i < n - 1 + spline::INDEX_START_WITH; i++) {
        const double ds = s[i + 1] - s[i];
        if (fabs(ds) < 1.e-10) {
            xs[i] = xs_former;
            is_duplicate = true;
            continue;
        } 
        
        const double dx = x[i + 1] - x[i];
        xs_later = dx / ds;
        if (is_duplicate) {
            xs[i] = xs_later;
            is_duplicate = false;
        } 
        else {
            xs[i] = (xs_former + xs_later) / 2;
        }
        xs_former = xs_later;
    }
    xs[n - 1 + INDEX_START_WITH] = xs_former;
    return xs;
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
bool spline::splind(double x[], double xs[], double s[], int n, double xs1, double xs2) {
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

  if (xs1 >= 998.0) {
    //----- set zero second derivative end condition
    matrixA(0, 0) = 2.0;
    matrixA(0, 1) = 1.0;
    vectorD(0) = 3.0 * (x[1 + INDEX_START_WITH] - x[0 + INDEX_START_WITH]) / (s[1 + INDEX_START_WITH] - s[0 + INDEX_START_WITH]);
  }
  else if (xs1 <= -998.0) {
      //----- set zero third derivative end condition
      matrixA(0, 0) = 1.0;
      matrixA(0, 1) = 1.0;
      vectorD(0) = 2.0 * (x[1 + INDEX_START_WITH] - x[0 + INDEX_START_WITH]) / (s[1 + INDEX_START_WITH] - s[0 + INDEX_START_WITH]);
  }
  else {
      //----- set specified first derivative end condition
      matrixA(0, 0) = 1.0;
      matrixA(0, 1) = 0.0;
      vectorD(0) = xs1;
  }

  if (xs2 >= 998.0) {
    matrixA(n - 1, n - 2) = 1.0;
    matrixA(n - 1, n - 1) = 2.0;
    vectorD(n - 1) = 3.0 * (x[n - 1 + INDEX_START_WITH] - x[n - 2 + INDEX_START_WITH]) / (s[n - 1 + INDEX_START_WITH] - s[n - 2 + INDEX_START_WITH]);
  }
  else if (xs2 <= -998.0) {
    matrixA(n - 1, n - 2) = 1.0;
    matrixA(n - 1, n - 1) = 1.0;
    vectorD(n - 1) = 2.0 * (x[n - 1 + INDEX_START_WITH] - x[n - 2 + INDEX_START_WITH]) / (s[n - 1  + INDEX_START_WITH] - s[n - 2 + INDEX_START_WITH]);
  } 
  else {
    matrixA(n - 1, n - 2) = 1.0;
    matrixA(n - 1, n - 1) = 0.0;
    vectorD(n - 1) = xs2;
  }
  

  if (n == 2 && xs1 <= -998.0 && xs2 <= -998.0) {
    matrixA(n - 1, n - 2) = 1.0;
    matrixA(n - 1, n - 1) = 2.0;
    vectorD(n - 1) = 3.0 * (x[n - 1 + INDEX_START_WITH] - x[n - 2 + INDEX_START_WITH]) / (s[n - 1  + INDEX_START_WITH] - s[n - 2 + INDEX_START_WITH]);
  }

  //---- solve for derivative array xs
  Eigen::VectorXd vectorXs = matrix::tridiagonalSolve(matrixA, vectorD).x;
  //FIXME xsの0とn移行が0でないと結果がおかしくなるバグが存在
  xs[0] = 0;
  xs[n] = 0;
  for (int i=0; i<n; i++) {
    xs[i + 1] = vectorXs(i);
  }
  return true;
}

/** -----------------------------------------------
 *      Splines x(s) array just like spline,      |
 *      but allows derivative discontinuities     |
 *      at segment joints.  Segment joints are    |
 *      defined by identical successive s values. |
 * ----------------------------------------------- */
Eigen::VectorXd spline::segspl(Eigen::VectorXd x, Eigen::VectorXd spline_length, int n) {
  return segspl(x, spline_length, n, -999.0, -999.0);
}

/** -----------------------------------------------
 *     splines x(s) array just like splind,      |
 *     but allows derivative discontinuities     |
 *     at segment joints.  segment joints are    |
 *     defined by identical successive s values. |
 * ----------------------------------------------- */
Eigen::VectorXd spline::segspl(Eigen::VectorXd x, Eigen::VectorXd spline_length, int n, double xs1, double xs2) {
  //TODO xsはxより1つ要素が多いのでそれに対応する
  Eigen::VectorXd xs = Eigen::VectorXd::Zero(x.size());

  int nseg, iseg, iseg0;

  iseg0 = 1;
  for (iseg = 2; iseg <= n - 2; iseg++) {
    if (spline_length[iseg] == spline_length[iseg + 1]) {
      nseg = iseg - iseg0 + 1;
      spline::splind(x.data() + iseg0 - 1, xs.data() + iseg0 - 1, spline_length.data() + iseg0 - 1, nseg, xs1, xs2);
      iseg0 = iseg + 1;
    }
  }
  nseg = n - iseg0 + 1;
  spline::splind(x.data() + iseg0 - 1, xs.data() + iseg0 - 1, spline_length.data() + iseg0 - 1, nseg, xs1, xs2);
  return xs;
}