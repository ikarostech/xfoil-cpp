#include "spline.hpp"
#include "matrix.hpp"
#include "../Eigen/Core"
#include "../Eigen/Dense"
#include "../Eigen/StdVector"
#include <iostream>
std::vector<double> spline::scalc(const double x[], const double y[], int n, const int s_size) {
    //TODO 引数のVector2d化
    std::vector<Eigen::Vector2d> pos(n + INDEX_START_WITH);
    for (int i=0; i<n+INDEX_START_WITH; i++) {
        pos[i][0] = x[i];
        pos[i][1] = y[i];
        
    }
    std::vector<double> s(s_size, 0);
    
    for (int i = 1 + INDEX_START_WITH; i < n + INDEX_START_WITH; i++) {
        s[i] = s[i - 1] + (pos[i] - pos[i - 1]).norm();
    }

    return s;
}
/**	  Calculates x(ss)
 *	   xs array must have been calculated by spline */
double spline::seval(double ss, const double x[], const double xs[], const double s[], int n) {

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
double spline::deval(double ss, const double x[], const double xs[], const double s[], int n) {

    //FIXME 引数のVector化
    std::vector<double> s_vector(s, s + n - 1 + INDEX_START_WITH);
    int i = std::distance(s_vector.begin(), std::lower_bound(s_vector.begin(), s_vector.end(), ss));

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
double spline::d2val(double ss, const double x[], const double xs[], const double s[], int n) {
    //FIXME 引数のVector化
    std::vector<double> s_vector(s, s + n - 1 + INDEX_START_WITH);
    int i = std::distance(s_vector.begin(), std::lower_bound(s_vector.begin(), s_vector.end(), ss));

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
std::vector<double> spline::splina(const double x[], const double s[], int n, int xs_size) {
    std::vector<double> xs(xs_size, 0);
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
  int nmax = 600;
  double a[n + 1], b[n + 1], c[n + 1];

  for (int i = 2; i <= n - 1; i++) {
    const double dsm = s[i] - s[i - 1];
    const double dsp = s[i + 1] - s[i];
    b[i] = dsp;
    a[i] = 2.0 * (dsm + dsp);
    c[i] = dsm;
    xs[i] =
        3.0 * ((x[i + 1] - x[i]) * dsm / dsp + (x[i] - x[i - 1]) * dsp / dsm);
  }

  if (xs1 >= 998.0) {
    //----- set zero second derivative end condition
    a[1] = 2.0;
    c[1] = 1.0;
    xs[1] = 3.0 * (x[2] - x[1]) / (s[2] - s[1]);
  } else {
    if (xs1 <= -998.0) {
      //----- set zero third derivative end condition
      a[1] = 1.0;
      c[1] = 1.0;
      xs[1] = 2.0 * (x[2] - x[1]) / (s[2] - s[1]);
    } else {
      //----- set specified first derivative end condition
      a[1] = 1.0;
      c[1] = 0.0;
      xs[1] = xs1;
    }
  }

  if (xs2 >= 998.0) {
    b[n] = 1.0;
    a[n] = 2.0;
    xs[n] = 3.0 * (x[n] - x[n - 1]) / (s[n] - s[n - 1]);
  } else {
    if (xs2 <= -998.0) {
      b[n] = 1.0;
      a[n] = 1.0;
      xs[n] = 2.0 * (x[n] - x[n - 1]) / (s[n] - s[n - 1]);
    } else {
      a[n] = 1.0;
      b[n] = 0.0;
      xs[n] = xs2;
    }
  }

  if (n == 2 && xs1 <= -998.0 && xs2 <= -998.0) {
    b[n] = 1.0;
    a[n] = 2.0;
    xs[n] = 3.0 * (x[n] - x[n - 1]) / (s[n] - s[n - 1]);
  }

  //---- solve for derivative array xs
  Eigen::MatrixXd matrixA = Eigen::MatrixXd(n, n);
  Eigen::VectorXd vectorD = Eigen::VectorXd(n);
  for (int i=0; i<n; i++) {
    matrixA(i,i) = a[i + 1];
    vectorD(i) = xs[i + 1];
  }
  for (int i=0; i<n-1; i++) {
    matrixA(i, i + 1) = c[i + 1];
    matrixA(i + 1, i) = b[i + 2];
  }
  Eigen::VectorXd vectorXs = matrix::tridiagonalSolve(matrixA, vectorD).x;
  //FIXME xsの0とn移行が0でないと結果がおかしくなるバグが存在
  for (int i=0; i<n; i++) {
    xs[i + 1] = vectorXs(i);
  }
  return true;
}