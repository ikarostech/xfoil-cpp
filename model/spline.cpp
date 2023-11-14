#include "spline.hpp"
#include "../Eigen/Core"
#include "../Eigen/Dense"
#include "../Eigen/StdVector"

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
    int ilow, i;
    double ds, t, cx1, cx2, deval;

    ilow = 1;
    //	i = nc;
    i = n;  /// techwinder modified

    while (i - ilow > 1) {
        int imid = (i + ilow) / 2;
        if (ss < s[imid])
        i = imid;
        else
        ilow = imid;
    }

    ds = s[i] - s[i - 1];
    t = (ss - s[i - 1]) / ds;
    cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    cx2 = ds * xs[i] - x[i] + x[i - 1];
    deval = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 +
            t * (3.0 * t - 2.0) * cx2;
    deval = deval / ds;
    return deval;
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