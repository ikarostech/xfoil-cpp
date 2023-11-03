#include "spline.hpp"
#include "../Eigen/Core"
#include "../Eigen/Dense"
#include "../Eigen/StdVector"
#include <iostream>

const int INDEX_START_WITH = 1;

std::vector<double> spline::scalc(const double x[], const double y[], int n, const int s_size) {
    //TODO 引数のVector2d化
    std::vector<Eigen::Vector2d> pos(n + INDEX_START_WITH);
    std::cout<<n<<std::endl;
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