#include "spline.hpp"

#include <iostream>


/**	  Calculates x(ss)
 *	   xs array must have been calculated by spline */
double spline::seval(double ss, const double x[], const double xs[], const double s[], int n) {
    int ilow, i;
    double ds, t, cx1, cx2;

    ilow = 1;
    i = n;

    while (i - ilow > 1) {
        int imid = (int)((i + ilow) / 2);
        if (ss < s[imid])
        i = imid;
        else
        ilow = imid;
    }

    ds = s[i] - s[i - 1];
    t = (ss - s[i - 1]) / ds;
    cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    cx2 = ds * xs[i] - x[i] + x[i - 1];
    std::cout<<"ss:"<<ss<<std::endl;
    std::cout<<"x:";
    for (int k=1; k<=n; k++) {
        std::cout<<x[k]<<",";
    }
    std::cout<<std::endl;
    std::cout<<"s:";
    for (int k=1; k<=n; k++) {
        std::cout<<s[k]<<",";
    }
    std::cout<<std::endl;
    std::cout<<"xs:";
    for (int k=1; k<=n; k++) {
        std::cout<<xs[k]<<",";
    }
    std::cout<<std::endl;
    double seval = t * x[i] + (1.0 - t) * x[i - 1] +
            (t - t * t) * ((1.0 - t) * cx1 - t * cx2);
    std::cout<<"seval:"<<seval<<std::endl;
    return seval;
}