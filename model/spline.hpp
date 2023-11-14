#include <vector>

class spline {
    public:
    const static int INDEX_START_WITH = 1;
    /*
    /**
     * @brief x,y座標から翼尾からの累計長さを求めます
     * 
     * @param x x座標のデータ
     * @param y y座標のデータ
     * @param n @deprecated データ数
     * @return std::vector<double> 始点からの累計長さのデータ
     */
    static std::vector<double> scalc(const double x[], const double y[], int n, const int s_size);
    static double seval(double ss, const double x[], const double xs[], const double s[], int n);

    /** --------------------------------------------------
    *	   calculates dx/ds(ss)                         |
    *	   xs array must have been calculated by spline |
    * -------------------------------------------------- */
    static double deval(double ss, const double x[], const double xs[], const double s[], int n);

    /** -------------------------------------------------------
     *      Calculates spline coefficients for x(s).          |
     *       A simple averaging of adjacent segment slopes    |
     *      is used to achieve non-oscillatory curve.         |
     *      End conditions are set by end segment slope.      |
     *      To evaluate the spline at some value of s,        |
     *      use spline::seval and/or deval.                           |
     *                                                        |
     *      s        independent variable array (input)       |
     *      x        dependent variable array   (input)       |
     *      xs       dx/ds array                (calculated)  |
     *      n        number of points           (input)       |
     *                                                        |
     * -------------------------------------------------------*/
    static std::vector<double> splina(const double x[], const double s[], int n, int xs_size);
};