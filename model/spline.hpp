#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/StdVector>

class spline {
    public:
    const static int INDEX_START_WITH = 1;
    /**
     * @brief x,y座標から翼尾からの累計長さを求めます
     * 
     * @param x x座標のデータ
     * @param y y座標のデータ
     * @param n @deprecated データ数
     * @return std::vector<double> 始点からの累計長さのデータ
     */
    static Eigen::VectorXd scalc(const Eigen::Matrix2Xd points, int n);
    static double seval(double ss, Eigen::VectorXd x, Eigen::VectorXd xs, Eigen::VectorXd s, int n);

    /** --------------------------------------------------
    *	   calculates dx/ds(ss)                         |
    *	   xs array must have been calculated by spline |
    * -------------------------------------------------- */
    static double deval(double ss, Eigen::VectorXd x, Eigen::VectorXd xs, Eigen::VectorXd s, int n);

    /** --------------------------------------------------
     *      calculates d2x/ds2(ss)                          /
     *      xs array must have been calculated by spline    /
     * --------------------------------------------------- */
    static double d2val(double ss, Eigen::VectorXd x, Eigen::VectorXd xs, Eigen::VectorXd s, int n);

    static double sinvrt(double si, double xi, Eigen::VectorXd x, Eigen::VectorXd xs, Eigen::VectorXd s, int n);

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
     *                                                        |
     * ------------------------------------------------------- */
    static Eigen::VectorXd splind(Eigen::VectorXd x, Eigen::VectorXd s, int n);
};