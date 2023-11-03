#include <vector>

class spline {
    public:
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
};