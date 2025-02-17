#include <Eigen/Core>

class ThomasAlgorithmResult {
    public:
    Eigen::MatrixXd A;
    Eigen::VectorXd x;
};
class MathUtil {
    public:
    static ThomasAlgorithmResult tridiagonalSolve(Eigen::MatrixXd A, Eigen::VectorXd d);
    static double pow(double a, int b);
};