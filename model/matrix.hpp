#include "../Eigen/Core"

class ThomasAlgorithmResult {
    public:
    Eigen::MatrixXd A;
    Eigen::VectorXd d;
    Eigen::VectorXd x;
};
class matrix {
    public:
    static bool trisol(double a[], double b[], double c[], double d[], int kk);
    static ThomasAlgorithmResult tridiagonalSolve(Eigen::MatrixXd A, Eigen::VectorXd d);
};