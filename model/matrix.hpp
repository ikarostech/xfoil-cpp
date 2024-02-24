#include "../Eigen/Core"

class ThomasAlgorithmResult {
    public:
    Eigen::MatrixXd A;
    Eigen::VectorXd x;
};
class matrix {
    public:
    static ThomasAlgorithmResult tridiagonalSolve(Eigen::MatrixXd A, Eigen::VectorXd d);
};