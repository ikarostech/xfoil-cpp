#include <Eigen/Core>

class ThomasAlgorithmResult {
public:
  Eigen::MatrixXd A;
  Eigen::VectorXd x;
};
class MathUtil {
public:
  static ThomasAlgorithmResult tridiagonalSolve(Eigen::MatrixXd A,
                                                Eigen::VectorXd d);
  static double pow(double a, int b);
  static double sign(double a, double b);
  static double atanc(double y, double x, double thold);
  static double cross2(const Eigen::Vector2d &a, const Eigen::Vector2d &b);
  static Eigen::Matrix2d getRotateMatrix(double alpha);
};
