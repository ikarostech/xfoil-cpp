#pragma once
#include <Eigen/Core>
#include "spline.hpp"
class FoilShape {
    private:
    Eigen::VectorXd calcSpline() const;
    Eigen::Matrix2Xd calcDPointsDs() const;
    Eigen::Matrix2Xd calcNormalVector(const Eigen::Matrix2Xd &dpoints_ds) const;

    public:
    Eigen::Matrix2Xd points;  ///< Airfoil and wake points.
    int n = 0;                ///< Number of airfoil points.
    Eigen::VectorXd spline_length;  ///< Spline lengths for airfoil points.
    Eigen::Matrix2Xd dpoints_ds;    ///< First derivative of airfoil coordinates with respect to arc length.
    Eigen::Matrix2Xd normal_vector; ///< Normal vectors for airfoil points.

    FoilShape() = default;
    FoilShape(const Eigen::Matrix2Xd& points, int n);

    void setFoilShape(const Eigen::Matrix2Xd& points, int n);
};
