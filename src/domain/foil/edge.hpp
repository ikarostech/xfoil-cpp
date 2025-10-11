#pragma once

#include <Eigen/Core>

class FoilShape;

class Edge {
  public:
    Edge();
    explicit Edge(const FoilShape& foilShape);

    Eigen::Vector2d point_le;
    Eigen::Vector2d point_te;
    double chord;
    double sle;

  private:
    static double lefind(const Eigen::Matrix2Xd& points,
                         const Eigen::Matrix2Xd& dpoints_ds,
                         const Eigen::VectorXd& s,
                         int n);
};
