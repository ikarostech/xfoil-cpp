#pragma once

#include "Eigen/Core"
#include "Eigen/Dense"

using namespace Eigen;

class Foil {
 public:
  int n;
  Matrix2Xd points;
  Matrix2Xd normal_vectors;
  Matrix2Xd dpoints_ds;
  VectorXd spline_length;
  Vector2d point_le;
  Vector2d point_te;
  double chord;
  double sle;
  Vector2d cmref;
  Matrix2Xd surface_vortex;
  VectorXd apanel;

  Foil()
      : n(0),
        chord(0.0),
        sle(0.0),
        cmref(Vector2d{0.25, 0.0}) {}
};

