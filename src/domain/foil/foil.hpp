#pragma once

#include <Eigen/Core>

#include "core/spline.hpp"
#include "domain/foil/foil_shape.hpp"
#include "domain/foil/edge.hpp"

/**
 * @brief Holds the geometric definition of an airfoil.
 */
class Foil {
  private:
  public:
    FoilShape foil_shape;  // geometric domain (airfoil points and count)
    FoilShape wake_shape;  // geometric domain (wake points and count)

    Edge edge;

    Foil() = default;
    Foil(const Eigen::Matrix2Xd& points, int n) {
      foil_shape = FoilShape(points, n);
      edge = Edge(foil_shape);
      // xyWake();
    }
};
