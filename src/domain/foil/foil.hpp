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
      wake_shape.n = foil_shape.n / 8 + 2;
      edge = Edge(foil_shape);
    }
    const FoilShape& getFoilShape(int index) const {
      if (0 <= index && index < foil_shape.n) {
        return foil_shape;
      } else if (foil_shape.n <= index &&
                 index < foil_shape.n + wake_shape.n) {
        return wake_shape;
      }
      throw std::out_of_range("Foil::getFoilShape: index out of range");
    }

    bool xyWake(int wake_point_count,
                const Eigen::Matrix2Xd &gamu,
                const Eigen::Matrix2Xd &surface_vortex,
                double alfa, double qinf);
};
