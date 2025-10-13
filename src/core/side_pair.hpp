#pragma once

#include <stdexcept>

/**
 * @brief Convenience container for holding values on the airfoil's
 *        top and bottom surfaces.
 *
 * Historically this template lived inside the XFoil class.  Moving it
 * to a standalone header allows the boundary-layer infrastructure to
 * reference the same abstraction without depending on XFoil directly.
 */
template <class T>
struct SidePair {
  T top;
  T bottom;

  T& get(int side) {
    switch (side) {
      case 1:
        return top;
      case 2:
        return bottom;
      default:
        throw std::invalid_argument("invalid side type");
    }
  }

  const T& get(int side) const {
    switch (side) {
      case 1:
        return top;
      case 2:
        return bottom;
      default:
        throw std::invalid_argument("invalid side type");
    }
  }
};
