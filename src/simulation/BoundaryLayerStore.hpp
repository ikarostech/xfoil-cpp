#pragma once
#include <array>
#include <boundary_layer.hpp>

class BoundaryLayerStore {
 public:
  // Saves the current boundary layer station data snapshot for later restore.
  void saveblData(blData blData, int icom);

  // Restores boundary layer station data from the previously saved snapshot.
  blData restoreblData(int icom);
  private:
  std::array<blData, 3> blsav{};
};
