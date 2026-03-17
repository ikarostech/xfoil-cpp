#include "solver/boundary_layer/boundary_layer_store.hpp"
std::array<blData, 3> blsav{};
void BoundaryLayerStore::saveblData(blData blData, int icom) {
  blsav[icom] = blData;
}

blData BoundaryLayerStore::restoreblData(int icom) { return blsav[icom]; }
