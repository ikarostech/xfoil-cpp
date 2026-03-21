#pragma once

#include "application/xfoil/XFoilInternalState.hpp"
#include "model/foil/foil.hpp"
#include "solver/boundary_layer/boundary_layer.hpp"

struct XFoilWorkspace {
  XFoilInternalState state;
  Foil foil;
  BoundaryLayer boundaryLayer;
};
