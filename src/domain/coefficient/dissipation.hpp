#pragma once

#include "domain/flow_regime.hpp"
#include "domain/boundary_layer.hpp"

class dissipation {
 public:
  class DissipationResult {
  public:
    double di;
    double di_hk;
    double di_rt;
  };

  static DissipationResult getDissipation(double hk, double rt,
                                          FlowRegimeEnum flowRegimeType);
  static DissipationResult dil(double hk, double rt);
  static DissipationResult dilw(double hk, double rt);
};
