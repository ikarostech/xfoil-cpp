#include "../../XFoil.h"

class dissipation {
public:
  class DissipationResult {
  public:
    double di;
    double di_hk;
    double di_rt;
  };
  static DissipationResult dil(double hk, double rt);
  static DissipationResult dilw(double hk, double rt);
  static DissipationResult getDissipation(double hk, double rt,
                                          XFoil::FlowRegimeEnum flowRegimeType);
};

