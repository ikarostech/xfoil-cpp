#include "../../enum/flowRegimeEnum.hpp"

class skin_friction {
public:
   class C_f {
    public:
    double cf = 0.0;
    double hk = 0.0;
    double rt = 0.0;
    /** squared freestream mach number at current cl*/
    double msq = 0.0;
  };
  static C_f getSkinFriction(double hk, double rt, double msq, FlowRegimeEnum flowRegimeType);
private:
  static C_f cfl(double hk, double rt);
  static C_f cft(double hk, double rt, double msq);
};