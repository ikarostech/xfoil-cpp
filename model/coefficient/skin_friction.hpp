#pragma once

class skin_friction {
public:
   class C_f {
    public:
    double cf;
    double hk;
    double rt;
    /** squared freestream mach number at current cl*/
    double msq;
  };
  static C_f cfl(double hk, double rt);
  static C_f cft(double hk, double rt, double msq);
};