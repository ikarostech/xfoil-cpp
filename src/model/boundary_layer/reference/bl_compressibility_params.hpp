#pragma once

class BlCompressibilityParams {
  public:
    BlCompressibilityParams() = default;
    BlCompressibilityParams(double qinf, double currentMach);

    double gm1bl() const { return gm1bl_; }
    double qinfbl() const { return qinfbl_; }
    double tkbl() const { return tkbl_; }
    double tkbl_ms() const { return tkbl_ms_; }
    double rstbl() const { return rstbl_; }
    double rstbl_ms() const { return rstbl_ms_; }
    double hstinv() const { return hstinv_; }
    double hstinv_ms() const { return hstinv_ms_; }

  private:
    double gm1bl_     = 0.0;
    double qinfbl_    = 0.0;
    double tkbl_      = 0.0;
    double tkbl_ms_   = 0.0;
    double rstbl_     = 0.0;
    double rstbl_ms_  = 0.0;
    double hstinv_    = 0.0;
    double hstinv_ms_ = 0.0;
};
