#pragma once

class StagnationFeature {
  public:
    StagnationFeature() = default;
    StagnationFeature(int index, double sst, double sst_go, double sst_gp, bool found);

    int index() const {
        return index_;
    }
    double sst() const {
        return sst_;
    }
    double sst_go() const {
        return sst_go_;
    }
    double sst_gp() const {
        return sst_gp_;
    }
    bool found() const {
        return found_;
    }

  private:
    int index_      = 0;
    double sst_     = 0.0;
    double sst_go_  = 0.0;
    double sst_gp_  = 0.0;
    bool found_     = true;
};
