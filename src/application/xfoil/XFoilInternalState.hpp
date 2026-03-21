#pragma once

#include <vector>

#include "Eigen/Core"
#include "model/foil/FoilAerodynamicCache.hpp"
#include "solver/boundary_layer/boundary_layer_geometry.hpp"

class XFoilInternalState {
  public:
    struct InitState {
        double amax = 0.0;
        std::vector<double> qf0;
        std::vector<double> qf1;
        std::vector<double> qf2;
        std::vector<double> qf3;
    };

    struct OperatingPointCouplingState {
        double machPerLift     = 0.0;
        double reynoldsPerLift = 0.0;
    };

    struct InviscidFieldState {
        Eigen::Matrix2Xd surfaceVortex;
        FoilAerodynamicCache cache;
        double gammaTe = 0.0;
        double sigmaTe = 0.0;
        Eigen::Matrix2Xd qinvu;
        Eigen::Matrix2Xd qinvMatrix;
        Eigen::VectorXd qgamm;
    };

    struct ViscousIterationState {
        StagnationResult stagnation;
        Eigen::VectorXd qvis;
        double convergedAlpha = 0.0;
        double convergedMach  = 0.0;
    };

    InitState init;
    OperatingPointCouplingState operatingPointCoupling;
    InviscidFieldState inviscid;
    ViscousIterationState viscous;
};
