#pragma once

#include "model/boundary_layer/physics.hpp"
#include "solver/boundary_layer/workflow/workflow.hpp"

class BoundaryLayerSolverOps {
  public:
    struct Context {
        BoundaryLayerVariablesSolver &boundaryLayerVariablesSolver;
        BlDiffSolver &blDiffSolver;
        BoundaryLayerTransitionSolver &transitionSolver;
        FlowRegimeEnum &flowRegime;
        const BlCompressibilityParams &blCompressibility;
        const BlReynoldsParams &blReynolds;
        const BlTransitionParams &blTransition;
        BoundaryLayerStationWindow &state;
        BlSystemCoeffs &blc;
        const SidePair<BoundaryLayerLattice> &lattice;
    };

    explicit BoundaryLayerSolverOps(Context context) : context_(context) {}

    bool blkin(BoundaryLayerStationWindow &state) const;
    SkinFrictionCoefficients blmid(FlowRegimeEnum flowRegimeType) const;
    blData blprv(blData data, double xsi, double ami, double cti, double thi, double dsi, double dswaki,
                 double uei) const;
    bool blsys() const;
    bool tesys(const BoundaryLayerSideState &top_profiles, const BoundaryLayerSideState &bottom_profiles,
               const Edge &edge) const;

  private:
    Context context_;
};
