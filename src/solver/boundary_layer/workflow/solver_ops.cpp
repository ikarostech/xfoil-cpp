#include "solver/boundary_layer/workflow/solver_ops.hpp"

bool BoundaryLayerSolverOps::blkin(BoundaryLayerStationWindow &state) const {
    return BoundaryLayerPhysics::blkin(state, context_.blCompressibility, context_.blReynolds);
}

SkinFrictionCoefficients BoundaryLayerSolverOps::blmid(FlowRegimeEnum flowRegimeType) const {
    return BoundaryLayerPhysics::blmid(context_.state, flowRegimeType);
}

blData BoundaryLayerSolverOps::blprv(blData data, double xsi, double ami, double cti, double thi, double dsi,
                                     double dswaki, double uei) const {
    return BoundaryLayerPhysics::blprv(data, context_.blCompressibility, xsi, ami, cti, thi, dsi, dswaki, uei);
}

bool BoundaryLayerSolverOps::blsys() const {
    blData &previous = context_.state.previous();
    blData &current  = context_.state.current();

    SkinFrictionCoefficients skinFriction = blmid(context_.flowRegime);
    current                               = context_.boundaryLayerVariablesSolver.solve(current, context_.flowRegime);

    if (context_.flowRegime == FlowRegimeEnum::Similarity) {
        context_.state.stepbl();
    }

    if (context_.flowRegime == FlowRegimeEnum::Transition) {
        context_.transitionSolver.trdif();
    } else {
        context_.blc = context_.blDiffSolver.solve(context_.flowRegime, context_.state, skinFriction,
                                                   context_.blTransition.amcrit);
    }

    if (context_.flowRegime == FlowRegimeEnum::Similarity) {
        context_.blc.a2 += context_.blc.a1;
        context_.blc.a1.setZero();
    }

    for (int k = 0; k < 4; ++k) {
        double res_u1 = context_.blc.a1(k, 3);
        double res_u2 = context_.blc.a2(k, 3);
        double res_ms = context_.blc.d_msq[k];

        context_.blc.a1(k, 3) *= previous.param.uz_uei;
        context_.blc.a2(k, 3) *= current.param.uz_uei;
        context_.blc.d_msq[k] = res_u1 * previous.param.uz_ms + res_u2 * current.param.uz_ms + res_ms;
    }

    return true;
}

bool BoundaryLayerSolverOps::tesys(const BoundaryLayerSideState &top_profiles,
                                   const BoundaryLayerSideState &bottom_profiles, const Edge &edge) const {
    context_.blc.clear();

    context_.state.station2 =
        context_.boundaryLayerVariablesSolver.solve(context_.state.station2, FlowRegimeEnum::Wake);

    const int top_te    = context_.lattice.top.trailingEdgeIndex;
    const int bottom_te = context_.lattice.bottom.trailingEdgeIndex;
    const double tte    = top_profiles.momentumThickness[top_te] + bottom_profiles.momentumThickness[bottom_te];
    const double dte =
        top_profiles.displacementThickness[top_te] + bottom_profiles.displacementThickness[bottom_te] + edge.ante;
    const double cte = (top_profiles.skinFrictionCoeff[top_te] * top_profiles.momentumThickness[top_te] +
                        bottom_profiles.skinFrictionCoeff[bottom_te] * bottom_profiles.momentumThickness[bottom_te]) /
                       tte;

    context_.blc.a1(0, 0) = -1.0;
    context_.blc.a2(0, 0) = 1.0;
    context_.blc.rhs[0]   = cte - context_.state.station2.param.sz;

    context_.blc.a1(1, 1) = -1.0;
    context_.blc.a2(1, 1) = 1.0;
    context_.blc.rhs[1]   = tte - context_.state.station2.param.tz;

    context_.blc.a1(2, 2) = -1.0;
    context_.blc.a2(2, 2) = 1.0;
    context_.blc.rhs[2]   = dte - context_.state.station2.param.dz - context_.state.station2.param.dwz;

    return true;
}
