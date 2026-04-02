#include "solver/xfoil/xfoil_flowfield.hpp"

#include <algorithm>
#include <cmath>

#include "Eigen/Dense"
#include "infrastructure/logger.hpp"
#include "numerics/math_util.hpp"
#include "solver/inviscid/InviscidSolver.hpp"
#include "solver/inviscid/psi.hpp"

namespace xfoil_flowfield {

void updateTrailingEdgeState(const Foil &foil, XFoilInternalState::InviscidFieldState &state) {
    const int node_count = foil.foil_shape.n;

    if (node_count < 2) {
        state.sigmaTe = 0.0;
        state.gammaTe = 0.0;
        return;
    }

    double scs = 0.0;
    double sds = 0.0;
    if (foil.edge.sharp) {
        scs = 1.0;
        sds = 0.0;
    } else if (foil.edge.dste != 0.0) {
        const double inv_dste = 1.0 / foil.edge.dste;
        scs                   = foil.edge.ante * inv_dste;
        sds                   = foil.edge.aste * inv_dste;
    }

    if (state.surfaceVortex.rows() > 0 && state.surfaceVortex.cols() >= node_count) {
        const double surface_delta = state.surfaceVortex(0, 0) - state.surfaceVortex(0, node_count - 1);
        state.sigmaTe              = 0.5 * surface_delta * scs;
        state.gammaTe              = -0.5 * surface_delta * sds;
        return;
    }

    state.sigmaTe = 0.0;
    state.gammaTe = 0.0;
}

Eigen::Matrix2Xd buildSurfaceVortex(const Foil &foil, const XFoilInternalState &state) {
    const int point_count = foil.foil_shape.n;
    Eigen::Matrix2Xd updated_surface_vortex(2, point_count);
    for (int i = 0; i < point_count; i++) {
        updated_surface_vortex(0, i) = state.viscous.qvis[i];
        updated_surface_vortex(1, i) = state.inviscid.qinvMatrix(1, i);
    }
    return updated_surface_vortex;
}

FoilAerodynamicCache buildUnitVorticityDistributions(const Foil &foil, const FlowState &analysis_state,
                                                     const XFoilInternalState::InviscidFieldState &state) {
    FoilAerodynamicCache cache = state.cache;

    Logger::instance().write("   Calculating unit vorticity distributions ...\n");

    Eigen::MatrixXd dpsi_dgam = Eigen::MatrixXd::Zero(foil.foil_shape.n + 1, foil.foil_shape.n + 1);
    Eigen::Matrix2Xd psi      = Eigen::Matrix2Xd::Zero(2, foil.foil_shape.n + 1);

    for (int i = 0; i < foil.foil_shape.n; i++) {
        PsiResult psi_result =
            psilin(foil, i, foil.foil_shape.points.col(i), foil.foil_shape.normal_vector.col(i), true, cache.gamu,
                   state.surfaceVortex, analysis_state.alpha, analysis_state.qinf, foil.foil_shape.angle_panel);
        dpsi_dgam.row(i).head(foil.foil_shape.n) = psi_result.dzdg;
        cache.bij.row(i).head(foil.foil_shape.n) = -psi_result.dzdm;
    }
    psi.leftCols(foil.foil_shape.n) =
        -analysis_state.qinf * Eigen::Matrix2d{{0.0, 1.0}, {-1.0, 0.0}} * foil.foil_shape.points;
    dpsi_dgam.col(foil.foil_shape.n)                    = -Eigen::VectorXd::Ones(foil.foil_shape.n);
    dpsi_dgam.row(foil.foil_shape.n)                    = Eigen::VectorXd::Zero(foil.foil_shape.n + 1);
    dpsi_dgam(foil.foil_shape.n, 0)                     = 1;
    dpsi_dgam(foil.foil_shape.n, foil.foil_shape.n - 1) = 1;

    cache.bij.row(foil.foil_shape.n).head(foil.foil_shape.n) = Eigen::VectorXd::Zero(foil.foil_shape.n);

    if (foil.edge.sharp) {
        const double ag1  = atan2(-foil.foil_shape.dpoints_ds.col(0).y(), -foil.foil_shape.dpoints_ds.col(0).x());
        const double ag2  = MathUtil::atanc(foil.foil_shape.dpoints_ds.col(foil.foil_shape.n - 1).y(),
                                            foil.foil_shape.dpoints_ds.col(foil.foil_shape.n - 1).x(), ag1);
        const double abis = 0.5 * (ag1 + ag2);

        Eigen::Vector2d bis_vector{cos(abis), sin(abis)};
        const double dsmin = std::min(
            (foil.foil_shape.points.col(0) - foil.foil_shape.points.col(1)).norm(),
            (foil.foil_shape.points.col(foil.foil_shape.n - 1) - foil.foil_shape.points.col(foil.foil_shape.n - 2))
                .norm());

        const Eigen::Vector2d bis = foil.edge.point_te - 0.1 * dsmin * bis_vector;
        const Eigen::Vector2d normal_bis{-bis_vector.y(), bis_vector.x()};

        PsiResult psi_result = psilin(foil, -1, bis, normal_bis, true, cache.gamu, state.surfaceVortex,
                                      analysis_state.alpha, analysis_state.qinf, foil.foil_shape.angle_panel);

        dpsi_dgam.row(foil.foil_shape.n - 1).head(foil.foil_shape.n) = psi_result.dzdg;
        cache.bij.row(foil.foil_shape.n - 1).head(foil.foil_shape.n) = -psi_result.dzdm;
        psi.col(foil.foil_shape.n - 1)                               = -bis_vector;
    }

    cache.psi_gamma_lu = Eigen::FullPivLU<Eigen::MatrixXd>(dpsi_dgam);
    cache.gamu.resize(2, foil.foil_shape.n + 1);
    cache.gamu.row(0) = cache.psi_gamma_lu.solve(psi.row(0).transpose()).transpose();
    cache.gamu.row(1) = cache.psi_gamma_lu.solve(psi.row(1).transpose()).transpose();
    return cache;
}

bool buildSourceInfluenceMatrix(const Foil &foil, const FlowState &analysis_state,
                                XFoilInternalState::InviscidFieldState &state) {
    Logger::instance().write("   Calculating source influence matrix ...\n");
    const int point_count = foil.foil_shape.n;

    auto hasAirfoilInfluenceMatrix = [&]() {
        const int total_nodes = point_count + foil.wake_shape.n;
        if (state.cache.dij.rows() < total_nodes || state.cache.dij.cols() < total_nodes) {
            return false;
        }
        const auto block = state.cache.dij.block(0, 0, point_count, point_count);
        return block.allFinite() && block.cwiseAbs().maxCoeff() > 0.0;
    };

    if (!hasAirfoilInfluenceMatrix()) {
        state.cache.bij.block(0, 0, point_count + 1, point_count) =
            state.cache.psi_gamma_lu.solve(state.cache.bij.block(0, 0, point_count + 1, point_count)).eval();
        state.cache.dij.block(0, 0, point_count, point_count) = state.cache.bij.block(0, 0, point_count, point_count);
    }

    for (int i = 0; i < point_count; i++) {
        PsiResult psi_result = pswlin(foil, i, foil.foil_shape.points.col(i), foil.foil_shape.normal_vector.col(i),
                                      foil.wake_shape.angle_panel);
        state.cache.bij.row(i).segment(point_count, foil.wake_shape.n) =
            -psi_result.dzdm.segment(point_count, foil.wake_shape.n).transpose();
    }

    state.cache.bij.row(point_count).segment(point_count, foil.wake_shape.n).setZero();

    state.cache.bij.block(0, point_count, point_count + 1, foil.wake_shape.n) =
        state.cache.psi_gamma_lu.solve(state.cache.bij.block(0, point_count, point_count + 1, foil.wake_shape.n))
            .eval();
    state.cache.dij.block(0, point_count, point_count, foil.wake_shape.n) =
        state.cache.bij.block(0, point_count, point_count, foil.wake_shape.n);

    Eigen::MatrixXd cij = Eigen::MatrixXd::Zero(foil.wake_shape.n, point_count);
    for (int i = point_count; i < point_count + foil.wake_shape.n; i++) {
        const int iw = i - point_count;
        PsiResult psi_result =
            psilin(foil, i, foil.wake_shape.points.col(i), foil.wake_shape.normal_vector.col(i), true, state.cache.gamu,
                   state.surfaceVortex, analysis_state.alpha, analysis_state.qinf, foil.wake_shape.angle_panel);
        cij.row(iw)                              = psi_result.dqdg.head(point_count).transpose();
        state.cache.dij.row(i).head(point_count) = psi_result.dqdm.head(point_count).transpose();

        psi_result = pswlin(foil, i, foil.wake_shape.points.col(i), foil.wake_shape.normal_vector.col(i),
                            foil.wake_shape.angle_panel);
        state.cache.dij.row(i).segment(point_count, foil.wake_shape.n) =
            psi_result.dqdm.segment(point_count, foil.wake_shape.n).transpose();
    }

    state.cache.dij.block(point_count, 0, foil.wake_shape.n, point_count) +=
        cij * state.cache.dij.topLeftCorner(point_count, point_count);
    state.cache.dij.block(point_count, point_count, foil.wake_shape.n, foil.wake_shape.n) +=
        cij * state.cache.bij.block(0, point_count, point_count, foil.wake_shape.n);
    state.cache.dij.row(point_count) = state.cache.dij.row(point_count - 1);
    return true;
}

Eigen::VectorXd buildViscousTangentialVelocity(const BoundaryLayer &boundary_layer, const Eigen::VectorXd &base_qvis) {
    Eigen::VectorXd updated_qvis = base_qvis;
    for (int is = 1; is <= 2; is++) {
        const int limit = boundary_layer.readSideStationCount(is) - 1;
        for (int ibl = 0; ibl < limit; ++ibl) {
            const int i = boundary_layer.stationToPanel(is, ibl);
            updated_qvis[i] =
                boundary_layer.panelInfluenceFactor(is, ibl) * boundary_layer.readStationModel(is, ibl).edgeVelocity;
        }
    }
    return updated_qvis;
}

Eigen::Matrix2Xd buildWakeInviscidVelocity(const Foil &foil, const Eigen::Matrix2Xd &base_qinvu,
                                           const Eigen::Matrix2Xd &gamu, const Eigen::Matrix2Xd &surface_vortex,
                                           double alpha, double qinf) {
    const int point_count          = foil.foil_shape.n;
    Eigen::Matrix2Xd updated_qinvu = base_qinvu;

    if (point_count >= 1 && point_count < updated_qinvu.cols()) {
        updated_qinvu.col(point_count) = updated_qinvu.col(point_count - 1);
    }

    for (int i = point_count + 1; i < point_count + foil.wake_shape.n; i++) {
        updated_qinvu.col(i) = psilin(foil, i, foil.wake_shape.points.col(i), foil.wake_shape.normal_vector.col(i),
                                      false, gamu, surface_vortex, alpha, qinf, foil.wake_shape.angle_panel)
                                   .qtan;
    }
    return updated_qinvu;
}

void prepareViscousGeometry(XFoilInternalState &state, const FlowState &analysis_state, Foil &foil,
                            BoundaryLayer &boundary_layer, bool has_panel_map) {
    foil.xyWake(foil.wake_shape.n, state.inviscid.cache.gamu, state.inviscid.surfaceVortex, analysis_state.alpha,
                analysis_state.qinf);
    state.inviscid.qinvu =
        buildWakeInviscidVelocity(foil, state.inviscid.qinvu, state.inviscid.cache.gamu, state.inviscid.surfaceVortex,
                                  analysis_state.alpha, analysis_state.qinf);
    state.inviscid.qinvMatrix = InviscidSolver::qiset(analysis_state.alpha, state.inviscid.qinvu);

    const auto new_stagnation = StagnationFeature(state.inviscid.surfaceVortex, foil.foil_shape.spline_length);

    if (!has_panel_map) {
        if (!new_stagnation.found()) {
            Logger::instance().write("stfind: Stagnation point not found. Continuing ...\n");
        }
        boundary_layer.setStagnationState(new_stagnation);
        state.viscous.stagnation = new_stagnation;
        boundary_layer.buildPanelMap(foil.foil_shape.n, foil.wake_shape.n);
        boundary_layer.rebuildArcLengthCoordinates(foil);
        boundary_layer.buildSystemMapping();
    }
}

void assignCurrentInviscidEdgeVelocity(const XFoilInternalState &state, BoundaryLayer &boundary_layer) {
    const auto inviscid_edge_velocity = boundary_layer.computeInviscidEdgeVelocity(state.inviscid.qinvMatrix);
    boundary_layer.assignInviscidEdgeVelocity(inviscid_edge_velocity);
}

void ensureBoundaryLayerEdgeSeed(const BoundaryLayer &boundary_layer, bool bl_initialized) {
    if (!bl_initialized) {
        const_cast<BoundaryLayer &>(boundary_layer).seedEdgeVelocityFromInviscid();
    }
}

void restoreConvergedOperatingPoint(XFoilResult &result, const XFoilInternalState &state,
                                    const FlowState &analysis_state, const BoundaryLayer &boundary_layer,
                                    int point_count, int total_nodes_with_wake) {
    if (analysis_state.viscous) {
        result.viscousCp  = InviscidSolver::cpcalc(total_nodes_with_wake, state.viscous.qvis, analysis_state.qinf,
                                                   analysis_state.currentMach);
        result.inviscidCp = InviscidSolver::cpcalc(total_nodes_with_wake, state.inviscid.qinvMatrix.row(0).transpose(),
                                                   analysis_state.qinf, analysis_state.currentMach);
        return;
    }

    result.inviscidCp = InviscidSolver::cpcalc(point_count, state.inviscid.qinvMatrix.row(0).transpose(),
                                               analysis_state.qinf, analysis_state.currentMach);
}

void refreshViscousFlowFields(XFoilInternalState &state, const FlowState &analysis_state, const Foil &foil,
                              BoundaryLayer &boundary_layer) {
    state.viscous.qvis           = buildViscousTangentialVelocity(boundary_layer, state.viscous.qvis);
    state.inviscid.surfaceVortex = buildSurfaceVortex(foil, state);
    boundary_layer.moveStagnation(state.inviscid.surfaceVortex, foil.foil_shape.spline_length, foil,
                                  state.inviscid.qinvMatrix, state.viscous.stagnation);
}

} // namespace xfoil_flowfield
