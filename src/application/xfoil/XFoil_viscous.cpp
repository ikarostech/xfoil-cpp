#include "infrastructure/logger.hpp"
#include "numerics/math_util.hpp"
#include "solver/boundary_layer/blsolve.hpp"
#include "solver/boundary_layer/boundary_layer_builder.hpp"
#include "solver/boundary_layer/initialization/viscous_initializer.hpp"
#include "solver/inviscid/InviscidSolver.hpp"
#include "solver/inviscid/psi.hpp"
#include "solver/march/march.hpp"
#include "solver/xfoil/viscous_update.hpp"
#include "application/xfoil/XFoil.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <numbers>
#include <utility>
using Eigen::Matrix2Xd;
using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::VectorXd;

namespace {
constexpr double kAngleTolerance = 40.0;
} // namespace

/** Loads the Foil's geometry in XFoil,
 *  calculates the normal vectors,
 *  and sets the results in current foil */
bool XFoil::initXFoilGeometry(int fn, const double *fx, const double *fy) {

    Matrix2Xd buffer_points = Matrix2Xd::Zero(2, fn);
    for (int i = 0; i < fn; i++) {
        buffer_points.col(i).x() = fx[i];
        buffer_points.col(i).y() = fy[i];
    }

    if (!isValidFoilPointSize(buffer_points) || !isValidFoilAngles(buffer_points)) {
        Logger::instance().write("Unrecognized foil format");
        return false;
    }

    abcopy(buffer_points);
    state_.inviscid.cache = ggcalc();
    state_.inviscid.qinvu.leftCols(foil.foil_shape.n + 1) = state_.inviscid.cache.gamu;
    return true;
}

bool XFoil::initXFoilAnalysis(double Re, double alpha, double Mach, double NCrit, double XtrTop, double XtrBot,
                              ReynoldsType reType, MachType maType, bool bViscous) {
    setBLInitialized(false);
    invalidatePanelMap();
    invalidateWakeGeometry();
    invalidateConvergedSolution();

    FlowState &state     = analysis_state_;
    state.referenceRe    = Re;
    state.alpha          = alpha * std::numbers::pi / 180.0;
    state.referenceMach  = Mach;
    state.reynoldsType   = reType;
    state.machType       = maType;
    state.controlByAlpha = true;
    state.qinf           = 1.0;
    state.viscous        = bViscous;

    acrit = NCrit;
    boundaryLayer.setTransitionLocations(XtrTop, XtrBot);

    if (Mach > 0.000001) {
        setMach();
    }

    return true;
}

/** --------------------------------------------------------------------
 *	   Calculates current streamfunction psi and tangential velocity
 *	   qtan at panel node or wake node i due to freestream and wake
 *	   sources sig.  also calculates sensitivity vectors dpsi/dsig
 *	   (dzdm) and dqtan/dsig (dqdm).
 *
 *			airfoil:  1   < i < n
 *			wake:	  n+1 < i < n+foil.wake_shape.n
 *-------------------------------------------------------------------- */
/** ------------------------------------------------------
 *         calculates source panel influence coefficient
 * 	   matrix for current airfoil and wake geometry.
 * ------------------------------------------------------ */
bool XFoil::qdcalc() {
    // TRACE("calculating source influence matrix ...\n");
    Logger::instance().write("   Calculating source influence matrix ...\n");
    const int point_count = foil.foil_shape.n;

    if (!hasAirfoilInfluenceMatrix()) {
        //----- calculate source influence matrix for airfoil surface
        state_.inviscid.cache.bij.block(0, 0, point_count + 1, point_count) =
            state_.inviscid.cache.psi_gamma_lu
                .solve(state_.inviscid.cache.bij.block(0, 0, point_count + 1, point_count))
                .eval();

        //------- store resulting dgam/dsig = dqtan/dsig vector
        state_.inviscid.cache.dij.block(0, 0, point_count, point_count) =
            state_.inviscid.cache.bij.block(0, 0, point_count, point_count);
    }

    //---- set up coefficient matrix of dpsi/dm on airfoil surface
    for (int i = 0; i < point_count; i++) {
        PsiResult psi_result = pswlin(foil, i, foil.foil_shape.points.col(i), foil.foil_shape.normal_vector.col(i),
                                      foil.wake_shape.angle_panel);
        state_.inviscid.cache.bij.row(i).segment(point_count, foil.wake_shape.n) =
            -psi_result.dzdm.segment(point_count, foil.wake_shape.n).transpose();
    }

    //---- set up kutta condition (no direct source influence)

    state_.inviscid.cache.bij.row(point_count).segment(point_count, foil.wake_shape.n).setZero();

    //---- multiply by inverse of factored dpsi/dgam matrix
    state_.inviscid.cache.bij.block(0, point_count, point_count + 1, foil.wake_shape.n) =
        state_.inviscid.cache.psi_gamma_lu
            .solve(state_.inviscid.cache.bij.block(0, point_count, point_count + 1, foil.wake_shape.n))
            .eval();
    //---- set the source influence matrix for the wake sources
    state_.inviscid.cache.dij.block(0, point_count, point_count, foil.wake_shape.n) =
        state_.inviscid.cache.bij.block(0, point_count, point_count, foil.wake_shape.n);

    //**** now we need to calculate the influence of sources on the wake
    // velocities

    //---- calculate dqtan/dgam and dqtan/dsig at the wake points
    MatrixXd cij = MatrixXd::Zero(foil.wake_shape.n, point_count);
    for (int i = point_count; i < point_count + foil.wake_shape.n; i++) {
        int iw = i - point_count;
        //------ airfoil contribution at wake panel node
        PsiResult psi_result = psilin(foil, i, foil.wake_shape.points.col(i), foil.wake_shape.normal_vector.col(i),
                                      true, state_.inviscid.cache.gamu, state_.inviscid.surfaceVortex,
                                      analysis_state_.alpha, analysis_state_.qinf, foil.wake_shape.angle_panel);
        cij.row(iw)          = psi_result.dqdg.head(point_count).transpose();
        state_.inviscid.cache.dij.row(i).head(point_count) = psi_result.dqdm.head(point_count).transpose();
        //------ wake contribution
        psi_result = pswlin(foil, i, foil.wake_shape.points.col(i), foil.wake_shape.normal_vector.col(i),
                            foil.wake_shape.angle_panel);
        state_.inviscid.cache.dij.row(i).segment(point_count, foil.wake_shape.n) =
            psi_result.dqdm.segment(point_count, foil.wake_shape.n).transpose();
    }

    //---- add on effect of all sources on airfoil vorticity which effects wake
    // qtan
    state_.inviscid.cache.dij.block(point_count, 0, foil.wake_shape.n, point_count) +=
        cij * state_.inviscid.cache.dij.topLeftCorner(point_count, point_count);

    state_.inviscid.cache.dij.block(point_count, point_count, foil.wake_shape.n, foil.wake_shape.n) +=
        cij * state_.inviscid.cache.bij.block(0, point_count, point_count, foil.wake_shape.n);

    //---- make sure first wake point has same velocity as trailing edge
    state_.inviscid.cache.dij.row(point_count) = state_.inviscid.cache.dij.row(point_count - 1);

    return true;
}

/** -------------------------------------------------------------
 *     sets panel viscous tangential velocity from viscous ue
 * -------------------------------------------------------------- */
VectorXd XFoil::qvfue(const VectorXd &base_qvis) const {
    VectorXd updated_qvis = base_qvis;
    for (int is = 1; is <= 2; is++) {
        const int limit = boundaryLayer.readSideStationCount(is) - 1;
        for (int ibl = 0; ibl < limit; ++ibl) {
            int i = boundaryLayer.stationToPanel(is, ibl);
            updated_qvis[i] =
                boundaryLayer.panelInfluenceFactor(is, ibl) * boundaryLayer.readStationModel(is, ibl).edgeVelocity;
        }
    }

    return updated_qvis;
}

/** ---------------------------------------------------------------
 *      sets inviscid tangential velocity for alpha = 0, 90
 *      on wake due to freestream and airfoil surface vorticity.
 * --------------------------------------------------------------- */
Matrix2Xd XFoil::qwcalc(const Foil &foil, const Matrix2Xd &base_qinvu, const Matrix2Xd &gamu,
                        const Matrix2Xd &surface_vortex, double alpha, double qinf) const {
    const int point_count   = foil.foil_shape.n;
    Matrix2Xd updated_qinvu = base_qinvu;

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

XFoil::UpdateResult XFoil::update(const XFoil::Matrix3x2dVector &vdel) const {
    return BoundaryLayerViscousUpdate::run(
        ViscousUpdateInput{boundaryLayer,
                           analysis_state_,
                           aero_coeffs_,
                           state_.operatingPointCoupling.machPerLift,
                           BoundaryLayerAerodynamicContext{
                               state_.inviscid.cache.dij,
                               foil.foil_shape.points,
                               analysis_state_.alpha,
                               analysis_state_.qinf,
                               analysis_state_.currentMach,
                               analysis_state_.controlByAlpha}},
        vdel);
}

bool XFoil::viscal() {
    ////--------------------------------------
    //     converges viscous operating point
    ////--------------------------------------
    const int point_count           = foil.foil_shape.n;
    const int total_nodes_with_wake = point_count + foil.wake_shape.n;

    ensureWakeTrajectoryAndInviscidVelocity();
    ensurePanelMapAndBoundaryLayerGeometry();
    assignCurrentInviscidEdgeVelocity();
    ensureBoundaryLayerEdgeSeed();
    if (hasConvergedSolution()) {
        restoreConvergedOperatingPoint(point_count, total_nodes_with_wake);
    }
    ensureSourceInfluenceMatrix();

    return true;
}

XFoil::ViscalEndResult XFoil::ViscalEnd() {
    ViscalEndResult result;
    const int total_nodes_with_wake = foil.foil_shape.n + foil.wake_shape.n;
    result.inviscidCp = InviscidSolver::cpcalc(total_nodes_with_wake, state_.inviscid.qinvMatrix.row(0).transpose(),
                                               analysis_state_.qinf, analysis_state_.currentMach);
    result.viscousCp  = InviscidSolver::cpcalc(total_nodes_with_wake, state_.viscous.qvis, analysis_state_.qinf,
                                               analysis_state_.currentMach);
    return result;
}

bool XFoil::ViscousIter() {
    //	Performs one iteration
    double eps1 = 0.0001;

    auto setbl_output        = initializeBoundaryLayerNewtonSystem();
    auto vdel                = solveBoundaryLayerNewtonStep(setbl_output);
    const auto update_result = update(vdel);
    applyBoundaryLayerIterationUpdate(update_result);
    updateFreestreamForIteration(update_result);
    refreshViscousFlowFields();
    finalizeViscousIteration(update_result, eps1);

    return true;
}

void XFoil::ensureWakeTrajectoryAndInviscidVelocity() {
    foil.xyWake(foil.wake_shape.n, state_.inviscid.cache.gamu, state_.inviscid.surfaceVortex, analysis_state_.alpha,
                analysis_state_.qinf);
    state_.inviscid.qinvu =
        qwcalc(foil, state_.inviscid.qinvu, state_.inviscid.cache.gamu, state_.inviscid.surfaceVortex,
               analysis_state_.alpha, analysis_state_.qinf);
    state_.inviscid.qinvMatrix = InviscidSolver::qiset(analysis_state_.alpha, state_.inviscid.qinvu);
}

void XFoil::ensurePanelMapAndBoundaryLayerGeometry() {
    if (hasPanelMap()) {
        return;
    }

    const auto new_stagnation =
        boundaryLayer.findStagnation(state_.inviscid.surfaceVortex, foil.foil_shape.spline_length);
    if (!new_stagnation.found) {
        Logger::instance().write("stfind: Stagnation point not found. Continuing ...\n");
    }

    boundaryLayer.setStagnationState(new_stagnation);
    state_.viscous.stagnation = new_stagnation;
    boundaryLayer.buildPanelMap(foil.foil_shape.n, foil.wake_shape.n);
    boundaryLayer.rebuildArcLengthCoordinates(foil);
    boundaryLayer.buildSystemMapping();
}

void XFoil::assignCurrentInviscidEdgeVelocity() {
    const auto inviscid_edge_velocity = boundaryLayer.computeInviscidEdgeVelocity(state_.inviscid.qinvMatrix);
    boundaryLayer.assignInviscidEdgeVelocity(inviscid_edge_velocity);
}

void XFoil::ensureBoundaryLayerEdgeSeed() {
    if (!isBLInitialized()) {
        boundaryLayer.seedEdgeVelocityFromInviscid();
    }
}

void XFoil::restoreConvergedOperatingPoint(int point_count, int total_nodes_with_wake) {
    state_.viscous.qvis = qvfue(state_.viscous.qvis);

    if (analysis_state_.viscous) {
        cpv = InviscidSolver::cpcalc(total_nodes_with_wake, state_.viscous.qvis, analysis_state_.qinf,
                                     analysis_state_.currentMach);
        cpi = InviscidSolver::cpcalc(total_nodes_with_wake, state_.inviscid.qinvMatrix.row(0).transpose(),
                                     analysis_state_.qinf, analysis_state_.currentMach);
    } else {
        cpi = InviscidSolver::cpcalc(point_count, state_.inviscid.qinvMatrix.row(0).transpose(), analysis_state_.qinf,
                                     analysis_state_.currentMach);
    }

    const auto cl_result = clcalc(cmref);
    applyClComputation(cl_result);
    aero_coeffs_.cd = cdcalc();
}

void XFoil::ensureSourceInfluenceMatrix() {
    if (!hasAirfoilInfluenceMatrix() || !hasWakeInfluenceMatrix()) {
        qdcalc();
    }
}

SetblOutputView XFoil::initializeBoundaryLayerNewtonSystem() {
    auto setbl_output =
        BoundaryLayerInitializer::run(boundaryLayer, analysis_state_, aero_coeffs_, acrit, foil,
                                      state_.viscous.stagnation, state_.inviscid.cache.dij, isBLInitialized());
    BoundaryLayerInitializer::applyOutput(boundaryLayer, setbl_output);
    return setbl_output;
}

XFoil::Matrix3x2dVector XFoil::solveBoundaryLayerNewtonStep(const SetblOutputView &setbl_output) const {
    Blsolve solver;
    SidePair<int> ivte{boundaryLayer.trailingEdgeSystemIndex(1), boundaryLayer.trailingEdgeSystemIndex(2)};
    return solver.solve(boundaryLayer.systemSize(), ivte, VAccel(), setbl_output.bl_newton_system);
}

void XFoil::applyBoundaryLayerIterationUpdate(const UpdateResult &update_result) {
    boundaryLayer.applyProfiles(update_result.profiles);
}

void XFoil::updateFreestreamForIteration(const UpdateResult &update_result) {
    if (update_result.analysis_state.controlByAlpha) {
        state_.operatingPointCoupling.machPerLift =
            getActualMach(update_result.aero_coeffs.cl, update_result.analysis_state.machType);
        state_.operatingPointCoupling.reynoldsPerLift =
            getActualReynolds(update_result.aero_coeffs.cl, update_result.analysis_state.reynoldsType);
        return;
    }

    state_.inviscid.qinvMatrix = InviscidSolver::qiset(update_result.analysis_state.alpha, state_.inviscid.qinvu);
    assignCurrentInviscidEdgeVelocity();
}

void XFoil::refreshViscousFlowFields() {
    state_.viscous.qvis           = qvfue(state_.viscous.qvis);
    state_.inviscid.surfaceVortex = gamqv();
    boundaryLayer.moveStagnation(state_.inviscid.surfaceVortex, foil.foil_shape.spline_length, foil,
                                 state_.inviscid.qinvMatrix, state_.viscous.stagnation);

    const auto cl_result = clcalc(cmref);
    applyClComputation(cl_result);
    aero_coeffs_.cd = cdcalc();
}

void XFoil::finalizeViscousIteration(const UpdateResult &update_result, double eps1) {
    if (update_result.rmsbl >= eps1) {
        return;
    }

    state_.viscous.convergedAlpha = update_result.analysis_state.alpha;
    state_.viscous.convergedMach  = update_result.analysis_state.currentMach;
    Logger::instance().write("----------CONVERGED----------\n\n");
}

bool XFoil::isValidFoilAngles(Matrix2Xd points) {
    auto cang = [&](Matrix2Xd points) {
        double max_angle = 0;
        //---- go over each point, calculating corner angle
        for (int i = 1; i < points.cols() - 1; i++) {
            Vector2d delta_former = points.col(i) - points.col(i - 1);
            Vector2d delta_later  = points.col(i) - points.col(i + 1);
            double sin         = MathUtil::cross2(delta_later, delta_former) / delta_former.norm() / delta_later.norm();
            double delta_angle = asin(sin) * 180.0 / std::numbers::pi;
            max_angle          = std::max(fabs(delta_angle), max_angle);
        }
        return max_angle;
    };
    double max_angle = cang(points);
    return max_angle <= kAngleTolerance;
}

bool XFoil::isValidFoilPointSize(Matrix2Xd points) {
    return points.cols() >= 3;
}
