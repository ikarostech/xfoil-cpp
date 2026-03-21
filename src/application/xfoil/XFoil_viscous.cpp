#include "infrastructure/logger.hpp"
#include "numerics/math_util.hpp"
#include "application/xfoil/XFoilAnalysis.hpp"
#include "model/coefficient/xfoil_postprocess.hpp"
#include "solver/boundary_layer/boundary_layer_builder.hpp"
#include "solver/boundary_layer/viscous_iteration_ops.hpp"
#include "solver/inviscid/InviscidSolver.hpp"
#include "solver/march/march.hpp"
#include "solver/xfoil/xfoil_flowfield.hpp"
#include "solver/xfoil/viscous_update.hpp"
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
bool XFoilAnalysis::initXFoilGeometry(int fn, const double *fx, const double *fy) {

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

bool XFoilAnalysis::initXFoilAnalysis(double Re, double alpha, double Mach, double NCrit, double XtrTop, double XtrBot,
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

    acrit_ = NCrit;
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
bool XFoilAnalysis::qdcalc() {
    return xfoil_flowfield::buildSourceInfluenceMatrix(foil, analysis_state_,
                                                       state_.inviscid);
}

/** -------------------------------------------------------------
 *     sets panel viscous tangential velocity from viscous ue
 * -------------------------------------------------------------- */
VectorXd XFoilAnalysis::qvfue(const VectorXd &base_qvis) const {
    return xfoil_flowfield::buildViscousTangentialVelocity(boundaryLayer,
                                                           base_qvis);
}

/** ---------------------------------------------------------------
 *      sets inviscid tangential velocity for alpha = 0, 90
 *      on wake due to freestream and airfoil surface vorticity.
 * --------------------------------------------------------------- */
Matrix2Xd XFoilAnalysis::qwcalc(const Foil &foil, const Matrix2Xd &base_qinvu, const Matrix2Xd &gamu,
                        const Matrix2Xd &surface_vortex, double alpha, double qinf) const {
    return xfoil_flowfield::buildWakeInviscidVelocity(
        foil, base_qinvu, gamu, surface_vortex, alpha, qinf);
}

XFoilAnalysis::UpdateResult XFoilAnalysis::update(const XFoilAnalysis::Matrix3x2dVector &vdel) const {
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

bool XFoilAnalysis::viscal() {
    ////--------------------------------------
    //     converges viscous operating point
    ////--------------------------------------
    const int point_count           = foil.foil_shape.n;
    const int total_nodes_with_wake = point_count + foil.wake_shape.n;

    ensurePanelMapAndBoundaryLayerGeometry();
    assignCurrentInviscidEdgeVelocity();
    ensureBoundaryLayerEdgeSeed();
    if (hasConvergedSolution()) {
        restoreConvergedOperatingPoint(point_count, total_nodes_with_wake);
    }
    ensureSourceInfluenceMatrix();

    return true;
}

XFoilAnalysis::ViscalEndResult XFoilAnalysis::ViscalEnd() {
    ViscalEndResult result;
    const int total_nodes_with_wake = foil.foil_shape.n + foil.wake_shape.n;
    result.inviscidCp = InviscidSolver::cpcalc(total_nodes_with_wake, state_.inviscid.qinvMatrix.row(0).transpose(),
                                               analysis_state_.qinf, analysis_state_.currentMach);
    result.viscousCp  = InviscidSolver::cpcalc(total_nodes_with_wake, state_.viscous.qvis, analysis_state_.qinf,
                                               analysis_state_.currentMach);
    return result;
}

bool XFoilAnalysis::ViscousIter() {
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

void XFoilAnalysis::ensureWakeTrajectoryAndInviscidVelocity() {
    xfoil_flowfield::prepareViscousGeometry(state_, analysis_state_, foil,
                                            boundaryLayer, hasPanelMap());
}

void XFoilAnalysis::ensurePanelMapAndBoundaryLayerGeometry() {
    xfoil_flowfield::prepareViscousGeometry(state_, analysis_state_, foil,
                                            boundaryLayer, hasPanelMap());
}

void XFoilAnalysis::assignCurrentInviscidEdgeVelocity() {
    xfoil_flowfield::assignCurrentInviscidEdgeVelocity(state_, boundaryLayer);
}

void XFoilAnalysis::ensureBoundaryLayerEdgeSeed() {
    xfoil_flowfield::ensureBoundaryLayerEdgeSeed(boundaryLayer,
                                                 isBLInitialized());
}

void XFoilAnalysis::restoreConvergedOperatingPoint(int point_count, int total_nodes_with_wake) {
    state_.viscous.qvis = qvfue(state_.viscous.qvis);

    xfoil_flowfield::restoreConvergedOperatingPoint(
        result_, state_, analysis_state_, boundaryLayer, point_count,
        total_nodes_with_wake);
    cpi_ = result_.inviscidCp;
    cpv_ = result_.viscousCp;

    const auto cl_result = clcalc(cmref);
    applyClComputation(cl_result);
    aero_coeffs_.cd = xfoil_postprocess::computeCd(analysis_state_, boundaryLayer,
                                                   isBLInitialized());
}

void XFoilAnalysis::ensureSourceInfluenceMatrix() {
    if (!hasAirfoilInfluenceMatrix() || !hasWakeInfluenceMatrix()) {
        qdcalc();
    }
}

SetblOutputView XFoilAnalysis::initializeBoundaryLayerNewtonSystem() {
    return boundary_layer_iteration_ops::initializeNewtonSystem(
        boundaryLayer, analysis_state_, aero_coeffs_, acrit_, foil,
        state_.viscous.stagnation, state_.inviscid.cache.dij,
        isBLInitialized());
}

XFoilAnalysis::Matrix3x2dVector XFoilAnalysis::solveBoundaryLayerNewtonStep(const SetblOutputView &setbl_output) const {
    return boundary_layer_iteration_ops::solveNewtonStep(boundaryLayer,
                                                         VAccel(),
                                                         setbl_output);
}

void XFoilAnalysis::applyBoundaryLayerIterationUpdate(const UpdateResult &update_result) {
    boundary_layer_iteration_ops::applyIterationUpdate(boundaryLayer,
                                                       update_result);
}

void XFoilAnalysis::updateFreestreamForIteration(const UpdateResult &update_result) {
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

void XFoilAnalysis::refreshViscousFlowFields() {
    xfoil_flowfield::refreshViscousFlowFields(state_, analysis_state_, foil,
                                              boundaryLayer);

    const auto cl_result = clcalc(cmref);
    applyClComputation(cl_result);
    aero_coeffs_.cd = xfoil_postprocess::computeCd(analysis_state_, boundaryLayer,
                                                   isBLInitialized());
}

void XFoilAnalysis::finalizeViscousIteration(const UpdateResult &update_result, double eps1) {
    if (update_result.rmsbl >= eps1) {
        return;
    }

    state_.viscous.convergedAlpha = update_result.analysis_state.alpha;
    state_.viscous.convergedMach  = update_result.analysis_state.currentMach;
    Logger::instance().write("----------CONVERGED----------\n\n");
}

bool XFoilAnalysis::isValidFoilAngles(Matrix2Xd points) {
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

bool XFoilAnalysis::isValidFoilPointSize(Matrix2Xd points) {
    return points.cols() >= 3;
}
