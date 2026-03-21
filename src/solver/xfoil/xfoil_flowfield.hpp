#pragma once

#include "Eigen/Core"

#include "application/xfoil/XFoilInternalState.hpp"
#include "application/xfoil/XFoilResult.hpp"
#include "model/flow_state.hpp"
#include "model/foil/foil.hpp"
#include "solver/boundary_layer/boundary_layer.hpp"

namespace xfoil_flowfield {

void updateTrailingEdgeState(const Foil &foil,
                             XFoilInternalState::InviscidFieldState &state);

Eigen::Matrix2Xd buildSurfaceVortex(const Foil &foil,
                                    const XFoilInternalState &state);

FoilAerodynamicCache buildUnitVorticityDistributions(
    const Foil &foil, const FlowState &analysis_state,
    const XFoilInternalState::InviscidFieldState &state);

bool buildSourceInfluenceMatrix(const Foil &foil, const FlowState &analysis_state,
                                XFoilInternalState::InviscidFieldState &state);

Eigen::VectorXd buildViscousTangentialVelocity(const BoundaryLayer &boundary_layer,
                                               const Eigen::VectorXd &base_qvis);

Eigen::Matrix2Xd buildWakeInviscidVelocity(const Foil &foil,
                                           const Eigen::Matrix2Xd &base_qinvu,
                                           const Eigen::Matrix2Xd &gamu,
                                           const Eigen::Matrix2Xd &surface_vortex,
                                           double alpha, double qinf);

void prepareViscousGeometry(XFoilInternalState &state,
                            const FlowState &analysis_state, Foil &foil,
                            BoundaryLayer &boundary_layer, bool has_panel_map);

void assignCurrentInviscidEdgeVelocity(const XFoilInternalState &state,
                                       BoundaryLayer &boundary_layer);

void ensureBoundaryLayerEdgeSeed(const BoundaryLayer &boundary_layer,
                                 bool bl_initialized);

void restoreConvergedOperatingPoint(XFoilResult &result,
                                    const XFoilInternalState &state,
                                    const FlowState &analysis_state,
                                    const BoundaryLayer &boundary_layer,
                                    int point_count,
                                    int total_nodes_with_wake);

void refreshViscousFlowFields(XFoilInternalState &state, const FlowState &analysis_state,
                              const Foil &foil, BoundaryLayer &boundary_layer);

} // namespace xfoil_flowfield
