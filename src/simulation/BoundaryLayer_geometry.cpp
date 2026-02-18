#include "simulation/boundary_layer_geometry.hpp"

#include <algorithm>
#include <iostream>

#include "XFoil.h"
#include "core/math_util.hpp"
#include "infrastructure/logger.hpp"

using Eigen::Matrix2Xd;
using Eigen::VectorXd;

BoundaryLayerGeometry::BoundaryLayerGeometry(
    SidePair<BoundaryLayerLattice>& lattice, Eigen::VectorXd& wgap,
    int& stagnationIndex, double& stagnationSst)
    : lattice_(lattice),
      wgap_(wgap),
      stagnationIndex_(stagnationIndex),
      stagnationSst_(stagnationSst) {}

bool BoundaryLayerGeometry::iblpan(int point_count, int wake_point_count) {

  const int lattice_size = point_count + wake_point_count;
  lattice_.top.resize(lattice_size);
  lattice_.bottom.resize(lattice_size);

  for (int i = 0; i <= stagnationIndex_; i++) {
    lattice_.top.stationToPanel[i] = stagnationIndex_ - i;
    lattice_.top.panelInfluenceFactor[i] = 1.0;
  }

  lattice_.top.trailingEdgeIndex = stagnationIndex_;
  lattice_.top.stationCount = lattice_.top.trailingEdgeIndex + 2;

  for (int index = 0; index <= point_count - stagnationIndex_; ++index) {
    lattice_.bottom.stationToPanel[index] = stagnationIndex_ + 1 + index;
    lattice_.bottom.panelInfluenceFactor[index] = -1.0;
  }

  lattice_.bottom.trailingEdgeIndex = point_count - stagnationIndex_ - 2;

  for (int iw = 0; iw < wake_point_count; iw++) {
    const int panel = point_count + iw;
    const int index = lattice_.bottom.trailingEdgeIndex + iw + 2;
    lattice_.bottom.stationToPanel[index - 1] = panel;
    lattice_.bottom.panelInfluenceFactor[index - 1] = -1.0;
  }

  lattice_.bottom.stationCount =
      lattice_.bottom.trailingEdgeIndex + wake_point_count + 2;

  for (int iw = 0; iw < wake_point_count; iw++) {
    lattice_.top.stationToPanel[lattice_.top.trailingEdgeIndex + iw + 1] =
        lattice_.bottom.stationToPanel[lattice_.bottom.trailingEdgeIndex + iw + 1];
    lattice_.top.panelInfluenceFactor[lattice_.top.trailingEdgeIndex + iw + 1] =
        1.0;
  }

  return true;
}

bool BoundaryLayerGeometry::iblsys(int& nsys) {
  int iv = 0;
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 0; ibl < lattice_.get(is).stationCount - 1; ++ibl) {
      ++iv;
      lattice_.get(is).stationToSystem[ibl] = iv;
    }
  }

  nsys = iv + 1;

  return true;
}

StagnationResult BoundaryLayerGeometry::stfind(
    const Eigen::Matrix2Xd& surface_vortex,
    const Eigen::VectorXd& spline_length) const {
  int stagnation_index = 0;
  bool found = false;
  const int point_count = static_cast<int>(surface_vortex.cols());
  for (int i = 0; i < point_count - 1; ++i) {
    if (surface_vortex(0, i) >= 0.0 && surface_vortex(0, i + 1) < 0.0) {
      stagnation_index = i;
      found = true;
      break;
    }
  }

  if (!found) {
    stagnation_index = point_count / 2;
  }

  StagnationResult result;
  result.stagnationIndex = stagnation_index;
  result.found = found;
  const double dgam = surface_vortex(0, stagnation_index + 1) -
                      surface_vortex(0, stagnation_index);
  const double ds = spline_length[stagnation_index + 1] -
                    spline_length[stagnation_index];

  if (surface_vortex(0, stagnation_index) <
      -surface_vortex(0, stagnation_index + 1)) {
    result.sst = spline_length[stagnation_index] -
                 ds * (surface_vortex(0, stagnation_index) / dgam);
  } else {
    result.sst =
        spline_length[stagnation_index + 1] -
        ds * (surface_vortex(0, stagnation_index + 1) / dgam);
  }

  if (result.sst <= spline_length[stagnation_index])
    result.sst = spline_length[stagnation_index] + 0.0000001;
  if (result.sst >= spline_length[stagnation_index + 1])
    result.sst = spline_length[stagnation_index + 1] - 0.0000001;

  result.sst_go =
      (result.sst - spline_length[stagnation_index + 1]) / dgam;
  result.sst_gp = (spline_length[stagnation_index] - result.sst) / dgam;

  return result;
}

bool BoundaryLayerGeometry::stmove(const Eigen::Matrix2Xd& surface_vortex,
                                   const Eigen::VectorXd& spline_length,
                                   const Foil& foil,
                                   const Eigen::Matrix2Xd& qinv_matrix,
                                   StagnationResult& stagnation,
                                   bool& lipan,
                                   int& nsys) {
  const int previous = stagnationIndex_;
  const auto stagnation_result = stfind(surface_vortex, spline_length);
  if (!stagnation_result.found) {
    Logger::instance().write(
        "stfind: Stagnation point not found. Continuing ...\n");
  }
  stagnationIndex_ = stagnation_result.stagnationIndex;
  stagnation = stagnation_result;
  stagnationSst_ = stagnation_result.sst;

  if (previous == stagnationIndex_) {
    xicalc(foil);
  } else {
    if (iblpan(foil.foil_shape.n, foil.wake_shape.n)) {
      lipan = true;
    }
    const auto inviscid_edge_velocity = uicalc(qinv_matrix);
    lattice_.top.inviscidEdgeVelocityMatrix = inviscid_edge_velocity.top;
    lattice_.bottom.inviscidEdgeVelocityMatrix = inviscid_edge_velocity.bottom;
    xicalc(foil);
    iblsys(nsys);

    if (stagnationIndex_ > previous) {
      const int delta = stagnationIndex_ - previous;

      lattice_.top.profiles.transitionIndex += delta;
      lattice_.bottom.profiles.transitionIndex -= delta;

      for (int ibl = lattice_.top.stationCount - 2; ibl >= delta; --ibl) {
        copyStationState(1, ibl, ibl - delta);
      }

      const double dudx =
          lattice_.top.profiles.edgeVelocity[delta] /
          lattice_.top.arcLengthCoordinates[delta];
      for (int ibl = delta; ibl >= 1; --ibl) {
        copyStationState(1, ibl - 1, delta);
        lattice_.top.profiles.edgeVelocity[ibl - 1] =
            dudx * lattice_.top.arcLengthCoordinates[ibl - 1];
      }

      for (int ibl = 0; ibl < lattice_.bottom.stationCount - 1; ++ibl) {
        copyStationState(2, ibl, ibl + delta);
      }
    } else {
      const int delta = previous - stagnationIndex_;

      lattice_.top.profiles.transitionIndex -= delta;
      lattice_.bottom.profiles.transitionIndex += delta;

      for (int ibl = lattice_.bottom.stationCount - 1; ibl >= delta + 1;
           --ibl) {
        copyStationState(2, ibl - 1, (ibl - delta) - 1);
      }

      const double dudx =
      lattice_.bottom.profiles.edgeVelocity[delta] /
          lattice_.bottom.arcLengthCoordinates[delta];
      for (int ibl = delta; ibl >= 1; --ibl) {
        copyStationState(2, ibl - 1, delta);
        lattice_.bottom.profiles.edgeVelocity[ibl - 1] =
            dudx * lattice_.bottom.arcLengthCoordinates[ibl - 1];
      }

      for (int ibl = 0; ibl < lattice_.top.stationCount - 1; ++ibl) {
        copyStationState(1, ibl, ibl + delta);
      }
    }
  }

  for (int is = 1; is <= 2; ++is) {
    for (int ibl = 0; ibl < lattice_.get(is).stationCount - 1; ++ibl) {
      lattice_.get(is).profiles.massFlux[ibl] =
          lattice_.get(is).profiles.displacementThickness[ibl] *
          lattice_.get(is).profiles.edgeVelocity[ibl];
    }
  }

  return true;
}

SidePair<Eigen::Matrix2Xd> BoundaryLayerGeometry::uicalc(
    const Eigen::Matrix2Xd& qinv_matrix) const {
  //--------------------------------------------------------------
  //     sets inviscid ue from panel inviscid tangential velocity
  //--------------------------------------------------------------
  SidePair<Matrix2Xd> inviscid_matrix;
  for (int side = 1; side <= 2; ++side) {
    inviscid_matrix.get(side) =
        Matrix2Xd::Zero(2, lattice_.get(side).stationCount);
    for (int stationIndex = 0; stationIndex < lattice_.get(side).stationCount - 1;
         ++stationIndex) {
      const int panelIndex = lattice_.get(side).stationToPanel[stationIndex];
      inviscid_matrix.get(side)(0, stationIndex) =
          lattice_.get(side).panelInfluenceFactor[stationIndex] *
          qinv_matrix(0, panelIndex);
      inviscid_matrix.get(side)(1, stationIndex) =
          lattice_.get(side).panelInfluenceFactor[stationIndex] *
          qinv_matrix(1, panelIndex);
    }
  }

  return inviscid_matrix;
}

bool BoundaryLayerGeometry::xicalc(const Foil& foil) {
  //-------------------------------------------------------------
  //     sets bl arc length array on each airfoil side and wake
  //-------------------------------------------------------------

  const auto arc_lengths =
      computeArcLengthCoordinates(foil, stagnationSst_, lattice_);
  lattice_.top.arcLengthCoordinates = arc_lengths.top;
  lattice_.bottom.arcLengthCoordinates = arc_lengths.bottom;
  wgap_ = computeWakeGap(foil, lattice_.bottom, arc_lengths.bottom);
  return true;
}

SidePair<VectorXd> BoundaryLayerGeometry::computeArcLengthCoordinates(
    const Foil& foil, double stagnationSst,
    const SidePair<BoundaryLayerLattice>& lattice) {
  SidePair<VectorXd> arc_lengths;
  arc_lengths.top = lattice.top.arcLengthCoordinates;
  arc_lengths.bottom = lattice.bottom.arcLengthCoordinates;

  for (int ibl = 0; ibl <= lattice.top.trailingEdgeIndex; ++ibl) {
    arc_lengths.top[ibl] =
        stagnationSst -
        foil.foil_shape.spline_length[lattice.get(1).stationToPanel[ibl]];
  }

  for (int ibl = 0; ibl <= lattice.bottom.trailingEdgeIndex; ++ibl) {
    arc_lengths.bottom[ibl] =
        foil.foil_shape.spline_length[lattice.get(2).stationToPanel[ibl]] -
        stagnationSst;
  }

  // Wake: start from TE, duplicate TE value at first wake station
  arc_lengths.bottom[lattice.bottom.trailingEdgeIndex + 1] =
      arc_lengths.bottom[lattice.bottom.trailingEdgeIndex];
  for (int ibl = lattice.bottom.trailingEdgeIndex + 2;
       ibl < lattice.bottom.stationCount; ++ibl) {
    arc_lengths.bottom[ibl] =
        arc_lengths.bottom[ibl - 1] +
        (foil.wake_shape.points.col(lattice.get(2).stationToPanel[ibl]) -
         foil.wake_shape.points.col(lattice.get(2).stationToPanel[ibl - 1]))
            .norm();
  }
  return arc_lengths;
}

VectorXd BoundaryLayerGeometry::computeWakeGap(
    const Foil& foil, const BoundaryLayerLattice& bottom,
    const VectorXd& bottomArcLengths) {
  //---- trailing edge flap length to te gap ratio
  const double telrat = 2.50;

  //---- set up parameters for te flap cubics

  const int point_count = foil.foil_shape.n;
  const double crosp = MathUtil::cross2(
      foil.foil_shape.dpoints_ds.col(point_count - 1).normalized(),
      foil.foil_shape.dpoints_ds.col(0).normalized());
  double dwdxte = crosp / sqrt(1.0 - crosp * crosp);

  //---- limit cubic to avoid absurd te gap widths
  dwdxte = std::max(dwdxte, -3.0 / telrat);
  dwdxte = std::min(dwdxte, 3.0 / telrat);

  const double aa = 3.0 + telrat * dwdxte;
  const double bb = -2.0 - telrat * dwdxte;

  VectorXd wgap_result = VectorXd::Zero(foil.wake_shape.n);
  if (foil.edge.sharp) {
    return wgap_result;
  }

  else {
    //----- set te flap (wake gap) array (0-based: iw0=0..foil.wake_shape.n-1)
    for (int iw0 = 0; iw0 < foil.wake_shape.n; iw0++) {
      const int te_bot_0b =
          bottom.trailingEdgeIndex;  // 0-based TE for array indexing
      const double zn =
          1.0 -
          (bottomArcLengths[te_bot_0b + (iw0 + 1)] - bottomArcLengths[te_bot_0b]) /
              (telrat * foil.edge.ante);
      if (zn >= 0.0)
        wgap_result[iw0] = foil.edge.ante * (aa + bb * zn) * zn * zn;
    }
  }
  return wgap_result;
}

void BoundaryLayerGeometry::copyStationState(
    int side, int destination, int source) {
  lattice_.get(side).profiles.skinFrictionCoeff[destination] =
      lattice_.get(side).profiles.skinFrictionCoeff[source];
  lattice_.get(side).profiles.momentumThickness[destination] =
      lattice_.get(side).profiles.momentumThickness[source];
  lattice_.get(side).profiles.displacementThickness[destination] =
      lattice_.get(side).profiles.displacementThickness[source];
  lattice_.get(side).profiles.edgeVelocity[destination] =
      lattice_.get(side).profiles.edgeVelocity[source];
}
