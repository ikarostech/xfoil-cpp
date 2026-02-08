#pragma once

#include <string>

#include "core/side_pair.hpp"
#include "domain/foil/foil.hpp"
#include "simulation/boundary_layer_state.hpp"

class XFoil;

struct StagnationResult {
  int stagnationIndex = 0;
  double sst = 0.0;
  double sst_go = 0.0;
  double sst_gp = 0.0;
  bool found = true;
};

class BoundaryLayerGeometry {
 public:
  BoundaryLayerGeometry(SidePair<BoundaryLayerLattice>& lattice,
                        Eigen::VectorXd& wgap,
                        int& stagnationIndex,
                        double& stagnationSst);

  bool iblpan(int point_count, int wake_point_count);
  bool iblsys(int& nsys);
  StagnationResult stfind(const Eigen::Matrix2Xd& surface_vortex,
                          const Eigen::VectorXd& spline_length) const;
  bool stmove(XFoil& xfoil);
  bool xicalc(const Foil& foil);
  SidePair<Eigen::Matrix2Xd> uicalc(
      const Eigen::Matrix2Xd& qinv_matrix) const;

 private:
  static SidePair<Eigen::VectorXd> computeArcLengthCoordinates(
      const Foil& foil, double stagnationSst,
      const SidePair<BoundaryLayerLattice>& lattice);
  static Eigen::VectorXd computeWakeGap(
      const Foil& foil, const BoundaryLayerLattice& bottom,
      const Eigen::VectorXd& bottomArcLengths);
  void copyStationState(int side, int destination, int source);

  SidePair<BoundaryLayerLattice>& lattice_;
  Eigen::VectorXd& wgap_;
  int& stagnationIndex_;
  double& stagnationSst_;
};
