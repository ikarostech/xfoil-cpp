#pragma once
#include "core/side_pair.hpp"
#include "Eigen/Core"
#include "infrastructure/xfoil_params.h"

struct SetblInputView {
  const bool& lblini;
  SidePairRef<const Eigen::VectorXd> uedg;
  SidePairRef<const Eigen::VectorXd> ctau;
  SidePairRef<const Eigen::VectorXd> thet;
  SidePairRef<const Eigen::VectorXd> dstr;
  SidePairRef<const Eigen::VectorXd> mass;
  SidePairRef<const Eigen::VectorXd> ctq;
  SidePairRef<const int> itran;
};

struct SetblOutputView {
  bool& lblini;
  double& gm1bl;
  double& qinfbl;
  double& tkbl;
  double& tkbl_ms;
  double& rstbl;
  double& rstbl_ms;
  double& hstinv;
  double& hstinv_ms;
  double& reybl;
  double& reybl_re;
  double& reybl_ms;
  double& amcrit;
  SidePairRef<Eigen::VectorXd> uedg;
  SidePairRef<Eigen::VectorXd> ctau;
  SidePairRef<Eigen::VectorXd> thet;
  SidePairRef<Eigen::VectorXd> dstr;
  SidePairRef<Eigen::VectorXd> mass;
  SidePairRef<Eigen::VectorXd> ctq;
  SidePairRef<int> itran;
  std::vector<Eigen::Matrix<double, 3, 2>>& va;
  std::vector<Eigen::Matrix<double, 3, 2>>& vb;
  std::vector<Eigen::Matrix<double, 3, 2>>& vdel;
  double (&vm)[3][IZX][IZX];
  double (&vz)[3][2];
  bool& tran;
  bool& turb;
  bool& wake;
  bool& simi;
  double& xiforc;
};
