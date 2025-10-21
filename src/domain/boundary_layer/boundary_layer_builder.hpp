#pragma once
#include "core/side_pair.hpp"
#include "Eigen/Core"
#include "infrastructure/xfoil_params.h"

struct SetblInputView {
    const bool& lblini;
    const SidePair<Eigen::VectorXd>& uedg;
    const SidePair<Eigen::VectorXd>& ctau;
    const SidePair<Eigen::VectorXd>& thet;
    const SidePair<Eigen::VectorXd>& dstr;
    const SidePair<Eigen::VectorXd>& mass;
    const SidePair<Eigen::VectorXd>& ctq;
    const SidePair<int>& itran;
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
    SidePair<Eigen::VectorXd>& uedg;
    SidePair<Eigen::VectorXd>& ctau;
    SidePair<Eigen::VectorXd>& thet;
    SidePair<Eigen::VectorXd>& dstr;
    SidePair<Eigen::VectorXd>& mass;
    SidePair<Eigen::VectorXd>& ctq;
    SidePair<int>& itran;
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