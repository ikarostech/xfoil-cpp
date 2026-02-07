#include "simulation/InviscidSolver.hpp"

#include "XFoil.h"
#include "infrastructure/logger.hpp"

#include <cmath>

using Eigen::Matrix2Xd;
using Eigen::VectorXd;

bool InviscidSolver::specConverge(XFoil &xfoil, SpecTarget target) {
  xfoil.surface_vortex = MathUtil::getRotateMatrix(xfoil.analysis_state_.alpha) *
                         xfoil.aerodynamicCache.gamu;

  if (target == SpecTarget::AngleOfAttack) {
    xfoil.updateTrailingEdgeState();
    xfoil.qinv_matrix =
        qiset(xfoil.analysis_state_.alpha, xfoil.aerodynamicCache.qinvu);
  } else {
    xfoil.minf_cl = xfoil.getActualMach(xfoil.analysis_state_.clspec,
                                        xfoil.analysis_state_.machType);
    xfoil.reinf_cl = xfoil.getActualReynolds(xfoil.analysis_state_.clspec,
                                             xfoil.analysis_state_.reynoldsType);
  }

  xfoil.applyClComputation(xfoil.clcalc(xfoil.cmref));

  bool bConv = false;

  if (target == SpecTarget::AngleOfAttack) {
    double clm = 1.0;
    double minf_clm = xfoil.getActualMach(clm, xfoil.analysis_state_.machType);

    for (int itcl = 1; itcl <= 20; itcl++) {
      const double msq_clm = 2.0 * xfoil.analysis_state_.currentMach * minf_clm;
      const double dclm = (xfoil.aero_coeffs_.cl - clm) /
                          (1.0 - xfoil.aero_coeffs_.cl_msq * msq_clm);

      const double clm1 = clm;
      double rlx = 1.0;

      //------ under-relaxation loop to avoid driving m(cl) above 1
      for (int irlx = 1; irlx <= 12; irlx++) {
        clm = clm1 + rlx * dclm;

        //-------- set new freestream mach m(clm)
        minf_clm = xfoil.getActualMach(clm, xfoil.analysis_state_.machType);

        //-------- if mach is ok, go do next newton iteration
        // FIXME double型の==比較
        if (xfoil.analysis_state_.machType == XFoil::MachType::CONSTANT ||
            xfoil.analysis_state_.currentMach == 0.0 ||
            minf_clm != 0.0)
          break;

        rlx *= 0.5;
      }

      //------ set new cl(m)
      xfoil.applyClComputation(xfoil.clcalc(xfoil.cmref));

      if (fabs(dclm) <= 1.0e-6) {
        bConv = true;
        break;
      }
    }

    if (!bConv) {
      Logger::instance().write("Specal:  MInf convergence failed\n");
      return false;
    }

    //---- set final mach, cl, cp distributions, and hinge moment
    xfoil.qinv_matrix =
        qiset(xfoil.analysis_state_.alpha, xfoil.aerodynamicCache.qinvu);
    xfoil.applyClComputation(xfoil.clcalc(xfoil.cmref));

    xfoil.cpi = cpcalc(xfoil.foil.foil_shape.n,
                             xfoil.qinv_matrix.row(0).transpose(),
                             xfoil.analysis_state_.qinf,
                             xfoil.analysis_state_.currentMach);
    if (xfoil.analysis_state_.viscous) {
      xfoil.cpv =
          cpcalc(xfoil.foil.foil_shape.n + xfoil.foil.wake_shape.n,
                       xfoil.qvis, xfoil.analysis_state_.qinf,
                       xfoil.analysis_state_.currentMach);
      xfoil.cpi =
          cpcalc(xfoil.foil.foil_shape.n + xfoil.foil.wake_shape.n,
                       xfoil.qinv_matrix.row(0).transpose(),
                       xfoil.analysis_state_.qinf,
                       xfoil.analysis_state_.currentMach);
    } else
      xfoil.cpi = cpcalc(xfoil.foil.foil_shape.n,
                               xfoil.qinv_matrix.row(0).transpose(),
                               xfoil.analysis_state_.qinf,
                               xfoil.analysis_state_.currentMach);

    for (int i = 0; i < xfoil.foil.foil_shape.n; i++) {
      xfoil.qgamm[i] = xfoil.surface_vortex(0, i);
    }

    return true;
  }

  for (int ital = 1; ital <= 20; ital++) {
    const double dalfa =
        (xfoil.analysis_state_.clspec - xfoil.aero_coeffs_.cl) /
        xfoil.aero_coeffs_.cl_alf;
    double rlx = 1.0;

    xfoil.analysis_state_.alpha =
        xfoil.analysis_state_.alpha + rlx * dalfa;

    //------ set new cl(alpha)
    xfoil.applyClComputation(xfoil.clcalc(xfoil.cmref));

    if (fabs(dalfa) <= 1.0e-6) {
      bConv = true;
      break;
    }
  }
  if (!bConv) {
    Logger::instance().write("Speccl:  cl convergence failed");
    return false;
  }

  //---- set final surface speed and cp distributions
  xfoil.updateTrailingEdgeState();
  xfoil.qinv_matrix =
      qiset(xfoil.analysis_state_.alpha, xfoil.aerodynamicCache.qinvu);

  if (xfoil.analysis_state_.viscous) {
    xfoil.cpv =
        cpcalc(xfoil.foil.foil_shape.n + xfoil.foil.wake_shape.n,
                     xfoil.qvis, xfoil.analysis_state_.qinf,
                     xfoil.analysis_state_.currentMach);
    xfoil.cpi =
        cpcalc(xfoil.foil.foil_shape.n + xfoil.foil.wake_shape.n,
                     xfoil.qinv_matrix.row(0).transpose(),
                     xfoil.analysis_state_.qinf,
                     xfoil.analysis_state_.currentMach);

  } else {
    xfoil.cpi = cpcalc(xfoil.foil.foil_shape.n,
                             xfoil.qinv_matrix.row(0).transpose(),
                             xfoil.analysis_state_.qinf,
                             xfoil.analysis_state_.currentMach);
  }

  return true;
}

bool InviscidSolver::specal(XFoil& xfoil) {
  return specConverge(xfoil, SpecTarget::AngleOfAttack);
}

bool InviscidSolver::speccl(XFoil& xfoil) {
  return specConverge(xfoil, SpecTarget::LiftCoefficient);
}

Matrix2Xd InviscidSolver::qiset(double alpha, const Matrix2Xd& qinvu) {
  return MathUtil::getRotateMatrix(alpha) * qinvu;
}

VectorXd InviscidSolver::cpcalc(int n, VectorXd q, double qinf, double minf) {
  VectorXd cp = VectorXd::Zero(n);
  bool denneg = false;
  const double beta = sqrt(1.0 - MathUtil::pow(minf, 2));
  const double prandtlGlauertFactor =
      0.5 * MathUtil::pow(minf, 2) / (1.0 + beta);

  for (int i = 0; i < n; i++) {
    const double cpinc = 1.0 - (q[i] / qinf) * (q[i] / qinf);
    const double den = beta + prandtlGlauertFactor * cpinc;
    cp[i] = cpinc / den;
    if (den <= 0.0)
      denneg = true;
  }

  if (denneg) {
    Logger::instance().write(
        "CpCalc: local speed too larger\n Compressibility corrections "
        "invalid\n");
  }

  return cp;
}
