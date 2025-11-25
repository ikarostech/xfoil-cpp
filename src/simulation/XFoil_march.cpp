#include "XFoil.h"
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "domain/boundary_layer/boundary_layer_builder.hpp"
using Eigen::Matrix;
using Eigen::Vector;
using Eigen::Vector2d;
using Eigen::VectorXd;

SetblInputView SetblInputView::fromXFoil(const XFoil& xfoil) {
  const auto& lattice = xfoil.boundaryLayerWorkflow.lattice;
  return SetblInputView{
      xfoil.lblini,
      {lattice.top.profiles.edgeVelocity, lattice.bottom.profiles.edgeVelocity},
      {lattice.top.profiles.skinFrictionCoeff, lattice.bottom.profiles.skinFrictionCoeff},
      {lattice.top.profiles.momentumThickness, lattice.bottom.profiles.momentumThickness},
      {lattice.top.profiles.displacementThickness, lattice.bottom.profiles.displacementThickness},
      {lattice.top.profiles.massFlux, lattice.bottom.profiles.massFlux},
      {lattice.top.skinFrictionCoeffHistory, lattice.bottom.skinFrictionCoeffHistory},
      {lattice.top.transitionIndex, lattice.bottom.transitionIndex}};
}

SetblOutputView SetblOutputView::fromXFoil(XFoil& xfoil) {
  auto& lattice = xfoil.boundaryLayerWorkflow.lattice;
  return SetblOutputView{
      xfoil.lblini,
      xfoil.gm1bl,
      xfoil.qinfbl,
      xfoil.tkbl,
      xfoil.tkbl_ms,
      xfoil.rstbl,
      xfoil.rstbl_ms,
      xfoil.hstinv,
      xfoil.hstinv_ms,
      xfoil.reybl,
      xfoil.reybl_re,
      xfoil.reybl_ms,
      xfoil.amcrit,
      {lattice.top.profiles.edgeVelocity, lattice.bottom.profiles.edgeVelocity},
      {lattice.top.profiles.skinFrictionCoeff, lattice.bottom.profiles.skinFrictionCoeff},
      {lattice.top.profiles.momentumThickness, lattice.bottom.profiles.momentumThickness},
      {lattice.top.profiles.displacementThickness, lattice.bottom.profiles.displacementThickness},
      {lattice.top.profiles.massFlux, lattice.bottom.profiles.massFlux},
      {lattice.top.skinFrictionCoeffHistory, lattice.bottom.skinFrictionCoeffHistory},
      {lattice.top.transitionIndex, lattice.bottom.transitionIndex},
      xfoil.va,
      xfoil.vb,
      xfoil.vdel,
      xfoil.vm,
      xfoil.vz,
      xfoil.tran,
      xfoil.turb,
      xfoil.wake,
      xfoil.simi,
      xfoil.xiforc};
}

XFoil::MixedModeStationContext XFoil::prepareMixedModeStation(int side, int ibl,
                                                              int itrold,
                                                              double& ami) {
  MixedModeStationContext ctx;

  ctx.simi = (ibl == 0);
  ctx.wake = ibl > boundaryLayerWorkflow.lattice.get(side).trailingEdgeIndex;
  ctx.xsi = boundaryLayerWorkflow.lattice.get(side).arcLengthCoordinates[ibl];
  ctx.uei = boundaryLayerWorkflow.lattice.get(side).profiles.edgeVelocity[ibl];
  ctx.thi = boundaryLayerWorkflow.lattice.get(side).profiles.momentumThickness[ibl];
  ctx.dsi = boundaryLayerWorkflow.lattice.get(side).profiles.displacementThickness[ibl];

  if (ibl < itrold) {
    ami = boundaryLayerWorkflow.lattice.get(side).profiles.skinFrictionCoeff[ibl];
    ctx.cti = 0.03;
  } else {
    ctx.cti = boundaryLayerWorkflow.lattice.get(side).profiles.skinFrictionCoeff[ibl];
    if (ctx.cti <= 0.0) {
      ctx.cti = 0.03;
    }
  }
  ctx.ami = ami;

  if (ctx.wake) {
    int iw = ibl - boundaryLayerWorkflow.lattice.get(side).trailingEdgeIndex;
    ctx.dswaki = wgap[iw - 1];
  } else {
    ctx.dswaki = 0.0;
  }

  double thickness_limit = (ibl <= boundaryLayerWorkflow.lattice.get(side).trailingEdgeIndex) ? 1.02 : 1.00005;
  ctx.dsi = std::max(ctx.dsi - ctx.dswaki, thickness_limit * ctx.thi) + ctx.dswaki;

  simi = ctx.simi;
  wake = ctx.wake;

  return ctx;
}

void XFoil::checkTransitionIfNeeded(int side, int ibl, bool skipCheck,
                                    int laminarAdvance, double& ami) {
  if (skipCheck || turb) {
    return;
  }

  trchek();
  ami = boundaryLayerWorkflow.state.station2.param.amplz;
  if (tran) {
    boundaryLayerWorkflow.lattice.get(side).transitionIndex = ibl;
  } else {
    boundaryLayerWorkflow.lattice.get(side).transitionIndex = ibl + laminarAdvance;
  }
}

bool XFoil::performMixedModeNewtonIteration(int side, int ibl, int itrold,
                                            MixedModeStationContext& ctx,
                                            double deps, double senswt,
                                            double& sens, double& sennew,
                                            double& ami) {
  bool converged = false;
  double ueref = 0.0;
  double hkref = 0.0;

  for (int itbl = 1; itbl <= 25; ++itbl) {
    {
  blData updatedCurrent =
      boundaryLayerWorkflow.blprv(*this, boundaryLayerWorkflow.state.current(),
                                  ctx.xsi, ami, ctx.cti, ctx.thi, ctx.dsi,
                                  ctx.dswaki, ctx.uei);
      boundaryLayerWorkflow.state.current() = updatedCurrent;
    }
    blkin(boundaryLayerWorkflow.state);

    checkTransitionIfNeeded(side, ibl, ctx.simi, 1, ami);

    const bool startOfWake =
        boundaryLayerWorkflow.isStartOfWake(*this, side, ibl);
    boundaryLayerWorkflow.updateSystemMatricesForStation(*this, side, ibl,
                                                          ctx);

    if (itbl == 1) {
      boundaryLayerWorkflow.initializeFirstIterationState(
          *this, side, ibl, itrold, ctx, ueref, hkref, ami);
    }

    if (ctx.simi || startOfWake) {
      boundaryLayerWorkflow.configureSimilarityRow(ueref);
    } else {
      const bool resetSensitivity = (itbl <= 5);
      const bool averageSensitivity = (itbl > 5 && itbl <= 15);
      boundaryLayerWorkflow.configureViscousRow(
          hkref, ueref, senswt, resetSensitivity, averageSensitivity,
          sens, sennew);
    }

    if (boundaryLayerWorkflow.applyMixedModeNewtonStep(*this, side, ibl, deps,
                                                       ami, ctx)) {
      converged = true;
      break;
    }
  }

  return converged;
}

void XFoil::handleMixedModeNonConvergence(int side, int ibl,
                                          MixedModeStationContext& ctx,
                                          double& ami) {
  std::stringstream ss;
  ss << "     mrchdu: convergence failed at " << ibl << " ,  side " << side
     << ", res=" << std::setw(4) << std::fixed << std::setprecision(3)
     << ctx.dmax << "\n";
  writeString(ss.str());

  boundaryLayerWorkflow.resetStationKinematicsAfterFailure(
      side, ibl, ctx,
      BoundaryLayerWorkflow::EdgeVelocityFallbackMode::UsePreviousStation);

  {
    blData updatedCurrent =
        boundaryLayerWorkflow.blprv(
            *this, boundaryLayerWorkflow.state.current(), ctx.xsi, ami,
            ctx.cti, ctx.thi, ctx.dsi, ctx.dswaki, ctx.uei);
    boundaryLayerWorkflow.state.current() = updatedCurrent;
  }
  blkin(boundaryLayerWorkflow.state);

  checkTransitionIfNeeded(side, ibl, ctx.simi, 2, ami);

  boundaryLayerWorkflow.syncStationRegimeStates(side, ibl, ctx.wake, *this);

  ctx.ami = ami;
}

double XFoil::calcHtarg(int ibl, int is, bool wake) {
  if (ibl < boundaryLayerWorkflow.lattice.get(is).transitionIndex) {
    return boundaryLayerWorkflow.state.station1.hkz.scalar +
           0.03 * (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz) / boundaryLayerWorkflow.state.station1.param.tz;
  } else if (ibl == boundaryLayerWorkflow.lattice.get(is).transitionIndex) {
    return boundaryLayerWorkflow.state.station1.hkz.scalar +
           (0.03 * (xt - boundaryLayerWorkflow.state.station1.param.xz) - 0.15 * (boundaryLayerWorkflow.state.station2.param.xz - xt)) /
               boundaryLayerWorkflow.state.station1.param.tz;
  } else if (wake) {
    const double cst =
        0.03 * (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz) / boundaryLayerWorkflow.state.station1.param.tz;
    auto euler = [](double hk2, double hk1, double cst) {
      return hk2 - (hk2 + cst * pow(hk2 - 1, 3) - hk1) /
                       (1 + 3 * cst * pow(hk2 - 1, 2));
    };
    boundaryLayerWorkflow.state.station2.hkz.scalar = boundaryLayerWorkflow.state.station1.hkz.scalar;
    for (int i = 0; i < 3; i++) {
      boundaryLayerWorkflow.state.station2.hkz.scalar = euler(boundaryLayerWorkflow.state.station2.hkz.scalar, boundaryLayerWorkflow.state.station1.hkz.scalar, cst);
    }
    return boundaryLayerWorkflow.state.station2.hkz.scalar;
  } else {
    return boundaryLayerWorkflow.state.station1.hkz.scalar -
           0.15 * (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz) / boundaryLayerWorkflow.state.station1.param.tz;
  }
}

SetblOutputView XFoil::setbl(const SetblInputView& input,
                                    SetblOutputView output) {
  //-------------------------------------------------
  //	   sets up the bl newton system coefficients for the current bl
  // variables
  //     and the edge velocities received from setup. the local bl system
  //     coefficients are then incorporated into the global newton system.
  //-------------------------------------------------

  std::stringstream ss;
  int jvte1 = 0, jvte2 = 0;
  VectorXd u1_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd u2_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd d1_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd d2_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd ule1_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd ule2_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd ute1_m = VectorXd::Zero(2 * IVX + 1);
  VectorXd ute2_m = VectorXd::Zero(2 * IVX + 1);

  double msq_clmr = 0.0, mdi;
  double herat = 0.0, herat_ms = 0.0;

  double clmr = 0.0, ma_clmr = 0.0, re_clmr = 0.0;
  double ule1_a = 0.0, ule2_a = 0.0, u2_a, due2, dds2;
  double xsi, cti = 0.0, uei, thi, dsi, dswaki;
  double d2_a, d2_m2, d2_u2, dte_mte1, dte_ute1, dte_mte2, dte_ute2;
  double tte, cte, dte, dule1 = 0.0, dule2 = 0.0;
  double xi_ule1, xi_ule2;
  double ami = 0.0, tte_tte1 = 0.0, tte_tte2 = 0.0, cte_tte1 = 0.0,
         cte_tte2 = 0.0, cte_cte1 = 0.0, cte_cte2 = 0.0;

  //---- set the cl used to define mach, reynolds numbers
  if (analysis_state_.controlByAlpha)
    clmr = cl;
  else
    clmr = analysis_state_.clspec;

  cti = 0.0; // techwinder added, otherwise variable is not initialized

  //---- set current minf(cl)
  ma_clmr = getActualMach(clmr, analysis_state_.machType);
  re_clmr = getActualReynolds(clmr, analysis_state_.reynoldsType);

  msq_clmr = 2.0 * analysis_state_.currentMach * ma_clmr;

  //---- set compressibility parameter tklam and derivative tk_msq
  const auto compressibility = buildCompressibilityParams();
  tklam = compressibility.karmanTsienFactor;
  tkl_msq = compressibility.karmanTsienFactor_msq;

  //---- set gas constant (= cp/cv)
  output.gm1bl = gamm1;

  //---- set parameters for compressibility correction
  output.qinfbl = analysis_state_.qinf;
  output.tkbl = tklam;
  output.tkbl_ms = tkl_msq;

  //---- stagnation density and 1/enthalpy
  output.rstbl =
      pow((1.0 + 0.5 * output.gm1bl * analysis_state_.currentMach * analysis_state_.currentMach),
          (1.0 / output.gm1bl));
  output.rstbl_ms =
      0.5 * output.rstbl /
      (1.0 + 0.5 * output.gm1bl * analysis_state_.currentMach * analysis_state_.currentMach);
  output.hstinv = output.gm1bl *
                  MathUtil::pow(analysis_state_.currentMach / output.qinfbl, 2) /
                  (1.0 + 0.5 * output.gm1bl * analysis_state_.currentMach * analysis_state_.currentMach);
  output.hstinv_ms =
      output.gm1bl * MathUtil::pow(1.0 / output.qinfbl, 2) /
          (1.0 + 0.5 * output.gm1bl * analysis_state_.currentMach * analysis_state_.currentMach) -
      0.5 * output.gm1bl * output.hstinv /
          (1.0 + 0.5 * output.gm1bl * analysis_state_.currentMach * analysis_state_.currentMach);

  //---- set reynolds number based on freestream density, velocity, viscosity
  herat = 1.0 - 0.5 * output.qinfbl * output.qinfbl * output.hstinv;
  herat_ms = -0.5 * output.qinfbl * output.qinfbl * output.hstinv_ms;

  output.reybl = analysis_state_.currentRe * sqrt(herat * herat * herat) *
                 (1.0 + hvrat) / (herat + hvrat);
  output.reybl_re = sqrt(herat * herat * herat) * (1.0 + hvrat) / (herat + hvrat);
  output.reybl_ms = output.reybl * (1.5 / herat - 1.0 / (herat + hvrat)) * herat_ms;

  output.amcrit = acrit;

  if (!input.lblini) {
    //----- initialize bl by marching with ue (fudge at separation)
    // TRACE(" initializing bl ...\n");
    writeString("   Initializing bl ...\n");

    boundaryLayerWorkflow.mrchue(*this);
    output.lblini = true;
  }

  //---- march bl with current ue and ds to establish transition
  boundaryLayerWorkflow.mrchdu(*this);

  SidePair<VectorXd> usav;
  usav.top = input.edgeVelocity.top;
  usav.bottom = input.edgeVelocity.bottom;

  ueset();
  const auto swapped_edge_velocities = swapEdgeVelocities(usav);
  usav = swapped_edge_velocities.swappedUsav;
  output.edgeVelocity.top = swapped_edge_velocities.restoredUedg.top;
  output.edgeVelocity.bottom = swapped_edge_velocities.restoredUedg.bottom;
  jvte1 = boundaryLayerWorkflow.lattice.top.stationToSystem[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex];
  jvte2 = boundaryLayerWorkflow.lattice.bottom.stationToSystem[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex];

  dule1 = output.edgeVelocity.top[0] - usav.top[0];
  dule2 = output.edgeVelocity.bottom[0] - usav.bottom[0];

  //---- set le and te ue sensitivities wrt all m values
  const auto le_te_sensitivities = computeLeTeSensitivities(
      boundaryLayerWorkflow.lattice.get(1).stationToPanel[0], boundaryLayerWorkflow.lattice.get(2).stationToPanel[0], boundaryLayerWorkflow.lattice.get(1).stationToPanel[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex],
      boundaryLayerWorkflow.lattice.get(2).stationToPanel[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex]);
  ule1_m = le_te_sensitivities.ule1_m;
  ule2_m = le_te_sensitivities.ule2_m;
  ute1_m = le_te_sensitivities.ute1_m;
  ute2_m = le_te_sensitivities.ute2_m;

  ule1_a = boundaryLayerWorkflow.lattice.get(1).inviscidEdgeVelocityDerivative[0];
  ule2_a = boundaryLayerWorkflow.lattice.get(2).inviscidEdgeVelocityDerivative[0];

  writeString(" \n");

  //*** process each boundary layer side
  for (int is = 1; is <= 2; is++) {
    //---- there is no station "1" at similarity, so zero everything out
    const auto cleared_derivatives = clearDerivativeVectors(u1_m, d1_m);
    u1_m = cleared_derivatives.u;
    d1_m = cleared_derivatives.d;
    double u1_a = 0.0;
    double d1_a = 0.0;

    double due1 = 0.0;
    double dds1 = 0.0;

    //---- set forced transition arc length position
    output.xiforc = xifset(is);

    output.tran = false;
    output.turb = false;

    //**** sweep downstream setting up bl equation linearizations
    for (int ibl = 0; ibl < boundaryLayerWorkflow.lattice.get(is).stationCount - 1; ++ibl) {
      
      int iv = boundaryLayerWorkflow.lattice.get(is).stationToSystem[ibl];

      output.simi = (ibl == 0);
      output.wake = (ibl > boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex);
      output.tran = (ibl == output.itran.get(is));
      output.turb = (ibl > output.itran.get(is));

      //---- set primary variables for current station
      xsi = boundaryLayerWorkflow.lattice.get(is).arcLengthCoordinates[ibl];
      if (ibl < output.itran.get(is))
        ami = output.skinFrictionCoeff.get(is)[ibl];
      else
        cti = output.skinFrictionCoeff.get(is)[ibl];
      uei = output.edgeVelocity.get(is)[ibl];
      thi = output.momentumThickness.get(is)[ibl];
      mdi = output.massFlux.get(is)[ibl];

      dsi = mdi / uei;

      if (output.wake) {
        int iw = ibl - boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex;
        dswaki = wgap[iw - 1];
      } else
        dswaki = 0.0;

      //---- set derivatives of dsi (= d2)
      d2_m2 = 1.0 / uei;
      d2_u2 = -dsi / uei;

      for (int js = 1; js <= 2; js++) {
        for (int jbl = 0; jbl < boundaryLayerWorkflow.lattice.get(js).stationCount - 1; ++jbl) {
          int jv = boundaryLayerWorkflow.lattice.get(js).stationToSystem[jbl];
          u2_m[jv] = -boundaryLayerWorkflow.lattice.get(is).panelInfluenceFactor[ibl] * boundaryLayerWorkflow.lattice.get(js).panelInfluenceFactor[jbl] *
                     dij(boundaryLayerWorkflow.lattice.get(is).stationToPanel[ibl], boundaryLayerWorkflow.lattice.get(js).stationToPanel[jbl]);
          d2_m[jv] = d2_u2 * u2_m[jv];
        }
      }
      d2_m[iv] = d2_m[iv] + d2_m2;

      u2_a = boundaryLayerWorkflow.lattice.get(is).inviscidEdgeVelocityDerivative[ibl];
      d2_a = d2_u2 * u2_a;

  //---- "forced" changes due to mismatch between edge velocities and
  // usav=inviscidEdgeVelocity+dij*output.massFlux
      due2 = output.edgeVelocity.get(is)[ibl] - usav.get(is)[ibl];
      dds2 = d2_u2 * due2;

  {
    blData updatedCurrent =
        boundaryLayerWorkflow.blprv(*this,
                                    boundaryLayerWorkflow.state.current(), xsi,
                                    ami, cti, thi, dsi, dswaki, uei);
    boundaryLayerWorkflow.state.current() = updatedCurrent;
  } // cti
      blkin(boundaryLayerWorkflow.state);

      //---- check for transition and set output.tran, xt, etc. if found
      if (output.tran) {
        trchek();
        ami = boundaryLayerWorkflow.state.station2.param.amplz;
      }

      if (ibl == output.itran.get(is) && !output.tran) {
        // TRACE("setbl: xtr???  n1=%d n2=%d: \n", ampl1, ampl2);

        ss << "setbl: xtr???  n1=" << boundaryLayerWorkflow.state.station1.param.amplz
           << " n2=" << boundaryLayerWorkflow.state.station2.param.amplz << ":\n";
        writeString(ss.str());
        ss.str("");
      }

      //---- assemble 10x4 linearized system for dskinFrictionCoeff, dth, dds, due, dxi
      //	   at the previous "1" station and the current "2" station

      if (ibl == boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex + 1) {
        //----- define quantities at start of output.wake, adding te base thickness to
        // dstar
        tte = output.momentumThickness.get(1)[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] +
              output.momentumThickness.get(2)[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex];
        dte = output.displacementThickness.get(1)[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] +
              output.displacementThickness.get(2)[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex] + foil.edge.ante;
        cte = (output.skinFrictionCoeff.get(1)[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] *
                   output.momentumThickness.get(1)[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] +
               output.skinFrictionCoeff.get(2)[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex] *
                   output.momentumThickness.get(2)[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex]) /
               tte;
        boundaryLayerWorkflow.tesys(*this, cte, tte, dte);

        tte_tte1 = 1.0;
        tte_tte2 = 1.0;
        dte_mte1 = 1.0 / output.edgeVelocity.top[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex];
        dte_ute1 = -output.displacementThickness.get(1)[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] /
                    output.edgeVelocity.top[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex];
        dte_mte2 = 1.0 / output.edgeVelocity.bottom[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex];
        dte_ute2 = -output.displacementThickness.get(2)[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex] /
                    output.edgeVelocity.bottom[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex];
        cte_cte1 = output.momentumThickness.get(1)[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] / tte;
        cte_cte2 = output.momentumThickness.get(2)[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex] / tte;
        cte_tte1 = (output.skinFrictionCoeff.get(1)[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] - cte) / tte;
        cte_tte2 = (output.skinFrictionCoeff.get(2)[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex] - cte) / tte;

        //----- re-define d1 sensitivities wrt m since d1 depends on both te ds
        // values
      for (int js = 1; js <= 2; js++) {
        for (int jbl = 0; jbl < boundaryLayerWorkflow.lattice.get(js).stationCount - 1; ++jbl) {
            int jv = boundaryLayerWorkflow.lattice.get(js).stationToSystem[jbl];
            d1_m[jv] = dte_ute1 * ute1_m[jv] + dte_ute2 * ute2_m[jv];
          }
        }
        d1_m[jvte1] = d1_m[jvte1] + dte_mte1;
        d1_m[jvte2] = d1_m[jvte2] + dte_mte2;

        //----- "forced" changes from  output.edgeVelocity --- usav=inviscidEdgeVelocity+dij*output.massFlux	mismatch
        due1 = 0.0;
        dds1 =
            dte_ute1 *
                (output.edgeVelocity.top[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex] -
                 usav.top[boundaryLayerWorkflow.lattice.top.trailingEdgeIndex]) +
            dte_ute2 *
                (output.edgeVelocity.bottom[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex] -
                 usav.bottom[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex]);
      } else {
        boundaryLayerWorkflow.blsys(*this);
      }

      //---- save wall shear and equil. max shear coefficient for plotting
      // output
      output.skinFrictionCoeffHistory.get(is)[ibl] = boundaryLayerWorkflow.state.station2.cqz.scalar;

      //---- set xi sensitivities wrt le ue changes
      if (is == 1) {
        xi_ule1 = sst_go;
        xi_ule2 = -sst_gp;
      } else {
        xi_ule1 = -sst_go;
        xi_ule2 = sst_gp;
      }

      //---- stuff bl system coefficients into main jacobian matrix

      for (int jv = 1; jv <= nsys; jv++) {
        output.vm[0][jv][iv] = boundaryLayerWorkflow.blc.a1(0, 2) * d1_m[jv] + boundaryLayerWorkflow.blc.a1(0, 3) * u1_m[jv] +
                        boundaryLayerWorkflow.blc.a2(0, 2) * d2_m[jv] + boundaryLayerWorkflow.blc.a2(0, 3) * u2_m[jv] +
                        (boundaryLayerWorkflow.blc.a1(0, 4) + boundaryLayerWorkflow.blc.a2(0, 4) + boundaryLayerWorkflow.blc.d_xi[0]) *
                            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
      }

      output.vb[iv](0, 0) = boundaryLayerWorkflow.blc.a1(0, 0);
      output.vb[iv](0, 1) = boundaryLayerWorkflow.blc.a1(0, 1);

      output.va[iv](0, 0) = boundaryLayerWorkflow.blc.a2(0, 0);
      output.va[iv](0, 1) = boundaryLayerWorkflow.blc.a2(0, 1);

      if (analysis_state_.controlByAlpha)
        output.vdel[iv](0, 1) = boundaryLayerWorkflow.blc.d_re[0] * re_clmr + boundaryLayerWorkflow.blc.d_msq[0] * msq_clmr;
      else
        output.vdel[iv](0, 1) = (boundaryLayerWorkflow.blc.a1(0, 3) * u1_a + boundaryLayerWorkflow.blc.a1(0, 2) * d1_a) +
                         (boundaryLayerWorkflow.blc.a2(0, 3) * u2_a + boundaryLayerWorkflow.blc.a2(0, 2) * d2_a) +
                         (boundaryLayerWorkflow.blc.a1(0, 4) + boundaryLayerWorkflow.blc.a2(0, 4) + boundaryLayerWorkflow.blc.d_xi[0]) *
                             (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

      output.vdel[iv](0, 0) = boundaryLayerWorkflow.blc.rhs[0] + (boundaryLayerWorkflow.blc.a1(0, 3) * due1 + boundaryLayerWorkflow.blc.a1(0, 2) * dds1) +
                       (boundaryLayerWorkflow.blc.a2(0, 3) * due2 + boundaryLayerWorkflow.blc.a2(0, 2) * dds2) +
                       (boundaryLayerWorkflow.blc.a1(0, 4) + boundaryLayerWorkflow.blc.a2(0, 4) + boundaryLayerWorkflow.blc.d_xi[0]) *
                           (xi_ule1 * dule1 + xi_ule2 * dule2);

      for (int jv = 1; jv <= nsys; jv++) {
        output.vm[1][jv][iv] = boundaryLayerWorkflow.blc.a1(1, 2) * d1_m[jv] + boundaryLayerWorkflow.blc.a1(1, 3) * u1_m[jv] +
                        boundaryLayerWorkflow.blc.a2(1, 2) * d2_m[jv] + boundaryLayerWorkflow.blc.a2(1, 3) * u2_m[jv] +
                        (boundaryLayerWorkflow.blc.a1(1, 4) + boundaryLayerWorkflow.blc.a2(1, 4) + boundaryLayerWorkflow.blc.d_xi[1]) *
                            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
      }
      output.vb[iv](1, 0) = boundaryLayerWorkflow.blc.a1(1, 0);
      output.vb[iv](1, 1) = boundaryLayerWorkflow.blc.a1(1, 1);

      output.va[iv](1, 0) = boundaryLayerWorkflow.blc.a2(1, 0);
      output.va[iv](1, 1) = boundaryLayerWorkflow.blc.a2(1, 1);

      if (analysis_state_.controlByAlpha)
        output.vdel[iv](1, 1) = boundaryLayerWorkflow.blc.d_re[1] * re_clmr + boundaryLayerWorkflow.blc.d_msq[1] * msq_clmr;
      else
        output.vdel[iv](1, 1) = (boundaryLayerWorkflow.blc.a1(1, 3) * u1_a + boundaryLayerWorkflow.blc.a1(1, 2) * d1_a) +
                         (boundaryLayerWorkflow.blc.a2(1, 3) * u2_a + boundaryLayerWorkflow.blc.a2(1, 2) * d2_a) +
                         (boundaryLayerWorkflow.blc.a1(1, 4) + boundaryLayerWorkflow.blc.a2(1, 4) + boundaryLayerWorkflow.blc.d_xi[1]) *
                             (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

      output.vdel[iv](1, 0) = boundaryLayerWorkflow.blc.rhs[1] + (boundaryLayerWorkflow.blc.a1(1, 3) * due1 + boundaryLayerWorkflow.blc.a1(1, 2) * dds1) +
                       (boundaryLayerWorkflow.blc.a2(1, 3) * due2 + boundaryLayerWorkflow.blc.a2(1, 2) * dds2) +
                       (boundaryLayerWorkflow.blc.a1(1, 4) + boundaryLayerWorkflow.blc.a2(1, 4) + boundaryLayerWorkflow.blc.d_xi[1]) *
                           (xi_ule1 * dule1 + xi_ule2 * dule2);

      // memory overlap problem
      for (int jv = 1; jv <= nsys; jv++) {
        output.vm[2][jv][iv] = boundaryLayerWorkflow.blc.a1(2, 2) * d1_m[jv] + boundaryLayerWorkflow.blc.a1(2, 3) * u1_m[jv] +
                        boundaryLayerWorkflow.blc.a2(2, 2) * d2_m[jv] + boundaryLayerWorkflow.blc.a2(2, 3) * u2_m[jv] +
                        (boundaryLayerWorkflow.blc.a1(2, 4) + boundaryLayerWorkflow.blc.a2(2, 4) + boundaryLayerWorkflow.blc.d_xi[2]) *
                            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
      }

      output.vb[iv](2, 0) = boundaryLayerWorkflow.blc.a1(2, 0);
      output.vb[iv](2, 1) = boundaryLayerWorkflow.blc.a1(2, 1);

      output.va[iv](2, 0) = boundaryLayerWorkflow.blc.a2(2, 0);
      output.va[iv](2, 1) = boundaryLayerWorkflow.blc.a2(2, 1);

      if (analysis_state_.controlByAlpha)
        output.vdel[iv](2, 1) = boundaryLayerWorkflow.blc.d_re[2] * re_clmr + boundaryLayerWorkflow.blc.d_msq[2] * msq_clmr;
      else
        output.vdel[iv](2, 1) = (boundaryLayerWorkflow.blc.a1(2, 3) * u1_a + boundaryLayerWorkflow.blc.a1(2, 2) * d1_a) +
                         (boundaryLayerWorkflow.blc.a2(2, 3) * u2_a + boundaryLayerWorkflow.blc.a2(2, 2) * d2_a) +
                         (boundaryLayerWorkflow.blc.a1(2, 4) + boundaryLayerWorkflow.blc.a2(2, 4) + boundaryLayerWorkflow.blc.d_xi[2]) *
                             (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

      output.vdel[iv](2, 0) = boundaryLayerWorkflow.blc.rhs[2] + (boundaryLayerWorkflow.blc.a1(2, 3) * due1 + boundaryLayerWorkflow.blc.a1(2, 2) * dds1) +
                       (boundaryLayerWorkflow.blc.a2(2, 3) * due2 + boundaryLayerWorkflow.blc.a2(2, 2) * dds2) +
                       (boundaryLayerWorkflow.blc.a1(2, 4) + boundaryLayerWorkflow.blc.a2(2, 4) + boundaryLayerWorkflow.blc.d_xi[2]) *
                           (xi_ule1 * dule1 + xi_ule2 * dule2);

      if (ibl == boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex + 1) {
        //----- redefine coefficients for tte, dte, etc
        output.vz[0][0] = boundaryLayerWorkflow.blc.a1(0, 0) * cte_cte1;
        output.vz[0][1] = boundaryLayerWorkflow.blc.a1(0, 0) * cte_tte1 + boundaryLayerWorkflow.blc.a1(0, 1) * tte_tte1;
        output.vb[iv](0, 0) = boundaryLayerWorkflow.blc.a1(0, 0) * cte_cte2;
        output.vb[iv](0, 1) = boundaryLayerWorkflow.blc.a1(0, 0) * cte_tte2 + boundaryLayerWorkflow.blc.a1(0, 1) * tte_tte2;

        output.vz[1][0] = boundaryLayerWorkflow.blc.a1(1, 0) * cte_cte1;
        output.vz[1][1] = boundaryLayerWorkflow.blc.a1(1, 0) * cte_tte1 + boundaryLayerWorkflow.blc.a1(1, 1) * tte_tte1;
        output.vb[iv](1, 0) = boundaryLayerWorkflow.blc.a1(1, 0) * cte_cte2;
        output.vb[iv](1, 1) = boundaryLayerWorkflow.blc.a1(1, 0) * cte_tte2 + boundaryLayerWorkflow.blc.a1(1, 1) * tte_tte2;

        output.vz[2][0] = boundaryLayerWorkflow.blc.a1(2, 0) * cte_cte1;
        output.vz[2][1] = boundaryLayerWorkflow.blc.a1(2, 0) * cte_tte1 + boundaryLayerWorkflow.blc.a1(2, 1) * tte_tte1;
        output.vb[iv](2, 0) = boundaryLayerWorkflow.blc.a1(2, 0) * cte_cte2;
        output.vb[iv](2, 1) = boundaryLayerWorkflow.blc.a1(2, 0) * cte_tte2 + boundaryLayerWorkflow.blc.a1(2, 1) * tte_tte2;
      }

      //---- turbulent intervals will follow if currently at transition interval
      if (output.tran) {
        output.turb = true;

        //------ save transition location
        output.itran.get(is) = ibl;
      }

      output.tran = false;

      if (ibl == boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex) {
        //----- set "2" variables at te to output.wake correlations for next station

        output.turb = true;
        output.wake = true;
        boundaryLayerWorkflow.state.station2 =
            boundaryLayerWorkflow.blvar(
                boundaryLayerWorkflow.state.station2,
                FlowRegimeEnum::Wake);
        boundaryLayerWorkflow.blmid(*this, FlowRegimeEnum::Wake);
      }
      u1_m = u2_m;
      d1_m = d2_m;

      u1_a = u2_a;
      d1_a = d2_a;

      due1 = due2;
      dds1 = dds2;

      //---- set bl variables for next station
      boundaryLayerWorkflow.state.stepbl();
    }
  }

  return output;
}

bool XFoil::trchek() {
  //----------------------------------------------------------------
  //     new second-order version:  december 1994.
  //
  //     checks if transition occurs in the current interval x1..x2.
  //     if transition occurs, then set transition location xt, and
  //     its sensitivities to "1" and "2" variables.  if no transition,
  //     set amplification ampl2.
  //
  //     solves the implicit amplification equation for n2:
  //
  //       n2 - n1     n'(xt,nt) + n'(x1,n1)
  //       -------  =  ---------------------
  //       x2 - x1               2
  //
  //     in effect, a 2-point central difference is used between
  //     x1..x2 (no transition), or x1..xt (transition).  the switch
  //     is done by defining xt,nt in the equation above depending
  //     on whether n2 exceeds ncrit.
  //
  //  if n2<ncrit:  nt=n2    , xt=x2                  (no transition)
  //
  //  if n2>ncrit:  nt=ncrit , xt=(ncrit-n1)/(n2-n1)  (transition)
  //
  //----------------------------------------------------------------
  double amplt, sfa, sfa_a1, sfa_a2, sfx;
  double sfx_x1, sfx_x2, sfx_xf;
  double tt, dt, ut, amsave;
  double res = 0.0, res_a2 = 0.0;
  double da2 = 0.0, dxt = 0.0, tt_t1 = 0.0, dt_d1 = 0.0, ut_u1 = 0.0;
  double tt_t2 = 0.0, dt_d2 = 0.0, ut_u2 = 0.0, tt_a1 = 0.0, dt_a1 = 0.0;
  double ut_a1 = 0.0, tt_x1 = 0.0, dt_x1 = 0.0, ut_x1 = 0.0, tt_x2 = 0.0,
         dt_x2 = 0.0, ut_x2 = 0.0;
  double z_ax = 0.0, z_a1 = 0.0, z_t1 = 0.0, z_d1 = 0.0, z_u1 = 0.0, z_x1 = 0.0,
         z_a2 = 0.0, z_t2 = 0.0, z_d2 = 0.0, z_u2 = 0.0, z_x2 = 0.0, z_ms = 0.0,
         z_re = 0.0;
  double amplt_a2, wf, wf_a1, wf_a2, wf_xf, wf_x1, wf_x2;
  double xt_a2, dt_a2, tt_a2;
  double ut_a2;
  double daeps = 0.00005;

  amplt_a2 = 0.0;
  xt_a2 = dt_a2 = tt_a2 = 0.0;
  ut_a2 = 0.0;

  //---- save variables and sensitivities at ibl ("2") for future restoration
  saveblData(2);

  //---- calculate average amplification rate ax over x1..x2 interval
  BoundaryLayerUtil::AxResult ax_result =
      BoundaryLayerUtil::axset(boundaryLayerWorkflow.state.station1.hkz.scalar, boundaryLayerWorkflow.state.station1.param.tz, boundaryLayerWorkflow.state.station1.rtz.scalar,
            boundaryLayerWorkflow.state.station1.param.amplz, boundaryLayerWorkflow.state.station2.hkz.scalar, boundaryLayerWorkflow.state.station2.param.tz,
            boundaryLayerWorkflow.state.station2.rtz.scalar, boundaryLayerWorkflow.state.station2.param.amplz, amcrit);

  //---- set initial guess for iterate n2 (ampl2) at x2
  boundaryLayerWorkflow.state.station2.param.amplz = boundaryLayerWorkflow.state.station1.param.amplz +
                        ax_result.ax * (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz);
  //---- solve implicit system for amplification ampl2
  auto iterateAmplification = [&]() -> bool {
    for (int itam = 0; itam < 30; itam++) {
      //---- define weighting factors wf1,wf2 for defining "t" quantities
      if (boundaryLayerWorkflow.state.station2.param.amplz <= amcrit) {
        //------ there is no transition yet,  "t" is the same as "2"
        amplt = boundaryLayerWorkflow.state.station2.param.amplz;
        amplt_a2 = 1.0;
        sfa = 1.0;
        sfa_a1 = 0.0;
        sfa_a2 = 0.0;
      } else {
        //------ there is transition in x1..x2, "t" is set from n1, n2
        amplt = amcrit;
        amplt_a2 = 0.0;
        sfa = (amplt - boundaryLayerWorkflow.state.station1.param.amplz) /
              (boundaryLayerWorkflow.state.station2.param.amplz - boundaryLayerWorkflow.state.station1.param.amplz);
        sfa_a1 = (sfa - 1.0) / (boundaryLayerWorkflow.state.station2.param.amplz - boundaryLayerWorkflow.state.station1.param.amplz);
        sfa_a2 = (-sfa) / (boundaryLayerWorkflow.state.station2.param.amplz - boundaryLayerWorkflow.state.station1.param.amplz);
      }

      if (xiforc < boundaryLayerWorkflow.state.station2.param.xz) {
        sfx = (xiforc - boundaryLayerWorkflow.state.station1.param.xz) /
              (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz);
        sfx_x1 = (sfx - 1.0) / (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz);
        sfx_x2 = (-sfx) / (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz);
        sfx_xf = 1.0 / (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz);
      } else {
        sfx = 1.0;
        sfx_x1 = 0.0;
        sfx_x2 = 0.0;
        sfx_xf = 0.0;
      }

      //---- set weighting factor from free or forced transition
      if (sfa < sfx) {
        wf = sfa;
        wf_a1 = sfa_a1;
        wf_a2 = sfa_a2;
        wf_x1 = 0.0;
        wf_x2 = 0.0;
        wf_xf = 0.0;
      } else {
        wf = sfx;
        wf_a1 = 0.0;
        wf_a2 = 0.0;
        wf_x1 = sfx_x1;
        wf_x2 = sfx_x2;
        wf_xf = sfx_xf;
      }

      //---- interpolate bl variables to xt
      xt = boundaryLayerWorkflow.state.station1.param.xz * (1 - wf) + boundaryLayerWorkflow.state.station2.param.xz * wf;
      tt = boundaryLayerWorkflow.state.station1.param.tz * (1 - wf) + boundaryLayerWorkflow.state.station2.param.tz * wf;
      dt = boundaryLayerWorkflow.state.station1.param.dz * (1 - wf) + boundaryLayerWorkflow.state.station2.param.dz * wf;
      ut = boundaryLayerWorkflow.state.station1.param.uz * (1 - wf) + boundaryLayerWorkflow.state.station2.param.uz * wf;

      xt_a2 = (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz) * wf_a2;
      tt_a2 = (boundaryLayerWorkflow.state.station2.param.tz - boundaryLayerWorkflow.state.station1.param.tz) * wf_a2;
      dt_a2 = (boundaryLayerWorkflow.state.station2.param.dz - boundaryLayerWorkflow.state.station1.param.dz) * wf_a2;
      ut_a2 = (boundaryLayerWorkflow.state.station2.param.uz - boundaryLayerWorkflow.state.station1.param.uz) * wf_a2;

      //---- temporarily set "2" variables from "t" for blkin
      boundaryLayerWorkflow.state.station2.param.xz = xt;
      boundaryLayerWorkflow.state.station2.param.tz = tt;
      boundaryLayerWorkflow.state.station2.param.dz = dt;
      boundaryLayerWorkflow.state.station2.param.uz = ut;

      //---- calculate laminar secondary "t" variables hkt, rtt
      blkin(boundaryLayerWorkflow.state);

      blData::blVector hkt = boundaryLayerWorkflow.state.station2.hkz;
      blData::blVector rtt = boundaryLayerWorkflow.state.station2.rtz;

      //---- restore clobbered "2" variables, except for ampl2
      amsave = boundaryLayerWorkflow.state.station2.param.amplz;

      restoreblData(2);

      boundaryLayerWorkflow.state.station2.param.amplz = amsave;

      //---- calculate amplification rate ax over current x1-xt interval
      ax_result = BoundaryLayerUtil::axset(boundaryLayerWorkflow.state.station1.hkz.scalar, boundaryLayerWorkflow.state.station1.param.tz,
                        boundaryLayerWorkflow.state.station1.rtz.scalar, boundaryLayerWorkflow.state.station1.param.amplz, hkt.scalar, tt, rtt.scalar,
                        amplt, amcrit);

      //---- punch out early if there is no amplification here
      if (ax_result.ax <= 0.0) {
        return true;
      }

      //---- set sensitivity of ax(a2)
      ax_result.ax_a2 =
          (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
           ax_result.ax_rt2 * rtt.t()) *
              tt_a2 +
          (ax_result.ax_hk2 * hkt.d()) * dt_a2 +
          (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_a2 +
          ax_result.ax_a2 * amplt_a2;

      //---- residual for implicit ampl2 definition (amplification equation)
      res = boundaryLayerWorkflow.state.station2.param.amplz - boundaryLayerWorkflow.state.station1.param.amplz -
            ax_result.ax * (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz);
      res_a2 = 1.0 - ax_result.ax_a2 * (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz);

      da2 = -res / res_a2;

      rlx = 1.0;
      dxt = xt_a2 * da2;

      if (rlx * fabs(dxt / (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz)) > 0.05) {
        rlx = 0.05 * fabs((boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz) / dxt);
      }

      if (rlx * fabs(da2) > 1.0) {
        rlx = 1.0 * fabs(1.0 / da2);
      }

      //---- check if converged
      if (fabs(da2) < daeps) {
        return true;
      }

      if ((boundaryLayerWorkflow.state.station2.param.amplz > amcrit &&
           boundaryLayerWorkflow.state.station2.param.amplz + rlx * da2 < amcrit) ||
          (boundaryLayerWorkflow.state.station2.param.amplz < amcrit &&
           boundaryLayerWorkflow.state.station2.param.amplz + rlx * da2 > amcrit)) {
        //------ limited newton step so ampl2 doesn't step across amcrit either
        // way
        boundaryLayerWorkflow.state.station2.param.amplz = amcrit;
      } else {
        //------ regular newton step
        boundaryLayerWorkflow.state.station2.param.amplz = boundaryLayerWorkflow.state.station2.param.amplz + rlx * da2;
      }
    }
    return false;
  };

  if (!iterateAmplification()) {
    // TRACE("trchek2 - n2 convergence failed\n");
    writeString("trchek2 - n2 convergence failed\n");
    if (isCancelled())
      return false;
  }

  //---- test for free or forced transition
  trfree = (boundaryLayerWorkflow.state.station2.param.amplz >= amcrit);
  trforc = (xiforc > boundaryLayerWorkflow.state.station1.param.xz) && (xiforc <= boundaryLayerWorkflow.state.station2.param.xz);

  //---- set transition interval flag
  tran = (trforc || trfree);

  if (!tran)
    return false;

  //---- resolve if both forced and free transition
  if (trfree && trforc) {
    trforc = xiforc < xt;
    trfree = xiforc >= xt;
  }

  if (trforc) {
    //----- if forced transition, then xt is prescribed,
    //-     no sense calculating the sensitivities, since we know them...
    xt = xiforc;
    xt_a1 = 0.0;
    xt_x1 = 0.0;
    xt_t1 = 0.0;
    xt_d1 = 0.0;
    xt_u1 = 0.0;
    xt_x2 = 0.0;
    xt_t2 = 0.0;
    xt_d2 = 0.0;
    xt_u2 = 0.0;
    xt_ms = 0.0;
    xt_re = 0.0;
    xt_xf = 1.0;
    return true;
  }

  //---- free transition ... set sensitivities of xt

  xt_x1 = (1 - wf);
  tt_t1 = (1 - wf);
  dt_d1 = (1 - wf);
  ut_u1 = (1 - wf);

  xt_x2 = wf;
  tt_t2 = wf;
  dt_d2 = wf;
  ut_u2 = wf;

  xt_a1 = (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz) * wf_a1;
  tt_a1 = (boundaryLayerWorkflow.state.station2.param.tz - boundaryLayerWorkflow.state.station1.param.tz) * wf_a1;
  dt_a1 = (boundaryLayerWorkflow.state.station2.param.dz - boundaryLayerWorkflow.state.station1.param.dz) * wf_a1;
  ut_a1 = (boundaryLayerWorkflow.state.station2.param.uz - boundaryLayerWorkflow.state.station1.param.uz) * wf_a1;

  xt_x1 += (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz) * wf_x1;
  tt_x1 = (boundaryLayerWorkflow.state.station2.param.tz - boundaryLayerWorkflow.state.station1.param.tz) * wf_x1;
  dt_x1 = (boundaryLayerWorkflow.state.station2.param.dz - boundaryLayerWorkflow.state.station1.param.dz) * wf_x1;
  ut_x1 = (boundaryLayerWorkflow.state.station2.param.uz - boundaryLayerWorkflow.state.station1.param.uz) * wf_x1;

  xt_x2 += (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz) * wf_x2;
  tt_x2 = (boundaryLayerWorkflow.state.station2.param.tz - boundaryLayerWorkflow.state.station1.param.tz) * wf_x2;
  dt_x2 = (boundaryLayerWorkflow.state.station2.param.dz - boundaryLayerWorkflow.state.station1.param.dz) * wf_x2;
  ut_x2 = (boundaryLayerWorkflow.state.station2.param.uz - boundaryLayerWorkflow.state.station1.param.uz) * wf_x2;

  xt_xf = (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz) * wf_xf;

  //---- at this point, ax = ax( hk1, t1, rt1, a1, hkt, tt, rtt, at )
  blData::blVector hkt = boundaryLayerWorkflow.state.station2.hkz;
  blData::blVector rtt = boundaryLayerWorkflow.state.station2.rtz;

  //---- set sensitivities of ax( t1 d1 u1 a1 t2 d2 u2 a2 ms re )
  double ax_t1 = ax_result.ax_hk1 * boundaryLayerWorkflow.state.station1.hkz.t() + ax_result.ax_t1 +
                 ax_result.ax_rt1 * boundaryLayerWorkflow.state.station1.rtz.t() +
                 (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
                  ax_result.ax_rt2 * rtt.t()) *
                     tt_t1;
  double ax_d1 =
      ax_result.ax_hk1 * boundaryLayerWorkflow.state.station1.hkz.d() + (ax_result.ax_hk2 * hkt.d()) * dt_d1;
  double ax_u1 =
      ax_result.ax_hk1 * boundaryLayerWorkflow.state.station1.hkz.u() + ax_result.ax_rt1 * boundaryLayerWorkflow.state.station1.rtz.u() +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_u1;
  double ax_a1 =
      ax_result.ax_a1 +
      (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
       ax_result.ax_rt2 * rtt.t()) *
          tt_a1 +
      (ax_result.ax_hk2 * hkt.d()) * dt_a1 +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_a1;
  double ax_x1 =
      (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
       ax_result.ax_rt2 * rtt.t()) *
          tt_x1 +
      (ax_result.ax_hk2 * hkt.d()) * dt_x1 +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_x1;

  double ax_t2 = (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
                  ax_result.ax_rt2 * rtt.t()) *
                 tt_t2;
  double ax_d2 = (ax_result.ax_hk2 * hkt.d()) * dt_d2;
  double ax_u2 =
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_u2;
  double ax_a2 =
      ax_result.ax_a2 * amplt_a2 +
      (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
       ax_result.ax_rt2 * rtt.t()) *
          tt_a2 +
      (ax_result.ax_hk2 * hkt.d()) * dt_a2 +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_a2;
  double ax_x2 =
      (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
       ax_result.ax_rt2 * rtt.t()) *
          tt_x2 +
      (ax_result.ax_hk2 * hkt.d()) * dt_x2 +
      (ax_result.ax_hk2 * hkt.u() + ax_result.ax_rt2 * rtt.u()) * ut_x2;

  double ax_ms = ax_result.ax_hk2 * hkt.ms() + ax_result.ax_rt2 * rtt.ms() +
                 ax_result.ax_hk1 * boundaryLayerWorkflow.state.station1.hkz.ms() +
                 ax_result.ax_rt1 * boundaryLayerWorkflow.state.station1.rtz.ms();
  double ax_re =
      ax_result.ax_rt2 * rtt.re() + ax_result.ax_rt1 * boundaryLayerWorkflow.state.station1.rtz.re();

  //---- set sensitivities of residual res
  z_ax = -(boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz);

  z_a1 = z_ax * ax_a1 - 1.0;
  z_t1 = z_ax * ax_t1;
  z_d1 = z_ax * ax_d1;
  z_u1 = z_ax * ax_u1;
  z_x1 = z_ax * ax_x1 + ax_result.ax;

  z_a2 = z_ax * ax_a2 + 1.0;
  z_t2 = z_ax * ax_t2;
  z_d2 = z_ax * ax_d2;
  z_u2 = z_ax * ax_u2;
  z_x2 = z_ax * ax_x2 - ax_result.ax;

  z_ms = z_ax * ax_ms;
  z_re = z_ax * ax_re;

  //---- set sensitivities of xt, with res being stationary for a2 constraint
  xt_a1 = xt_a1 - (xt_a2 / z_a2) * z_a1;
  xt_t1 = -(xt_a2 / z_a2) * z_t1;
  xt_d1 = -(xt_a2 / z_a2) * z_d1;
  xt_u1 = -(xt_a2 / z_a2) * z_u1;
  xt_x1 = xt_x1 - (xt_a2 / z_a2) * z_x1;
  xt_t2 = -(xt_a2 / z_a2) * z_t2;
  xt_d2 = -(xt_a2 / z_a2) * z_d2;
  xt_u2 = -(xt_a2 / z_a2) * z_u2;
  xt_x2 = xt_x2 - (xt_a2 / z_a2) * z_x2;
  xt_ms = -(xt_a2 / z_a2) * z_ms;
  xt_re = -(xt_a2 / z_a2) * z_re;
  xt_xf = 0.0;

  return true;
}


bool XFoil::trdif() {
  //-----------------------------------------------
  //     sets up the newton system governing the
  //     transition interval.  equations governing
  //     the  laminar  part  x1 < xi < xt  and
  //     the turbulent part  xt < xi < x2
  //     are simply summed.
  //-----------------------------------------------
  Matrix<double, 4, 5> bl1, bl2, bt1, bt2;
  Vector<double, 4> blrez, blm, blr, blx, btrez, btm, btr, btx;

  double tt, tt_a1, tt_x1, tt_x2, tt_t1, tt_t2, tt_d1, tt_d2, tt_u1, tt_u2;
  double tt_ms, tt_re, tt_xf, dt, dt_a1, dt_x1, dt_x2, dt_t1, dt_t2;
  double dt_d1, dt_d2, dt_u1, dt_u2, dt_ms, dt_re, dt_xf;
  double ut, ut_a1, ut_x1, ut_x2, ut_t1, ut_t2, ut_d1, ut_d2, ut_u1, ut_u2;
  double ut_ms, ut_re, ut_xf;
  double st, st_tt, st_dt, st_ut, st_ms, st_re, st_a1, st_x1, st_x2, st_t1,
      st_t2;
  double st_d1, st_d2, st_u1, st_u2, st_xf;
  double ctr, ctr_hk2;

  saveblData(1);
  saveblData(2);

  //---- weighting factors for linear interpolation to transition point
  double wf2 = (xt - boundaryLayerWorkflow.state.station1.param.xz) / (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz);
  double wf2_xt = 1.0 / (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz);

  double wf2_a1 = wf2_xt * xt_a1;
  double wf2_x1 =
      wf2_xt * xt_x1 + (wf2 - 1.0) / (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz);
  double wf2_x2 = wf2_xt * xt_x2 - wf2 / (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz);
  double wf2_t1 = wf2_xt * xt_t1;
  double wf2_t2 = wf2_xt * xt_t2;
  double wf2_d1 = wf2_xt * xt_d1;
  double wf2_d2 = wf2_xt * xt_d2;
  double wf2_u1 = wf2_xt * xt_u1;
  double wf2_u2 = wf2_xt * xt_u2;
  double wf2_ms = wf2_xt * xt_ms;
  double wf2_re = wf2_xt * xt_re;
  double wf2_xf = wf2_xt * xt_xf;

  double wf1 = 1.0 - wf2;
  double wf1_a1 = -wf2_a1;
  double wf1_x1 = -wf2_x1;
  double wf1_x2 = -wf2_x2;
  double wf1_t1 = -wf2_t1;
  double wf1_t2 = -wf2_t2;
  double wf1_d1 = -wf2_d1;
  double wf1_d2 = -wf2_d2;
  double wf1_u1 = -wf2_u1;
  double wf1_u2 = -wf2_u2;
  double wf1_ms = -wf2_ms;
  double wf1_re = -wf2_re;
  double wf1_xf = -wf2_xf;

  //-----interpolate primary variables to transition point
  tt = boundaryLayerWorkflow.state.station1.param.tz * wf1 + boundaryLayerWorkflow.state.station2.param.tz * wf2;
  tt_a1 = boundaryLayerWorkflow.state.station1.param.tz * wf1_a1 + boundaryLayerWorkflow.state.station2.param.tz * wf2_a1;
  tt_x1 = boundaryLayerWorkflow.state.station1.param.tz * wf1_x1 + boundaryLayerWorkflow.state.station2.param.tz * wf2_x1;
  tt_x2 = boundaryLayerWorkflow.state.station1.param.tz * wf1_x2 + boundaryLayerWorkflow.state.station2.param.tz * wf2_x2;
  tt_t1 = boundaryLayerWorkflow.state.station1.param.tz * wf1_t1 + boundaryLayerWorkflow.state.station2.param.tz * wf2_t1 + wf1;
  tt_t2 = boundaryLayerWorkflow.state.station1.param.tz * wf1_t2 + boundaryLayerWorkflow.state.station2.param.tz * wf2_t2 + wf2;
  tt_d1 = boundaryLayerWorkflow.state.station1.param.tz * wf1_d1 + boundaryLayerWorkflow.state.station2.param.tz * wf2_d1;
  tt_d2 = boundaryLayerWorkflow.state.station1.param.tz * wf1_d2 + boundaryLayerWorkflow.state.station2.param.tz * wf2_d2;
  tt_u1 = boundaryLayerWorkflow.state.station1.param.tz * wf1_u1 + boundaryLayerWorkflow.state.station2.param.tz * wf2_u1;
  tt_u2 = boundaryLayerWorkflow.state.station1.param.tz * wf1_u2 + boundaryLayerWorkflow.state.station2.param.tz * wf2_u2;
  tt_ms = boundaryLayerWorkflow.state.station1.param.tz * wf1_ms + boundaryLayerWorkflow.state.station2.param.tz * wf2_ms;
  tt_re = boundaryLayerWorkflow.state.station1.param.tz * wf1_re + boundaryLayerWorkflow.state.station2.param.tz * wf2_re;
  tt_xf = boundaryLayerWorkflow.state.station1.param.tz * wf1_xf + boundaryLayerWorkflow.state.station2.param.tz * wf2_xf;

  dt = boundaryLayerWorkflow.state.station1.param.dz * wf1 + boundaryLayerWorkflow.state.station2.param.dz * wf2;
  dt_a1 = boundaryLayerWorkflow.state.station1.param.dz * wf1_a1 + boundaryLayerWorkflow.state.station2.param.dz * wf2_a1;
  dt_x1 = boundaryLayerWorkflow.state.station1.param.dz * wf1_x1 + boundaryLayerWorkflow.state.station2.param.dz * wf2_x1;
  dt_x2 = boundaryLayerWorkflow.state.station1.param.dz * wf1_x2 + boundaryLayerWorkflow.state.station2.param.dz * wf2_x2;
  dt_t1 = boundaryLayerWorkflow.state.station1.param.dz * wf1_t1 + boundaryLayerWorkflow.state.station2.param.dz * wf2_t1;
  dt_t2 = boundaryLayerWorkflow.state.station1.param.dz * wf1_t2 + boundaryLayerWorkflow.state.station2.param.dz * wf2_t2;
  dt_d1 = boundaryLayerWorkflow.state.station1.param.dz * wf1_d1 + boundaryLayerWorkflow.state.station2.param.dz * wf2_d1 + wf1;
  dt_d2 = boundaryLayerWorkflow.state.station1.param.dz * wf1_d2 + boundaryLayerWorkflow.state.station2.param.dz * wf2_d2 + wf2;
  dt_u1 = boundaryLayerWorkflow.state.station1.param.dz * wf1_u1 + boundaryLayerWorkflow.state.station2.param.dz * wf2_u1;
  dt_u2 = boundaryLayerWorkflow.state.station1.param.dz * wf1_u2 + boundaryLayerWorkflow.state.station2.param.dz * wf2_u2;
  dt_ms = boundaryLayerWorkflow.state.station1.param.dz * wf1_ms + boundaryLayerWorkflow.state.station2.param.dz * wf2_ms;
  dt_re = boundaryLayerWorkflow.state.station1.param.dz * wf1_re + boundaryLayerWorkflow.state.station2.param.dz * wf2_re;
  dt_xf = boundaryLayerWorkflow.state.station1.param.dz * wf1_xf + boundaryLayerWorkflow.state.station2.param.dz * wf2_xf;

  ut = boundaryLayerWorkflow.state.station1.param.uz * wf1 + boundaryLayerWorkflow.state.station2.param.uz * wf2;
  ut_a1 = boundaryLayerWorkflow.state.station1.param.uz * wf1_a1 + boundaryLayerWorkflow.state.station2.param.uz * wf2_a1;
  ut_x1 = boundaryLayerWorkflow.state.station1.param.uz * wf1_x1 + boundaryLayerWorkflow.state.station2.param.uz * wf2_x1;
  ut_x2 = boundaryLayerWorkflow.state.station1.param.uz * wf1_x2 + boundaryLayerWorkflow.state.station2.param.uz * wf2_x2;
  ut_t1 = boundaryLayerWorkflow.state.station1.param.uz * wf1_t1 + boundaryLayerWorkflow.state.station2.param.uz * wf2_t1;
  ut_t2 = boundaryLayerWorkflow.state.station1.param.uz * wf1_t2 + boundaryLayerWorkflow.state.station2.param.uz * wf2_t2;
  ut_d1 = boundaryLayerWorkflow.state.station1.param.uz * wf1_d1 + boundaryLayerWorkflow.state.station2.param.uz * wf2_d1;
  ut_d2 = boundaryLayerWorkflow.state.station1.param.uz * wf1_d2 + boundaryLayerWorkflow.state.station2.param.uz * wf2_d2;
  ut_u1 = boundaryLayerWorkflow.state.station1.param.uz * wf1_u1 + boundaryLayerWorkflow.state.station2.param.uz * wf2_u1 + wf1;
  ut_u2 = boundaryLayerWorkflow.state.station1.param.uz * wf1_u2 + boundaryLayerWorkflow.state.station2.param.uz * wf2_u2 + wf2;
  ut_ms = boundaryLayerWorkflow.state.station1.param.uz * wf1_ms + boundaryLayerWorkflow.state.station2.param.uz * wf2_ms;
  ut_re = boundaryLayerWorkflow.state.station1.param.uz * wf1_re + boundaryLayerWorkflow.state.station2.param.uz * wf2_re;
  ut_xf = boundaryLayerWorkflow.state.station1.param.uz * wf1_xf + boundaryLayerWorkflow.state.station2.param.uz * wf2_xf;

  //---- set primary "t" variables at xt  (really placed into "2" variables)
  boundaryLayerWorkflow.state.station2.param.xz = xt;
  boundaryLayerWorkflow.state.station2.param.tz = tt;
  boundaryLayerWorkflow.state.station2.param.dz = dt;
  boundaryLayerWorkflow.state.station2.param.uz = ut;

  boundaryLayerWorkflow.state.station2.param.amplz = amcrit;
  boundaryLayerWorkflow.state.station2.param.sz = 0.0;

  //---- calculate laminar secondary "t" variables
  blkin(boundaryLayerWorkflow.state);
  boundaryLayerWorkflow.state.station2 =
      boundaryLayerWorkflow.blvar(boundaryLayerWorkflow.state.station2,
                                  FlowRegimeEnum::Laminar);

  //---- calculate x1-xt midpoint cfm value
  SkinFrictionCoefficients laminarSkinFriction =
      boundaryLayerWorkflow.blmid(*this, FlowRegimeEnum::Laminar);

  //=    at this point, all "2" variables are really "t" variables at xt

  //---- set up newton system for dam, dth, dds, due, dxi  at  x1 and xt
  boundaryLayerWorkflow.blc = blDiffSolver.solve(FlowRegimeEnum::Laminar, boundaryLayerWorkflow.state, laminarSkinFriction, amcrit);

  //---- the current newton system is in terms of "1" and "t" variables,
  //-    so calculate its equivalent in terms of "1" and "2" variables.
  //-    in other words, convert residual sensitivities wrt "t" variables
  //-    into sensitivities wrt "1" and "2" variables.  the amplification
  //-    equation is unnecessary here, so the k=1 row is left empty.
  for (int k = 1; k < 3; k++) {
    blrez[k] = boundaryLayerWorkflow.blc.rhs[k];
    blm[k] = boundaryLayerWorkflow.blc.d_msq[k] + boundaryLayerWorkflow.blc.a2(k, 1) * tt_ms + boundaryLayerWorkflow.blc.a2(k, 2) * dt_ms +
             boundaryLayerWorkflow.blc.a2(k, 3) * ut_ms + boundaryLayerWorkflow.blc.a2(k, 4) * xt_ms;
    blr[k] = boundaryLayerWorkflow.blc.d_re[k] + boundaryLayerWorkflow.blc.a2(k, 1) * tt_re + boundaryLayerWorkflow.blc.a2(k, 2) * dt_re +
             boundaryLayerWorkflow.blc.a2(k, 3) * ut_re + boundaryLayerWorkflow.blc.a2(k, 4) * xt_re;
    blx[k] = boundaryLayerWorkflow.blc.d_xi[k] + boundaryLayerWorkflow.blc.a2(k, 1) * tt_xf + boundaryLayerWorkflow.blc.a2(k, 2) * dt_xf +
             boundaryLayerWorkflow.blc.a2(k, 3) * ut_xf + boundaryLayerWorkflow.blc.a2(k, 4) * xt_xf;

    bl1(k, 0) = boundaryLayerWorkflow.blc.a1(k, 0) + boundaryLayerWorkflow.blc.a2(k, 1) * tt_a1 + boundaryLayerWorkflow.blc.a2(k, 2) * dt_a1 +
                boundaryLayerWorkflow.blc.a2(k, 3) * ut_a1 + boundaryLayerWorkflow.blc.a2(k, 4) * xt_a1;
    bl1(k, 1) = boundaryLayerWorkflow.blc.a1(k, 1) + boundaryLayerWorkflow.blc.a2(k, 1) * tt_t1 + boundaryLayerWorkflow.blc.a2(k, 2) * dt_t1 +
                boundaryLayerWorkflow.blc.a2(k, 3) * ut_t1 + boundaryLayerWorkflow.blc.a2(k, 4) * xt_t1;
    bl1(k, 2) = boundaryLayerWorkflow.blc.a1(k, 2) + boundaryLayerWorkflow.blc.a2(k, 1) * tt_d1 + boundaryLayerWorkflow.blc.a2(k, 2) * dt_d1 +
                boundaryLayerWorkflow.blc.a2(k, 3) * ut_d1 + boundaryLayerWorkflow.blc.a2(k, 4) * xt_d1;
    bl1(k, 3) = boundaryLayerWorkflow.blc.a1(k, 3) + boundaryLayerWorkflow.blc.a2(k, 1) * tt_u1 + boundaryLayerWorkflow.blc.a2(k, 2) * dt_u1 +
                boundaryLayerWorkflow.blc.a2(k, 3) * ut_u1 + boundaryLayerWorkflow.blc.a2(k, 4) * xt_u1;
    bl1(k, 4) = boundaryLayerWorkflow.blc.a1(k, 4) + boundaryLayerWorkflow.blc.a2(k, 1) * tt_x1 + boundaryLayerWorkflow.blc.a2(k, 2) * dt_x1 +
                boundaryLayerWorkflow.blc.a2(k, 3) * ut_x1 + boundaryLayerWorkflow.blc.a2(k, 4) * xt_x1;

    bl2(k, 0) = 0.0;
    bl2(k, 1) = boundaryLayerWorkflow.blc.a2(k, 1) * tt_t2 + boundaryLayerWorkflow.blc.a2(k, 2) * dt_t2 + boundaryLayerWorkflow.blc.a2(k, 3) * ut_t2 +
                boundaryLayerWorkflow.blc.a2(k, 4) * xt_t2;
    bl2(k, 2) = boundaryLayerWorkflow.blc.a2(k, 1) * tt_d2 + boundaryLayerWorkflow.blc.a2(k, 2) * dt_d2 + boundaryLayerWorkflow.blc.a2(k, 3) * ut_d2 +
                boundaryLayerWorkflow.blc.a2(k, 4) * xt_d2;
    bl2(k, 3) = boundaryLayerWorkflow.blc.a2(k, 1) * tt_u2 + boundaryLayerWorkflow.blc.a2(k, 2) * dt_u2 + boundaryLayerWorkflow.blc.a2(k, 3) * ut_u2 +
                boundaryLayerWorkflow.blc.a2(k, 4) * xt_u2;
    bl2(k, 4) = boundaryLayerWorkflow.blc.a2(k, 1) * tt_x2 + boundaryLayerWorkflow.blc.a2(k, 2) * dt_x2 + boundaryLayerWorkflow.blc.a2(k, 3) * ut_x2 +
                boundaryLayerWorkflow.blc.a2(k, 4) * xt_x2;
  }

  //**** second, set up turbulent part between xt and x2  ****

  //---- calculate equilibrium shear coefficient cqt at transition point
  boundaryLayerWorkflow.state.station2 =
      boundaryLayerWorkflow.blvar(boundaryLayerWorkflow.state.station2,
                                  FlowRegimeEnum::Turbulent);

  //---- set initial shear coefficient value st at transition point
  //-    ( note that cq2, cq2_t2, etc. are really "cqt", "cqt_tt", etc.)

  ctr = 1.8 * exp(-3.3 / (boundaryLayerWorkflow.state.station2.hkz.scalar - 1.0));
  ctr_hk2 = ctr * 3.3 / (boundaryLayerWorkflow.state.station2.hkz.scalar - 1.0) / (boundaryLayerWorkflow.state.station2.hkz.scalar - 1.0);

  st = ctr * boundaryLayerWorkflow.state.station2.cqz.scalar;
  st_tt =
      ctr * boundaryLayerWorkflow.state.station2.cqz.t() + boundaryLayerWorkflow.state.station2.cqz.scalar * ctr_hk2 * boundaryLayerWorkflow.state.station2.hkz.t();
  st_dt =
      ctr * boundaryLayerWorkflow.state.station2.cqz.d() + boundaryLayerWorkflow.state.station2.cqz.scalar * ctr_hk2 * boundaryLayerWorkflow.state.station2.hkz.d();
  st_ut =
      ctr * boundaryLayerWorkflow.state.station2.cqz.u() + boundaryLayerWorkflow.state.station2.cqz.scalar * ctr_hk2 * boundaryLayerWorkflow.state.station2.hkz.u();
  st_ms =
      ctr * boundaryLayerWorkflow.state.station2.cqz.ms() + boundaryLayerWorkflow.state.station2.cqz.scalar * ctr_hk2 * boundaryLayerWorkflow.state.station2.hkz.ms();
  st_re = ctr * boundaryLayerWorkflow.state.station2.cqz.re();

  //---- calculate st sensitivities wrt the actual "1" and "2" variables
  st_a1 = st_tt * tt_a1 + st_dt * dt_a1 + st_ut * ut_a1;
  st_x1 = st_tt * tt_x1 + st_dt * dt_x1 + st_ut * ut_x1;
  st_x2 = st_tt * tt_x2 + st_dt * dt_x2 + st_ut * ut_x2;
  st_t1 = st_tt * tt_t1 + st_dt * dt_t1 + st_ut * ut_t1;
  st_t2 = st_tt * tt_t2 + st_dt * dt_t2 + st_ut * ut_t2;
  st_d1 = st_tt * tt_d1 + st_dt * dt_d1 + st_ut * ut_d1;
  st_d2 = st_tt * tt_d2 + st_dt * dt_d2 + st_ut * ut_d2;
  st_u1 = st_tt * tt_u1 + st_dt * dt_u1 + st_ut * ut_u1;
  st_u2 = st_tt * tt_u2 + st_dt * dt_u2 + st_ut * ut_u2;
  st_ms = st_tt * tt_ms + st_dt * dt_ms + st_ut * ut_ms + st_ms;
  st_re = st_tt * tt_re + st_dt * dt_re + st_ut * ut_re + st_re;
  st_xf = st_tt * tt_xf + st_dt * dt_xf + st_ut * ut_xf;

  boundaryLayerWorkflow.state.station2.param.amplz = 0.0;
  boundaryLayerWorkflow.state.station2.param.sz = st;

  //---- recalculate turbulent secondary "t" variables using proper cti
  boundaryLayerWorkflow.state.station2 =
      boundaryLayerWorkflow.blvar(boundaryLayerWorkflow.state.station2,
                                  FlowRegimeEnum::Turbulent);

  boundaryLayerWorkflow.state.stepbl();
  restoreblData(2);

  //---- calculate xt-x2 midpoint cfm value
  SkinFrictionCoefficients turbulentSkinFriction =
      boundaryLayerWorkflow.blmid(*this, FlowRegimeEnum::Turbulent);

  //---- set up newton system for dct, dth, dds, due, dxi  at  xt and x2
  boundaryLayerWorkflow.blc = blDiffSolver.solve(FlowRegimeEnum::Turbulent, boundaryLayerWorkflow.state, turbulentSkinFriction, amcrit);

  //---- convert sensitivities wrt "t" variables into sensitivities
  //-    wrt "1" and "2" variables as done before for the laminar part
  Matrix<double, 5, 5> bt1_right =
      Matrix<double, 5, 5>{{st_a1, st_t1, st_d1, st_u1, st_x1},
                           {tt_a1, tt_t1, tt_d1, tt_u1, tt_x1},
                           {dt_a1, dt_t1, dt_d1, dt_u1, dt_x1},
                           {ut_a1, ut_t1, ut_d1, ut_u1, ut_x1},
                           {xt_a1, xt_t1, xt_d1, xt_u1, xt_x1}};

  Matrix<double, 5, 5> bt2_right =
      Matrix<double, 5, 5>{{0, st_t2, st_d2, st_u2, st_x2},
                           {0, tt_t2, tt_d2, tt_u2, tt_x2},
                           {0, dt_t2, dt_d2, dt_u2, dt_x2},
                           {0, ut_t2, ut_d2, ut_u2, ut_x2},
                           {0, xt_t2, xt_d2, xt_u2, xt_x2}};
  bt1.block(0, 0, 3, 5) = boundaryLayerWorkflow.blc.a1.block(0, 0, 3, 5) * bt1_right;
  bt2.block(0, 0, 3, 5) = boundaryLayerWorkflow.blc.a1.block(0, 0, 3, 5) * bt2_right;
  bt2 += boundaryLayerWorkflow.blc.a2;
  for (int k = 0; k < 3; k++) {
    btrez[k] = boundaryLayerWorkflow.blc.rhs[k];
    btm[k] = boundaryLayerWorkflow.blc.d_msq[k] + boundaryLayerWorkflow.blc.a1(k, 0) * st_ms + boundaryLayerWorkflow.blc.a1(k, 1) * tt_ms +
             boundaryLayerWorkflow.blc.a1(k, 2) * dt_ms + boundaryLayerWorkflow.blc.a1(k, 3) * ut_ms + boundaryLayerWorkflow.blc.a1(k, 4) * xt_ms;
    btr[k] = boundaryLayerWorkflow.blc.d_re[k] + boundaryLayerWorkflow.blc.a1(k, 0) * st_re + boundaryLayerWorkflow.blc.a1(k, 1) * tt_re +
             boundaryLayerWorkflow.blc.a1(k, 2) * dt_re + boundaryLayerWorkflow.blc.a1(k, 3) * ut_re + boundaryLayerWorkflow.blc.a1(k, 4) * xt_re;
    btx[k] = boundaryLayerWorkflow.blc.d_xi[k] + boundaryLayerWorkflow.blc.a1(k, 0) * st_xf + boundaryLayerWorkflow.blc.a1(k, 1) * tt_xf +
             boundaryLayerWorkflow.blc.a1(k, 2) * dt_xf + boundaryLayerWorkflow.blc.a1(k, 3) * ut_xf + boundaryLayerWorkflow.blc.a1(k, 4) * xt_xf;
  }

  //---- add up laminar and turbulent parts to get final system
  //-    in terms of honest-to-god "1" and "2" variables.
  boundaryLayerWorkflow.blc.rhs[0] = btrez[0];
  boundaryLayerWorkflow.blc.rhs[1] = blrez[1] + btrez[1];
  boundaryLayerWorkflow.blc.rhs[2] = blrez[2] + btrez[2];
  boundaryLayerWorkflow.blc.d_msq[0] = btm[0];
  boundaryLayerWorkflow.blc.d_msq[1] = blm[1] + btm[1];
  boundaryLayerWorkflow.blc.d_msq[2] = blm[2] + btm[2];
  boundaryLayerWorkflow.blc.d_re[0] = btr[0];
  boundaryLayerWorkflow.blc.d_re[1] = blr[1] + btr[1];
  boundaryLayerWorkflow.blc.d_re[2] = blr[2] + btr[2];
  boundaryLayerWorkflow.blc.d_xi[0] = btx[0];
  boundaryLayerWorkflow.blc.d_xi[1] = blx[1] + btx[1];
  boundaryLayerWorkflow.blc.d_xi[2] = blx[2] + btx[2];
  boundaryLayerWorkflow.blc.a1.row(0) = bt1.row(0);
  boundaryLayerWorkflow.blc.a2.row(0) = bt2.row(0);
  boundaryLayerWorkflow.blc.a1.middleRows(1, 2) = bl1.middleRows(1, 2) + bt1.middleRows(1, 2);
  boundaryLayerWorkflow.blc.a2.middleRows(1, 2) = bl2.middleRows(1, 2) + bt2.middleRows(1, 2);

  //---- to be sanitary, restore "1" quantities which got clobbered
  //-    in all of the numerical gymnastics above.  the "2" variables
  //-    were already restored for the xt-x2 differencing part.
  //	for (icom=1; icom<=ncom;icom++){
  //		com1[icom] = c1sav[icom];
  //	}
  restoreblData(1);

  return true;
}


bool XFoil::ueset() {
  //---------------------------------------------------------
  //     sets ue from inviscid ue plus all source influence
  //---------------------------------------------------------
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 0; ibl < boundaryLayerWorkflow.lattice.get(is).stationCount - 1; ++ibl) {
      double dui = 0.0;
      for (int js = 1; js <= 2; js++) {
        for (int jbl = 0; jbl < boundaryLayerWorkflow.lattice.get(js).stationCount - 1; ++jbl) {
          double ue_m = -boundaryLayerWorkflow.lattice.get(is).panelInfluenceFactor[ibl] * boundaryLayerWorkflow.lattice.get(js).panelInfluenceFactor[jbl] *
                        dij(boundaryLayerWorkflow.lattice.get(is).stationToPanel[ibl],
                            boundaryLayerWorkflow.lattice.get(js).stationToPanel[jbl]);
          dui += ue_m * boundaryLayerWorkflow.lattice.get(js).profiles.massFlux[jbl];
        }
      }
      boundaryLayerWorkflow.lattice.get(is).profiles.edgeVelocity[ibl] =
          boundaryLayerWorkflow.lattice.get(is).inviscidEdgeVelocity[ibl] + dui;
    }
  }
  return true;
}


bool XFoil::xicalc() {
  //-------------------------------------------------------------
  //     sets bl arc length array on each airfoil side and wake
  //-------------------------------------------------------------

  
    for (int ibl = 0; ibl <= boundaryLayerWorkflow.lattice.top.trailingEdgeIndex; ++ibl) {
      boundaryLayerWorkflow.lattice.top.arcLengthCoordinates[ibl] = sst - foil.foil_shape.spline_length[boundaryLayerWorkflow.lattice.get(1).stationToPanel[ibl]];
    }
  
    for (int ibl = 0; ibl <= boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex; ++ibl) {
      boundaryLayerWorkflow.lattice.bottom.arcLengthCoordinates[ibl] = foil.foil_shape.spline_length[boundaryLayerWorkflow.lattice.get(2).stationToPanel[ibl]] - sst;
    }

    // Wake: start from TE, duplicate TE value at first wake station
    boundaryLayerWorkflow.lattice.bottom.arcLengthCoordinates[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex + 1] = boundaryLayerWorkflow.lattice.bottom.arcLengthCoordinates[boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex];
    for (int ibl = boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex + 2; ibl < boundaryLayerWorkflow.lattice.bottom.stationCount; ++ibl) {
      boundaryLayerWorkflow.lattice.bottom.arcLengthCoordinates[ibl] = boundaryLayerWorkflow.lattice.bottom.arcLengthCoordinates[ibl - 1] +
                          (foil.wake_shape.points.col(boundaryLayerWorkflow.lattice.get(2).stationToPanel[ibl]) -
                           foil.wake_shape.points.col(boundaryLayerWorkflow.lattice.get(2).stationToPanel[ibl - 1]))
                              .norm();
    }
  

  //---- trailing edge flap length to te gap ratio
  const double telrat = 2.50;

  //---- set up parameters for te flap cubics

  const int point_count = foil.foil_shape.n;
  const double crosp = MathUtil::cross2(foil.foil_shape.dpoints_ds.col(point_count - 1).normalized(),
                                        foil.foil_shape.dpoints_ds.col(0).normalized());
  double dwdxte = crosp / sqrt(1.0 - crosp * crosp);

  //---- limit cubic to avoid absurd te gap widths
  dwdxte = std::max(dwdxte, -3.0 / telrat);
  dwdxte = std::min(dwdxte, 3.0 / telrat);

  const double aa = 3.0 + telrat * dwdxte;
  const double bb = -2.0 - telrat * dwdxte;

  if (foil.edge.sharp) {
    for (int iw0 = 0; iw0 < foil.wake_shape.n; iw0++)
      wgap[iw0] = 0.0;
  }

  else {
    //----- set te flap (wake gap) array (0-based: iw0=0..foil.wake_shape.n-1)
    for (int iw0 = 0; iw0 < foil.wake_shape.n; iw0++) {
      const int te_bot_0b = boundaryLayerWorkflow.lattice.bottom.trailingEdgeIndex; // 0-based TE for array indexing
      const double zn = 1.0 - (boundaryLayerWorkflow.lattice.bottom.arcLengthCoordinates[te_bot_0b + (iw0 + 1)] -
                               boundaryLayerWorkflow.lattice.bottom.arcLengthCoordinates[te_bot_0b]) /
                                (telrat * foil.edge.ante);
      wgap[iw0] = 0.0;
      if (zn >= 0.0)
        wgap[iw0] = foil.edge.ante * (aa + bb * zn) * zn * zn;
    }
  }
  return true;
}


/** -----------------------------------------------------
 * 	   sets forced-transition bl coordinate locations.
 * ----------------------------------------------------- */
double XFoil::xifset(int is) {
  std::stringstream ss;
  VectorXd w1 = VectorXd::Zero(6 * IQX);
  VectorXd w2 = VectorXd::Zero(6 * IQX);
  VectorXd w3 = VectorXd::Zero(6 * IQX);
  VectorXd w4 = VectorXd::Zero(6 * IQX);
  double str;
  const int point_count = foil.foil_shape.n;

  if (boundaryLayerWorkflow.lattice.get(is).transitionLocation >= 1.0) {
    return boundaryLayerWorkflow.lattice.get(is).arcLengthCoordinates[boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex];
  }

  Vector2d point_chord = foil.edge.point_te - foil.edge.point_le;

  //---- calculate chord-based x/c, y/c
  for (int i = 0; i < point_count; i++) {
    w1[i] = (foil.foil_shape.points.col(i) - foil.edge.point_le).dot(point_chord.normalized());
    w2[i] = MathUtil::cross2(foil.foil_shape.points.col(i) - foil.edge.point_le, point_chord.normalized());
  }

  w3 = spline::splind(w1, foil.foil_shape.spline_length.head(point_count));
  w4 = spline::splind(w2, foil.foil_shape.spline_length.head(point_count));

  if (is == 1) {
    str = foil.edge.sle + (foil.foil_shape.spline_length[0] - foil.edge.sle) * boundaryLayerWorkflow.lattice.top.transitionLocation;
  } else {
    str = foil.edge.sle + (foil.foil_shape.spline_length[foil.foil_shape.n - 1] - foil.edge.sle) * boundaryLayerWorkflow.lattice.bottom.transitionLocation;
  }
  str = spline::sinvrt(str, boundaryLayerWorkflow.lattice.get(is).transitionLocation, w1, w3, foil.foil_shape.spline_length.head(point_count), point_count);
  xiforc = std::min((str - sst), boundaryLayerWorkflow.lattice.get(is).arcLengthCoordinates[boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex]);
  if (xiforc < 0.0) {
    ss << " ***  stagnation point is past trip on side " << is << "\n";
    writeString(ss.str());

    xiforc = boundaryLayerWorkflow.lattice.get(is).arcLengthCoordinates[boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex];
  }

  return xiforc;
}
