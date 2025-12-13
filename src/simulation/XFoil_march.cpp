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
      xfoil.flowRegime,
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

  flowRegime = boundaryLayerWorkflow.determineRegimeForStation(side, ibl, ctx.simi, ctx.wake);

  return ctx;
}

void XFoil::checkTransitionIfNeeded(int side, int ibl, bool skipCheck,
                                    int laminarAdvance, double& ami) {
  if (skipCheck || flowRegime == FlowRegimeEnum::Turbulent || flowRegime == FlowRegimeEnum::Wake) {
    return;
  }

  trchek();
  ami = boundaryLayerWorkflow.state.station2.param.amplz;
  if (flowRegime == FlowRegimeEnum::Transition) {
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
           (0.03 * (xt.scalar - boundaryLayerWorkflow.state.station1.param.xz) - 0.15 * (boundaryLayerWorkflow.state.station2.param.xz - xt.scalar)) /
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

    //**** sweep downstream setting up bl equation linearizations
    for (int ibl = 0; ibl < boundaryLayerWorkflow.lattice.get(is).stationCount - 1; ++ibl) {
      
      int iv = boundaryLayerWorkflow.lattice.get(is).stationToSystem[ibl];

      const bool stationIsSimilarity = (ibl == 0);
      const bool stationIsWake = (ibl > boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex);
      const bool stationIsTransitionCandidate = (ibl == output.itran.get(is));
      const bool stationIsDownstreamOfTransition = (ibl > output.itran.get(is));
      output.flowRegime =
          boundaryLayerWorkflow.determineRegimeForStation(is, ibl,
                                                          stationIsSimilarity,
                                                          stationIsWake);

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

      if (stationIsWake) {
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
                     aerodynamicCache.dij(boundaryLayerWorkflow.lattice.get(is).stationToPanel[ibl], boundaryLayerWorkflow.lattice.get(js).stationToPanel[jbl]);
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

      //---- check for transition and set xt, etc. if found
      if (stationIsTransitionCandidate) {
        trchek();
        ami = boundaryLayerWorkflow.state.station2.param.amplz;
      }

      if (stationIsTransitionCandidate && flowRegime != FlowRegimeEnum::Transition) {
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
      if (flowRegime == FlowRegimeEnum::Transition) {
        //------ save transition location
        output.itran.get(is) = ibl;
        output.flowRegime = FlowRegimeEnum::Turbulent;
      }

      if (ibl == boundaryLayerWorkflow.lattice.get(is).trailingEdgeIndex) {
        //----- set "2" variables at te to output.wake correlations for next station

        output.flowRegime = FlowRegimeEnum::Wake;
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
  return boundaryLayerWorkflow.trchek(*this);
}


bool XFoil::trdif() {
  return boundaryLayerWorkflow.trdif(*this);
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
                        aerodynamicCache.dij(boundaryLayerWorkflow.lattice.get(is).stationToPanel[ibl],
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
