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

/** ----------------------------------------------------
 *      marches the bls and wake in mixed mode using
 *      the current ue and hk.  the calculated ue
 *      and hk lie along a line quasi-normal to the
 *      natural ue-hk characteristic line of the
 *      current bl so that the goldstein or levy-lees
 *      singularity is never encountered.  continuous
 *      checking of transition onset is performed.
 * ----------------------------------------------------- */
bool XFoil::mrchdu(BoundaryLayerState& state,
                   [[maybe_unused]] BoundaryLayerLattice& lattice) {
  const double deps = 0.000005;
  const double senswt = 1000.0;

  double sens = 0.0;
  double sennew = 0.0;
  double ami = 0.0;

  for (int is = 1; is <= 2; ++is) {
    xiforc = xifset(is);
    int itrold = boundaryLayerWorkflow.lattice.transitionIndex.get(is);

    tran = false;
    turb = false;
    boundaryLayerWorkflow.lattice.transitionIndex.get(is) = boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is);

    for (int ibl = 0; ibl < boundaryLayerWorkflow.lattice.stationCount.get(is) - 1; ++ibl) {
      MixedModeStationContext ctx = prepareMixedModeStation(is, ibl, itrold, ami);
      bool converged = performMixedModeNewtonIteration(is, ibl, itrold, ctx, deps,
                                                       senswt, sens, sennew, ami);
      if (!converged) {
        handleMixedModeNonConvergence(is, ibl, ctx, ami);
      }

      sens = sennew;

      if (ibl < boundaryLayerWorkflow.lattice.transitionIndex.get(is)) {
        boundaryLayerWorkflow.lattice.ctau.get(is)[ibl] = ctx.ami;
      } else {
        boundaryLayerWorkflow.lattice.ctau.get(is)[ibl] = ctx.cti;
      }
      boundaryLayerWorkflow.lattice.thet.get(is)[ibl] = ctx.thi;
      boundaryLayerWorkflow.lattice.dstr.get(is)[ibl] = ctx.dsi;
      boundaryLayerWorkflow.lattice.uedg.get(is)[ibl] = ctx.uei;
      boundaryLayerWorkflow.lattice.mass.get(is)[ibl] = ctx.dsi * ctx.uei;
      boundaryLayerWorkflow.lattice.ctq.get(is)[ibl] = boundaryLayerWorkflow.state.station2.cqz.scalar;

      {
        blData updatedCurrent =
            blprv(state.current(), ctx.xsi, ctx.ami, ctx.cti, ctx.thi,
                  ctx.dsi, ctx.dswaki, ctx.uei);
        state.current() = updatedCurrent;
      }
      blkin(state);

      state.stepbl();

      if (tran || ibl == boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is)) {
        turb = true;
      }

      tran = false;

      if (isCancelled()) {
        return false;
      }
    }
  }
  return true;
}

bool XFoil::mrchdu() { return mrchdu(boundaryLayerWorkflow.state, boundaryLayerWorkflow.lattice); }

XFoil::MixedModeStationContext XFoil::prepareMixedModeStation(int side, int ibl,
                                                              int itrold,
                                                              double& ami) {
  MixedModeStationContext ctx;

  ctx.simi = (ibl == 0);
  ctx.wake = ibl > boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(side);
  ctx.xsi = boundaryLayerWorkflow.lattice.xssi.get(side)[ibl];
  ctx.uei = boundaryLayerWorkflow.lattice.uedg.get(side)[ibl];
  ctx.thi = boundaryLayerWorkflow.lattice.thet.get(side)[ibl];
  ctx.dsi = boundaryLayerWorkflow.lattice.dstr.get(side)[ibl];

  if (ibl < itrold) {
    ami = boundaryLayerWorkflow.lattice.ctau.get(side)[ibl];
    ctx.cti = 0.03;
  } else {
    ctx.cti = boundaryLayerWorkflow.lattice.ctau.get(side)[ibl];
    if (ctx.cti <= 0.0) {
      ctx.cti = 0.03;
    }
  }
  ctx.ami = ami;

  if (ctx.wake) {
    int iw = ibl - boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(side);
    ctx.dswaki = wgap[iw - 1];
  } else {
    ctx.dswaki = 0.0;
  }

  double thickness_limit = (ibl <= boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(side)) ? 1.02 : 1.00005;
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
    boundaryLayerWorkflow.lattice.transitionIndex.get(side) = ibl;
  } else {
    boundaryLayerWorkflow.lattice.transitionIndex.get(side) = ibl + laminarAdvance;
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
      blData updatedCurrent = blprv(boundaryLayerWorkflow.state.current(), ctx.xsi, ami,
                                    ctx.cti, ctx.thi, ctx.dsi, ctx.dswaki,
                                    ctx.uei);
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
      boundaryLayerWorkflow.configureSimilarityRow(*this, ueref);
    } else {
      const bool resetSensitivity = (itbl <= 5);
      const bool averageSensitivity = (itbl > 5 && itbl <= 15);
      boundaryLayerWorkflow.configureViscousRow(
          *this, hkref, ueref, senswt, resetSensitivity, averageSensitivity,
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

  if (ctx.dmax > 0.1 && ibl >= 2) {
    if (ibl <= boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(side)) {
      ctx.thi = boundaryLayerWorkflow.lattice.thet.get(side)[ibl - 1] *
                sqrt(boundaryLayerWorkflow.lattice.xssi.get(side)[ibl] / boundaryLayerWorkflow.lattice.xssi.get(side)[ibl - 1]);
      ctx.dsi = boundaryLayerWorkflow.lattice.dstr.get(side)[ibl - 1] *
                sqrt(boundaryLayerWorkflow.lattice.xssi.get(side)[ibl] / boundaryLayerWorkflow.lattice.xssi.get(side)[ibl - 1]);
      ctx.uei = boundaryLayerWorkflow.lattice.uedg.get(side)[ibl - 1];
    } else {
      if (ibl == boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(side) + 1) {
        ctx.cti = ctx.cte;
        ctx.thi = ctx.tte;
        ctx.dsi = ctx.dte;
        ctx.uei = boundaryLayerWorkflow.lattice.uedg.get(side)[ibl - 1];
      } else {
        ctx.thi = boundaryLayerWorkflow.lattice.thet.get(side)[ibl - 1];
        double ratlen = (boundaryLayerWorkflow.lattice.xssi.get(side)[ibl] - boundaryLayerWorkflow.lattice.xssi.get(side)[ibl - 1]) /
                        (10.0 * boundaryLayerWorkflow.lattice.dstr.get(side)[ibl - 1]);
        ctx.dsi = (boundaryLayerWorkflow.lattice.dstr.get(side)[ibl - 1] + ctx.thi * ratlen) /
                  (1.0 + ratlen);
        ctx.uei = boundaryLayerWorkflow.lattice.uedg.get(side)[ibl - 1];
      }
    }

    if (ibl == boundaryLayerWorkflow.lattice.transitionIndex.get(side)) {
      ctx.cti = 0.05;
    }
    if (ibl > boundaryLayerWorkflow.lattice.transitionIndex.get(side)) {
      ctx.cti = boundaryLayerWorkflow.lattice.ctau.get(side)[ibl - 1];
    }
  }

  {
    blData updatedCurrent = blprv(boundaryLayerWorkflow.state.current(), ctx.xsi, ami,
                                  ctx.cti, ctx.thi, ctx.dsi, ctx.dswaki,
                                  ctx.uei);
    boundaryLayerWorkflow.state.current() = updatedCurrent;
  }
  blkin(boundaryLayerWorkflow.state);

  checkTransitionIfNeeded(side, ibl, ctx.simi, 2, ami);

  if (ibl < boundaryLayerWorkflow.lattice.transitionIndex.get(side)) {
    boundaryLayerWorkflow.state.station2 = blvar(boundaryLayerWorkflow.state.station2, FlowRegimeEnum::Laminar);
    blmid(boundaryLayerWorkflow.state, FlowRegimeEnum::Laminar);
  }
  if (ibl >= boundaryLayerWorkflow.lattice.transitionIndex.get(side)) {
    boundaryLayerWorkflow.state.station2 = blvar(boundaryLayerWorkflow.state.station2, FlowRegimeEnum::Turbulent);
    blmid(boundaryLayerWorkflow.state, FlowRegimeEnum::Turbulent);
  }
  if (ctx.wake) {
    boundaryLayerWorkflow.state.station2 = blvar(boundaryLayerWorkflow.state.station2, FlowRegimeEnum::Wake);
    blmid(boundaryLayerWorkflow.state, FlowRegimeEnum::Wake);
  }

  ctx.ami = ami;
}

double XFoil::calcHtarg(int ibl, int is, bool wake) {
  if (ibl < boundaryLayerWorkflow.lattice.transitionIndex.get(is)) {
    return boundaryLayerWorkflow.state.station1.hkz.scalar +
           0.03 * (boundaryLayerWorkflow.state.station2.param.xz - boundaryLayerWorkflow.state.station1.param.xz) / boundaryLayerWorkflow.state.station1.param.tz;
  } else if (ibl == boundaryLayerWorkflow.lattice.transitionIndex.get(is)) {
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


/** ----------------------------------------------------
 *      marches the bls and wake in direct mode using
 *      the uedg array. if separation is encountered,
 *      a plausible value of hk extrapolated from
 *      upstream is prescribed instead.  continuous
 *      checking of transition onset is performed.
 * ----------------------------------------------------*/
bool XFoil::mrchue() {
  std::stringstream ss;
  bool direct;

  double msq, ratlen, dsw, hklim;
  double xsi, uei, ucon, tsq, thi, ami, cti, dsi;
  double dswaki;
  double htest, hktest;
  double cte, dte, tte, dmax, hmax, htarg = 0.0;
  // double cte = dte = tte = dmax = hmax = htarg = 0.0;

  //---- shape parameters for separation criteria
  const double hlmax = 3.8;
  const double htmax = 2.5;

  for (int is = 1; is <= 2; is++) { // 2000
    ss << "    Side " << is << " ...\n";
    writeString(ss.str());

    //---- set forced transition arc length position
    xiforc = xifset(is);

    //---- initialize similarity station with thwaites' formula (keep ibl=1)
    // xssi is 0-based: first BL station arc length at index 0
    xsi = boundaryLayerWorkflow.lattice.xssi.get(is)[0];
    uei = boundaryLayerWorkflow.lattice.uedg.get(is)[0];

    ucon = uei / xsi;
    tsq = 0.45 / (ucon * 6.0 * reybl);
    thi = sqrt(tsq);
    dsi = 2.2 * thi;
    ami = 0.0;

    //---- initialize ctau for first turbulent station
    cti = 0.03;

    tran = false;
    turb = false;
    boundaryLayerWorkflow.lattice.transitionIndex.get(is) = boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is);

    //---- march downstream
    for (int ibl = 0; ibl < boundaryLayerWorkflow.lattice.stationCount.get(is) - 1; ++ibl) {
      //int ibl = ibl - 1;
      int ibm = ibl + 1;
      int iw = ibl - boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is);
      simi = (ibl == 0);
      wake = ibl > boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is);

      //------ prescribed quantities (xssi now 0-based)
      xsi = boundaryLayerWorkflow.lattice.xssi.get(is)[ibl];
      uei = boundaryLayerWorkflow.lattice.uedg.get(is)[ibl];

      if (wake) {
        iw = ibl - boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is);
        dswaki = wgap[iw - 1];
      } else
        dswaki = 0.0;

      direct = true;
      bool converged = false;

      //------ newton iteration loop for current station
      for (int itbl = 1; itbl <= 25; itbl++) { // 100

        //-------- assemble 10x3 linearized system for dctau, dth, dds, due, dxi
        //         at the previous "1" station and the current "2" station
        //         (the "1" station coefficients will be ignored)

        {
          blData updatedCurrent =
              blprv(boundaryLayerWorkflow.state.current(), xsi, ami, cti, thi, dsi,
                    dswaki, uei);
          boundaryLayerWorkflow.state.current() = updatedCurrent;
        }
        blkin(boundaryLayerWorkflow.state);

        //-------- check for transition and set appropriate flags and things
        if ((!simi) && (!turb)) {
          trchek();
          ami = boundaryLayerWorkflow.state.station2.param.amplz;

          //--------- fixed bug   md 7 jun 99
          if (tran) {
            boundaryLayerWorkflow.lattice.transitionIndex.get(is) = ibl;
            if (cti <= 0.0) {
              cti = 0.03;
              boundaryLayerWorkflow.state.station2.param.sz = cti;
            }
          } else
            boundaryLayerWorkflow.lattice.transitionIndex.get(is) = ibl + 2;
        }

        if (ibl == boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is) + 1) {
          tte = boundaryLayerWorkflow.lattice.thet.get(1)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top] +
                boundaryLayerWorkflow.lattice.thet.get(2)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom];
          dte = boundaryLayerWorkflow.lattice.dstr.get(1)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top] +
                boundaryLayerWorkflow.lattice.dstr.get(2)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom] + foil.edge.ante;
          cte = (boundaryLayerWorkflow.lattice.ctau.get(1)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top] *
                     boundaryLayerWorkflow.lattice.thet.get(1)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top] +
                 boundaryLayerWorkflow.lattice.ctau.get(2)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom] *
                     boundaryLayerWorkflow.lattice.thet.get(2)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom]) /
                tte;
          boundaryLayerWorkflow.tesys(*this, cte, tte, dte);
        } else
          blsys(boundaryLayerWorkflow.state, boundaryLayerWorkflow.lattice);

        if (direct) {
          //--------- try direct mode (set due = 0 in currently empty 4th line)
          blc.a2(3, 0) = 0.0;
          blc.a2(3, 1) = 0.0;
          blc.a2(3, 2) = 0.0;
          blc.a2(3, 3) = 1.0;
          blc.rhs[3] = 0.0;
          //--------- solve newton system for current "2" station
          blc.rhs = blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);
          //--------- determine max changes and underrelax if necessary
          dmax = std::max(fabs(blc.rhs[1] / thi), fabs(blc.rhs[2] / dsi));
          if (ibl < boundaryLayerWorkflow.lattice.transitionIndex.get(is))
            dmax = std::max(dmax, fabs(blc.rhs[0] / 10.0));
          if (ibl >= boundaryLayerWorkflow.lattice.transitionIndex.get(is))
            dmax = std::max(dmax, fabs(blc.rhs[0] / cti));

          rlx = 1.0;
          if (dmax > 0.3)
            rlx = 0.3 / dmax;
          //--------- see if direct mode is not applicable
          if (ibl != boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is) + 1) {
            //---------- calculate resulting kinematic shape parameter hk
            msq =
                uei * uei * hstinv / (gm1bl * (1.0 - 0.5 * uei * uei * hstinv));
            htest = (dsi + rlx * blc.rhs[2]) / (thi + rlx * blc.rhs[1]);
            boundary_layer::KineticShapeParameterResult hkin_result =
                boundary_layer::hkin(htest, msq);
            hktest = hkin_result.hk;

            //---------- decide whether to do direct or inverse problem based on
            // hk
            if (ibl < boundaryLayerWorkflow.lattice.transitionIndex.get(is))
              hmax = hlmax;
            if (ibl >= boundaryLayerWorkflow.lattice.transitionIndex.get(is))
              hmax = htmax;
            direct = (hktest < hmax);
          }
          if (direct) {
            //---------- update as usual
            if (ibl >= boundaryLayerWorkflow.lattice.transitionIndex.get(is))
              cti = cti + rlx * blc.rhs[0];
            thi = thi + rlx * blc.rhs[1];
            dsi = dsi + rlx * blc.rhs[2];
          } else {
            //---------- set prescribed hk for inverse calculation at the
            // current station
            htarg = calcHtarg(ibl, is, wake);

            //---------- limit specified hk to something reasonable
            if (wake)
              htarg = std::max(htarg, 1.01);
            else
              htarg = std::max(htarg, hmax);

            ss.str("");
            ss << "     mrchue: inverse mode at " << ibl
               << "    hk=" << std::fixed << std::setprecision(3) << htarg << "\n";
            writeString(ss.str());
            ss.str("");

            //---------- try again with prescribed hk

            continue;
          }
        } else {
          //-------- inverse mode (force hk to prescribed value htarg)
          blc.a2(3, 0) = 0.0;
          blc.a2(3, 1) = boundaryLayerWorkflow.state.station2.hkz.t();
          blc.a2(3, 2) = boundaryLayerWorkflow.state.station2.hkz.d();
          blc.a2(3, 3) = boundaryLayerWorkflow.state.station2.hkz.u();
          blc.rhs[3] = htarg - boundaryLayerWorkflow.state.station2.hkz.scalar;
          blc.rhs = blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);

          dmax = std::max(fabs(blc.rhs[1] / thi), fabs(blc.rhs[2] / dsi));
          if (ibl >= boundaryLayerWorkflow.lattice.transitionIndex.get(is))
            dmax = std::max(dmax, fabs(blc.rhs[0] / cti));
          rlx = 1.0;
          if (dmax > 0.3)
            rlx = 0.3 / dmax;
          //--------- update variables
          if (ibl >= boundaryLayerWorkflow.lattice.transitionIndex.get(is))
            cti = cti + rlx * blc.rhs[0];
          thi = thi + rlx * blc.rhs[1];
          dsi = dsi + rlx * blc.rhs[2];
          uei = uei + rlx * blc.rhs[3];
        }
        //-------- eliminate absurd transients

        if (ibl >= boundaryLayerWorkflow.lattice.transitionIndex.get(is)) {
          cti = std::min(cti, 0.30);
          cti = std::max(cti, 0.0000001);
        }
        if (ibl <= boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is))
          hklim = 1.02;
        else
          hklim = 1.00005;
        msq = uei * uei * hstinv / (gm1bl * (1.0 - 0.5 * uei * uei * hstinv));
        dsw = dsi - dswaki;
        dslim(dsw, thi, msq, hklim);
        dsi = dsw + dswaki;
        if (dmax <= 0.00001) {
          converged = true;
          break;
        }

      } // end itbl loop

      if (!converged) {
        ss.str("");
        ss << "     mrchue: convergence failed at " << ibl << ",  side " << is
           << ", res =" << std::fixed << std::setprecision(3) << dmax << "\n";
        writeString(ss.str());

        //------ the current unconverged solution might still be reasonable...
        if (dmax > 0.1) {
          //------- the current solution is garbage --> extrapolate values
          // instead
          if (ibl >= 2) {
            if (ibl <= boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is)) {
              thi = boundaryLayerWorkflow.lattice.thet.get(is)[ibl - 1] *
                    sqrt(boundaryLayerWorkflow.lattice.xssi.get(is)[ibl] / boundaryLayerWorkflow.lattice.xssi.get(is)[ibl - 1]);
              dsi = boundaryLayerWorkflow.lattice.dstr.get(is)[ibl - 1] *
                    sqrt(boundaryLayerWorkflow.lattice.xssi.get(is)[ibl] / boundaryLayerWorkflow.lattice.xssi.get(is)[ibl - 1]);
            } else {
              if (ibl == boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is) + 1) {
                cti = cte;
                thi = tte;
                dsi = dte;
              } else {
                thi = boundaryLayerWorkflow.lattice.thet.get(is)[ibl - 1];
                ratlen = (boundaryLayerWorkflow.lattice.xssi.get(is)[ibl] - boundaryLayerWorkflow.lattice.xssi.get(is)[ibl - 1]) /
                         (10.0 * boundaryLayerWorkflow.lattice.dstr.get(is)[ibl - 1]);
                dsi = (boundaryLayerWorkflow.lattice.dstr.get(is)[ibl - 1] +
                       thi * ratlen) /
                      (1.0 + ratlen);
              }
            }
            if (ibl == boundaryLayerWorkflow.lattice.transitionIndex.get(is))
              cti = 0.05;
            if (ibl > boundaryLayerWorkflow.lattice.transitionIndex.get(is))
              cti = boundaryLayerWorkflow.lattice.ctau.get(is)[ibl - 1];

            uei = boundaryLayerWorkflow.lattice.uedg.get(is)[ibl];

            if (ibl < boundaryLayerWorkflow.lattice.stationCount.get(is) - 1)
              uei = 0.5 *
                    (boundaryLayerWorkflow.lattice.uedg.get(is)[ibl - 1] +
                     boundaryLayerWorkflow.lattice.uedg.get(is)[ibl + 1]);
          }
        }
        // 109
        {
          blData updatedCurrent =
              blprv(boundaryLayerWorkflow.state.current(), xsi, ami, cti, thi, dsi,
                    dswaki, uei);
          boundaryLayerWorkflow.state.current() = updatedCurrent;
        }
        blkin(boundaryLayerWorkflow.state);
        //------- check for transition and set appropriate flags and things
        if ((!simi) && (!turb)) {
          trchek();
          ami = boundaryLayerWorkflow.state.station2.param.amplz;
          if (tran)
            boundaryLayerWorkflow.lattice.transitionIndex.get(is) = ibl;
          if (!tran)
            boundaryLayerWorkflow.lattice.transitionIndex.get(is) = ibl + 2;
        }
        //------- set all other extrapolated values for current station
        if (ibl < boundaryLayerWorkflow.lattice.transitionIndex.get(is))
          boundaryLayerWorkflow.state.station2 = blvar(boundaryLayerWorkflow.state.station2, FlowRegimeEnum::Laminar);
        if (ibl >= boundaryLayerWorkflow.lattice.transitionIndex.get(is))
          boundaryLayerWorkflow.state.station2 = blvar(boundaryLayerWorkflow.state.station2, FlowRegimeEnum::Turbulent);
        if (wake)
          boundaryLayerWorkflow.state.station2 = blvar(boundaryLayerWorkflow.state.station2, FlowRegimeEnum::Wake);
        if (ibl < boundaryLayerWorkflow.lattice.transitionIndex.get(is))
          blmid(boundaryLayerWorkflow.state, FlowRegimeEnum::Laminar);
        if (ibl >= boundaryLayerWorkflow.lattice.transitionIndex.get(is))
          blmid(boundaryLayerWorkflow.state, FlowRegimeEnum::Turbulent);
        if (wake)
          blmid(boundaryLayerWorkflow.state, FlowRegimeEnum::Wake);
      }
      //------ store primary variables
      if (ibl < boundaryLayerWorkflow.lattice.transitionIndex.get(is))
        boundaryLayerWorkflow.lattice.ctau.get(is)[ibl] = ami;
      if (ibl >= boundaryLayerWorkflow.lattice.transitionIndex.get(is))
        boundaryLayerWorkflow.lattice.ctau.get(is)[ibl] = cti;
      boundaryLayerWorkflow.lattice.thet.get(is)[ibl] = thi;
      boundaryLayerWorkflow.lattice.dstr.get(is)[ibl] = dsi;
      boundaryLayerWorkflow.lattice.uedg.get(is)[ibl] = uei;
      boundaryLayerWorkflow.lattice.mass.get(is)[ibl] = dsi * uei;
      boundaryLayerWorkflow.lattice.ctq.get(is)[ibl] = boundaryLayerWorkflow.state.station2.cqz.scalar;

      //------ set "1" variables to "2" variables for next streamwise station
      {
        blData updatedCurrent =
            blprv(boundaryLayerWorkflow.state.current(), xsi, ami, cti, thi, dsi,
                  dswaki, uei);
        boundaryLayerWorkflow.state.current() = updatedCurrent;
      }
      blkin(boundaryLayerWorkflow.state);

      boundaryLayerWorkflow.state.stepbl();

      //------ turbulent intervals will follow transition interval or te
      if (tran || ibl == boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is)) {
        turb = true;
      }

      tran = false;

      if (ibl == boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is)) {
        thi = boundaryLayerWorkflow.lattice.thet.get(1)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top] +
              boundaryLayerWorkflow.lattice.thet.get(2)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom];
        dsi = boundaryLayerWorkflow.lattice.dstr.get(1)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top] +
              boundaryLayerWorkflow.lattice.dstr.get(2)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom] + foil.edge.ante;
      }
    }
  }
  return true;
}

SetblInputView XFoil::makeSetblInputView() const {
  return SetblInputView{lblini,
                        boundaryLayerWorkflow.lattice.uedg,
                        boundaryLayerWorkflow.lattice.ctau,
                        boundaryLayerWorkflow.lattice.thet,
                        boundaryLayerWorkflow.lattice.dstr,
                        boundaryLayerWorkflow.lattice.mass,
                        boundaryLayerWorkflow.lattice.ctq,
                        boundaryLayerWorkflow.lattice.transitionIndex};
}

SetblOutputView XFoil::makeSetblOutputView() {
  return SetblOutputView{lblini,
                         gm1bl,
                         qinfbl,
                         tkbl,
                         tkbl_ms,
                         rstbl,
                         rstbl_ms,
                         hstinv,
                         hstinv_ms,
                         reybl,
                         reybl_re,
                         reybl_ms,
                         amcrit,
                         boundaryLayerWorkflow.lattice.uedg,
                         boundaryLayerWorkflow.lattice.ctau,
                         boundaryLayerWorkflow.lattice.thet,
                         boundaryLayerWorkflow.lattice.dstr,
                         boundaryLayerWorkflow.lattice.mass,
                         boundaryLayerWorkflow.lattice.ctq,
                         boundaryLayerWorkflow.lattice.transitionIndex,
                         va,
                         vb,
                         vdel,
                         vm,
                         vz,
                         tran,
                         turb,
                         wake,
                         simi,
                         xiforc};
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
  if (lalfa)
    clmr = cl;
  else
    clmr = clspec;

  cti = 0.0; // techwinder added, otherwise variable is not initialized

  //---- set current minf(cl)
  ma_clmr = getActualMach(clmr, mach_type);
  re_clmr = getActualReynolds(clmr, reynolds_type);

  msq_clmr = 2.0 * minf * ma_clmr;

  //---- set compressibility parameter tklam and derivative tk_msq
  comset();

  //---- set gas constant (= cp/cv)
  output.gm1bl = gamm1;

  //---- set parameters for compressibility correction
  output.qinfbl = qinf;
  output.tkbl = tklam;
  output.tkbl_ms = tkl_msq;

  //---- stagnation density and 1/enthalpy
  output.rstbl = pow((1.0 + 0.5 * output.gm1bl * minf * minf), (1.0 / output.gm1bl));
  output.rstbl_ms = 0.5 * output.rstbl / (1.0 + 0.5 * output.gm1bl * minf * minf);
  output.hstinv = output.gm1bl * (minf / output.qinfbl) * (minf / output.qinfbl) /
           (1.0 + 0.5 * output.gm1bl * minf * minf);
  output.hstinv_ms = output.gm1bl * (1.0 / output.qinfbl) * (1.0 / output.qinfbl) /
                  (1.0 + 0.5 * output.gm1bl * minf * minf) -
              0.5 * output.gm1bl * output.hstinv / (1.0 + 0.5 * output.gm1bl * minf * minf);

  //---- set reynolds number based on freestream density, velocity, viscosity
  herat = 1.0 - 0.5 * output.qinfbl * output.qinfbl * output.hstinv;
  herat_ms = -0.5 * output.qinfbl * output.qinfbl * output.hstinv_ms;

  output.reybl = reinf * sqrt(herat * herat * herat) * (1.0 + hvrat) / (herat + hvrat);
  output.reybl_re = sqrt(herat * herat * herat) * (1.0 + hvrat) / (herat + hvrat);
  output.reybl_ms = output.reybl * (1.5 / herat - 1.0 / (herat + hvrat)) * herat_ms;

  output.amcrit = acrit;

  if (!input.lblini) {
    //----- initialize bl by marching with ue (fudge at separation)
    // TRACE(" initializing bl ...\n");
    writeString("   Initializing bl ...\n");

    mrchue();
    output.lblini = true;
  }

  //---- march bl with current ue and ds to establish transition
  mrchdu();

  SidePair<VectorXd> usav = input.uedg;

  ueset();
  const auto swapped_edge_velocities = swapEdgeVelocities(usav);
  usav = swapped_edge_velocities.swappedUsav;
  output.uedg = swapped_edge_velocities.restoredUedg;
  jvte1 = boundaryLayerWorkflow.lattice.stationToSystem.top[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top];
  jvte2 = boundaryLayerWorkflow.lattice.stationToSystem.bottom[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom];

  dule1 = output.uedg.top[0] - usav.top[0];
  dule2 = output.uedg.bottom[0] - usav.bottom[0];

  //---- set le and te ue sensitivities wrt all m values
  const auto le_te_sensitivities = computeLeTeSensitivities(
      boundaryLayerWorkflow.lattice.stationToPanel.get(1)[0], boundaryLayerWorkflow.lattice.stationToPanel.get(2)[0], boundaryLayerWorkflow.lattice.stationToPanel.get(1)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top],
      boundaryLayerWorkflow.lattice.stationToPanel.get(2)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom]);
  ule1_m = le_te_sensitivities.ule1_m;
  ule2_m = le_te_sensitivities.ule2_m;
  ute1_m = le_te_sensitivities.ute1_m;
  ute2_m = le_te_sensitivities.ute2_m;

  ule1_a = boundaryLayerWorkflow.lattice.uinv_a.get(1)[0];
  ule2_a = boundaryLayerWorkflow.lattice.uinv_a.get(2)[0];

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
    for (int ibl = 0; ibl < boundaryLayerWorkflow.lattice.stationCount.get(is) - 1; ++ibl) {
      
      int iv = boundaryLayerWorkflow.lattice.stationToSystem.get(is)[ibl];

      output.simi = (ibl == 0);
      output.wake = (ibl > boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is));
      output.tran = (ibl == output.itran.get(is));
      output.turb = (ibl > output.itran.get(is));

      //---- set primary variables for current station
      xsi = boundaryLayerWorkflow.lattice.xssi.get(is)[ibl];
      if (ibl < output.itran.get(is))
        ami = output.ctau.get(is)[ibl];
      else
        cti = output.ctau.get(is)[ibl];
      uei = output.uedg.get(is)[ibl];
      thi = output.thet.get(is)[ibl];
      mdi = output.mass.get(is)[ibl];

      dsi = mdi / uei;

      if (output.wake) {
        int iw = ibl - boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is);
        dswaki = wgap[iw - 1];
      } else
        dswaki = 0.0;

      //---- set derivatives of dsi (= d2)
      d2_m2 = 1.0 / uei;
      d2_u2 = -dsi / uei;

      for (int js = 1; js <= 2; js++) {
        for (int jbl = 0; jbl < boundaryLayerWorkflow.lattice.stationCount.get(js) - 1; ++jbl) {
          int jv = boundaryLayerWorkflow.lattice.stationToSystem.get(js)[jbl];
          u2_m[jv] = -boundaryLayerWorkflow.lattice.vti.get(is)[ibl] * boundaryLayerWorkflow.lattice.vti.get(js)[jbl] *
                     dij(boundaryLayerWorkflow.lattice.stationToPanel.get(is)[ibl], boundaryLayerWorkflow.lattice.stationToPanel.get(js)[jbl]);
          d2_m[jv] = d2_u2 * u2_m[jv];
        }
      }
      d2_m[iv] = d2_m[iv] + d2_m2;

      u2_a = boundaryLayerWorkflow.lattice.uinv_a.get(is)[ibl];
      d2_a = d2_u2 * u2_a;

  //---- "forced" changes due to mismatch between edge velocities and
  // usav=uinv+dij*output.mass
      due2 = output.uedg.get(is)[ibl] - usav.get(is)[ibl];
      dds2 = d2_u2 * due2;

      {
        blData updatedCurrent =
            blprv(boundaryLayerWorkflow.state.current(), xsi, ami, cti, thi, dsi,
                  dswaki, uei);
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

      //---- assemble 10x4 linearized system for dctau, dth, dds, due, dxi
      //	   at the previous "1" station and the current "2" station

      if (ibl == boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is) + 1) {
        //----- define quantities at start of output.wake, adding te base thickness to
        // dstar
        tte = output.thet.get(1)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top] +
              output.thet.get(2)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom];
        dte = output.dstr.get(1)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top] +
              output.dstr.get(2)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom] + foil.edge.ante;
        cte = (output.ctau.get(1)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top] *
                   output.thet.get(1)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top] +
               output.ctau.get(2)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom] *
                   output.thet.get(2)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom]) /
               tte;
        boundaryLayerWorkflow.tesys(*this, cte, tte, dte);

        tte_tte1 = 1.0;
        tte_tte2 = 1.0;
        dte_mte1 = 1.0 / output.uedg.top[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top];
        dte_ute1 = -output.dstr.get(1)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top] /
                    output.uedg.top[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top];
        dte_mte2 = 1.0 / output.uedg.bottom[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom];
        dte_ute2 = -output.dstr.get(2)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom] /
                    output.uedg.bottom[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom];
        cte_cte1 = output.thet.get(1)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top] / tte;
        cte_cte2 = output.thet.get(2)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom] / tte;
        cte_tte1 = (output.ctau.get(1)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top] - cte) / tte;
        cte_tte2 = (output.ctau.get(2)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom] - cte) / tte;

        //----- re-define d1 sensitivities wrt m since d1 depends on both te ds
        // values
      for (int js = 1; js <= 2; js++) {
        for (int jbl = 0; jbl < boundaryLayerWorkflow.lattice.stationCount.get(js) - 1; ++jbl) {
            int jv = boundaryLayerWorkflow.lattice.stationToSystem.get(js)[jbl];
            d1_m[jv] = dte_ute1 * ute1_m[jv] + dte_ute2 * ute2_m[jv];
          }
        }
        d1_m[jvte1] = d1_m[jvte1] + dte_mte1;
        d1_m[jvte2] = d1_m[jvte2] + dte_mte2;

        //----- "forced" changes from  output.uedg --- usav=uinv+dij*output.mass	mismatch
        due1 = 0.0;
        dds1 =
            dte_ute1 *
                (output.uedg.top[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top] -
                 usav.top[boundaryLayerWorkflow.lattice.trailingEdgeIndex.top]) +
            dte_ute2 *
                (output.uedg.bottom[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom] -
                 usav.bottom[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom]);
      } else {
        blsys(boundaryLayerWorkflow.state, boundaryLayerWorkflow.lattice);
      }

      //---- save wall shear and equil. max shear coefficient for plotting
      // output
      output.ctq.get(is)[ibl] = boundaryLayerWorkflow.state.station2.cqz.scalar;

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
        output.vm[0][jv][iv] = blc.a1(0, 2) * d1_m[jv] + blc.a1(0, 3) * u1_m[jv] +
                        blc.a2(0, 2) * d2_m[jv] + blc.a2(0, 3) * u2_m[jv] +
                        (blc.a1(0, 4) + blc.a2(0, 4) + blc.d_xi[0]) *
                            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
      }

      output.vb[iv](0, 0) = blc.a1(0, 0);
      output.vb[iv](0, 1) = blc.a1(0, 1);

      output.va[iv](0, 0) = blc.a2(0, 0);
      output.va[iv](0, 1) = blc.a2(0, 1);

      if (lalfa)
        output.vdel[iv](0, 1) = blc.d_re[0] * re_clmr + blc.d_msq[0] * msq_clmr;
      else
        output.vdel[iv](0, 1) = (blc.a1(0, 3) * u1_a + blc.a1(0, 2) * d1_a) +
                         (blc.a2(0, 3) * u2_a + blc.a2(0, 2) * d2_a) +
                         (blc.a1(0, 4) + blc.a2(0, 4) + blc.d_xi[0]) *
                             (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

      output.vdel[iv](0, 0) = blc.rhs[0] + (blc.a1(0, 3) * due1 + blc.a1(0, 2) * dds1) +
                       (blc.a2(0, 3) * due2 + blc.a2(0, 2) * dds2) +
                       (blc.a1(0, 4) + blc.a2(0, 4) + blc.d_xi[0]) *
                           (xi_ule1 * dule1 + xi_ule2 * dule2);

      for (int jv = 1; jv <= nsys; jv++) {
        output.vm[1][jv][iv] = blc.a1(1, 2) * d1_m[jv] + blc.a1(1, 3) * u1_m[jv] +
                        blc.a2(1, 2) * d2_m[jv] + blc.a2(1, 3) * u2_m[jv] +
                        (blc.a1(1, 4) + blc.a2(1, 4) + blc.d_xi[1]) *
                            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
      }
      output.vb[iv](1, 0) = blc.a1(1, 0);
      output.vb[iv](1, 1) = blc.a1(1, 1);

      output.va[iv](1, 0) = blc.a2(1, 0);
      output.va[iv](1, 1) = blc.a2(1, 1);

      if (lalfa)
        output.vdel[iv](1, 1) = blc.d_re[1] * re_clmr + blc.d_msq[1] * msq_clmr;
      else
        output.vdel[iv](1, 1) = (blc.a1(1, 3) * u1_a + blc.a1(1, 2) * d1_a) +
                         (blc.a2(1, 3) * u2_a + blc.a2(1, 2) * d2_a) +
                         (blc.a1(1, 4) + blc.a2(1, 4) + blc.d_xi[1]) *
                             (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

      output.vdel[iv](1, 0) = blc.rhs[1] + (blc.a1(1, 3) * due1 + blc.a1(1, 2) * dds1) +
                       (blc.a2(1, 3) * due2 + blc.a2(1, 2) * dds2) +
                       (blc.a1(1, 4) + blc.a2(1, 4) + blc.d_xi[1]) *
                           (xi_ule1 * dule1 + xi_ule2 * dule2);

      // memory overlap problem
      for (int jv = 1; jv <= nsys; jv++) {
        output.vm[2][jv][iv] = blc.a1(2, 2) * d1_m[jv] + blc.a1(2, 3) * u1_m[jv] +
                        blc.a2(2, 2) * d2_m[jv] + blc.a2(2, 3) * u2_m[jv] +
                        (blc.a1(2, 4) + blc.a2(2, 4) + blc.d_xi[2]) *
                            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
      }

      output.vb[iv](2, 0) = blc.a1(2, 0);
      output.vb[iv](2, 1) = blc.a1(2, 1);

      output.va[iv](2, 0) = blc.a2(2, 0);
      output.va[iv](2, 1) = blc.a2(2, 1);

      if (lalfa)
        output.vdel[iv](2, 1) = blc.d_re[2] * re_clmr + blc.d_msq[2] * msq_clmr;
      else
        output.vdel[iv](2, 1) = (blc.a1(2, 3) * u1_a + blc.a1(2, 2) * d1_a) +
                         (blc.a2(2, 3) * u2_a + blc.a2(2, 2) * d2_a) +
                         (blc.a1(2, 4) + blc.a2(2, 4) + blc.d_xi[2]) *
                             (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

      output.vdel[iv](2, 0) = blc.rhs[2] + (blc.a1(2, 3) * due1 + blc.a1(2, 2) * dds1) +
                       (blc.a2(2, 3) * due2 + blc.a2(2, 2) * dds2) +
                       (blc.a1(2, 4) + blc.a2(2, 4) + blc.d_xi[2]) *
                           (xi_ule1 * dule1 + xi_ule2 * dule2);

      if (ibl == boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is) + 1) {
        //----- redefine coefficients for tte, dte, etc
        output.vz[0][0] = blc.a1(0, 0) * cte_cte1;
        output.vz[0][1] = blc.a1(0, 0) * cte_tte1 + blc.a1(0, 1) * tte_tte1;
        output.vb[iv](0, 0) = blc.a1(0, 0) * cte_cte2;
        output.vb[iv](0, 1) = blc.a1(0, 0) * cte_tte2 + blc.a1(0, 1) * tte_tte2;

        output.vz[1][0] = blc.a1(1, 0) * cte_cte1;
        output.vz[1][1] = blc.a1(1, 0) * cte_tte1 + blc.a1(1, 1) * tte_tte1;
        output.vb[iv](1, 0) = blc.a1(1, 0) * cte_cte2;
        output.vb[iv](1, 1) = blc.a1(1, 0) * cte_tte2 + blc.a1(1, 1) * tte_tte2;

        output.vz[2][0] = blc.a1(2, 0) * cte_cte1;
        output.vz[2][1] = blc.a1(2, 0) * cte_tte1 + blc.a1(2, 1) * tte_tte1;
        output.vb[iv](2, 0) = blc.a1(2, 0) * cte_cte2;
        output.vb[iv](2, 1) = blc.a1(2, 0) * cte_tte2 + blc.a1(2, 1) * tte_tte2;
      }

      //---- turbulent intervals will follow if currently at transition interval
      if (output.tran) {
        output.turb = true;

        //------ save transition location
        output.itran.get(is) = ibl;
      }

      output.tran = false;

      if (ibl == boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is)) {
        //----- set "2" variables at te to output.wake correlations for next station

        output.turb = true;
        output.wake = true;
        boundaryLayerWorkflow.state.station2 = blvar(boundaryLayerWorkflow.state.station2, FlowRegimeEnum::Wake);
        blmid(boundaryLayerWorkflow.state, FlowRegimeEnum::Wake);
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
  boundaryLayerWorkflow.state.station2 = blvar(boundaryLayerWorkflow.state.station2, FlowRegimeEnum::Laminar);

  //---- calculate x1-xt midpoint cfm value
  SkinFrictionCoefficients laminarSkinFriction =
      blmid(boundaryLayerWorkflow.state, FlowRegimeEnum::Laminar);

  //=    at this point, all "2" variables are really "t" variables at xt

  //---- set up newton system for dam, dth, dds, due, dxi  at  x1 and xt
  blc = blDiffSolver.solve(FlowRegimeEnum::Laminar, boundaryLayerWorkflow.state, laminarSkinFriction, amcrit);

  //---- the current newton system is in terms of "1" and "t" variables,
  //-    so calculate its equivalent in terms of "1" and "2" variables.
  //-    in other words, convert residual sensitivities wrt "t" variables
  //-    into sensitivities wrt "1" and "2" variables.  the amplification
  //-    equation is unnecessary here, so the k=1 row is left empty.
  for (int k = 1; k < 3; k++) {
    blrez[k] = blc.rhs[k];
    blm[k] = blc.d_msq[k] + blc.a2(k, 1) * tt_ms + blc.a2(k, 2) * dt_ms +
             blc.a2(k, 3) * ut_ms + blc.a2(k, 4) * xt_ms;
    blr[k] = blc.d_re[k] + blc.a2(k, 1) * tt_re + blc.a2(k, 2) * dt_re +
             blc.a2(k, 3) * ut_re + blc.a2(k, 4) * xt_re;
    blx[k] = blc.d_xi[k] + blc.a2(k, 1) * tt_xf + blc.a2(k, 2) * dt_xf +
             blc.a2(k, 3) * ut_xf + blc.a2(k, 4) * xt_xf;

    bl1(k, 0) = blc.a1(k, 0) + blc.a2(k, 1) * tt_a1 + blc.a2(k, 2) * dt_a1 +
                blc.a2(k, 3) * ut_a1 + blc.a2(k, 4) * xt_a1;
    bl1(k, 1) = blc.a1(k, 1) + blc.a2(k, 1) * tt_t1 + blc.a2(k, 2) * dt_t1 +
                blc.a2(k, 3) * ut_t1 + blc.a2(k, 4) * xt_t1;
    bl1(k, 2) = blc.a1(k, 2) + blc.a2(k, 1) * tt_d1 + blc.a2(k, 2) * dt_d1 +
                blc.a2(k, 3) * ut_d1 + blc.a2(k, 4) * xt_d1;
    bl1(k, 3) = blc.a1(k, 3) + blc.a2(k, 1) * tt_u1 + blc.a2(k, 2) * dt_u1 +
                blc.a2(k, 3) * ut_u1 + blc.a2(k, 4) * xt_u1;
    bl1(k, 4) = blc.a1(k, 4) + blc.a2(k, 1) * tt_x1 + blc.a2(k, 2) * dt_x1 +
                blc.a2(k, 3) * ut_x1 + blc.a2(k, 4) * xt_x1;

    bl2(k, 0) = 0.0;
    bl2(k, 1) = blc.a2(k, 1) * tt_t2 + blc.a2(k, 2) * dt_t2 + blc.a2(k, 3) * ut_t2 +
                blc.a2(k, 4) * xt_t2;
    bl2(k, 2) = blc.a2(k, 1) * tt_d2 + blc.a2(k, 2) * dt_d2 + blc.a2(k, 3) * ut_d2 +
                blc.a2(k, 4) * xt_d2;
    bl2(k, 3) = blc.a2(k, 1) * tt_u2 + blc.a2(k, 2) * dt_u2 + blc.a2(k, 3) * ut_u2 +
                blc.a2(k, 4) * xt_u2;
    bl2(k, 4) = blc.a2(k, 1) * tt_x2 + blc.a2(k, 2) * dt_x2 + blc.a2(k, 3) * ut_x2 +
                blc.a2(k, 4) * xt_x2;
  }

  //**** second, set up turbulent part between xt and x2  ****

  //---- calculate equilibrium shear coefficient cqt at transition point
  boundaryLayerWorkflow.state.station2 = blvar(boundaryLayerWorkflow.state.station2, FlowRegimeEnum::Turbulent);

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
  boundaryLayerWorkflow.state.station2 = blvar(boundaryLayerWorkflow.state.station2, FlowRegimeEnum::Turbulent);

  boundaryLayerWorkflow.state.stepbl();
  restoreblData(2);

  //---- calculate xt-x2 midpoint cfm value
  SkinFrictionCoefficients turbulentSkinFriction =
      blmid(boundaryLayerWorkflow.state, FlowRegimeEnum::Turbulent);

  //---- set up newton system for dct, dth, dds, due, dxi  at  xt and x2
  blc = blDiffSolver.solve(FlowRegimeEnum::Turbulent, boundaryLayerWorkflow.state, turbulentSkinFriction, amcrit);

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
  bt1.block(0, 0, 3, 5) = blc.a1.block(0, 0, 3, 5) * bt1_right;
  bt2.block(0, 0, 3, 5) = blc.a1.block(0, 0, 3, 5) * bt2_right;
  bt2 += blc.a2;
  for (int k = 0; k < 3; k++) {
    btrez[k] = blc.rhs[k];
    btm[k] = blc.d_msq[k] + blc.a1(k, 0) * st_ms + blc.a1(k, 1) * tt_ms +
             blc.a1(k, 2) * dt_ms + blc.a1(k, 3) * ut_ms + blc.a1(k, 4) * xt_ms;
    btr[k] = blc.d_re[k] + blc.a1(k, 0) * st_re + blc.a1(k, 1) * tt_re +
             blc.a1(k, 2) * dt_re + blc.a1(k, 3) * ut_re + blc.a1(k, 4) * xt_re;
    btx[k] = blc.d_xi[k] + blc.a1(k, 0) * st_xf + blc.a1(k, 1) * tt_xf +
             blc.a1(k, 2) * dt_xf + blc.a1(k, 3) * ut_xf + blc.a1(k, 4) * xt_xf;
  }

  //---- add up laminar and turbulent parts to get final system
  //-    in terms of honest-to-god "1" and "2" variables.
  blc.rhs[0] = btrez[0];
  blc.rhs[1] = blrez[1] + btrez[1];
  blc.rhs[2] = blrez[2] + btrez[2];
  blc.d_msq[0] = btm[0];
  blc.d_msq[1] = blm[1] + btm[1];
  blc.d_msq[2] = blm[2] + btm[2];
  blc.d_re[0] = btr[0];
  blc.d_re[1] = blr[1] + btr[1];
  blc.d_re[2] = blr[2] + btr[2];
  blc.d_xi[0] = btx[0];
  blc.d_xi[1] = blx[1] + btx[1];
  blc.d_xi[2] = blx[2] + btx[2];
  blc.a1.row(0) = bt1.row(0);
  blc.a2.row(0) = bt2.row(0);
  blc.a1.middleRows(1, 2) = bl1.middleRows(1, 2) + bt1.middleRows(1, 2);
  blc.a2.middleRows(1, 2) = bl2.middleRows(1, 2) + bt2.middleRows(1, 2);

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
    for (int ibl = 0; ibl < boundaryLayerWorkflow.lattice.stationCount.get(is) - 1; ++ibl) {
      double dui = 0.0;
      for (int js = 1; js <= 2; js++) {
        for (int jbl = 0; jbl < boundaryLayerWorkflow.lattice.stationCount.get(js) - 1; ++jbl) {
          double ue_m = -boundaryLayerWorkflow.lattice.vti.get(is)[ibl] * boundaryLayerWorkflow.lattice.vti.get(js)[jbl] *
                        dij(boundaryLayerWorkflow.lattice.stationToPanel.get(is)[ibl],
                            boundaryLayerWorkflow.lattice.stationToPanel.get(js)[jbl]);
          dui += ue_m * boundaryLayerWorkflow.lattice.mass.get(js)[jbl];
        }
      }
      boundaryLayerWorkflow.lattice.uedg.get(is)[ibl] = boundaryLayerWorkflow.lattice.uinv.get(is)[ibl] + dui;
    }
  }
  return true;
}


bool XFoil::uicalc() {
  //--------------------------------------------------------------
  //     sets inviscid ue from panel inviscid tangential velocity
  //--------------------------------------------------------------
  for (int is = 1; is <= 2; is++) {
    boundaryLayerWorkflow.lattice.uinv.get(is)[0] = 0.0;
    boundaryLayerWorkflow.lattice.uinv_a.get(is)[0] = 0.0;
    for (int ibl = 0; ibl < boundaryLayerWorkflow.lattice.stationCount.get(is) - 1; ++ibl) {
      int i = boundaryLayerWorkflow.lattice.stationToPanel.get(is)[ibl];
      boundaryLayerWorkflow.lattice.uinv.get(is)[ibl] = boundaryLayerWorkflow.lattice.vti.get(is)[ibl] * qinv[i];
      boundaryLayerWorkflow.lattice.uinv_a.get(is)[ibl] = boundaryLayerWorkflow.lattice.vti.get(is)[ibl] * qinv_a[i];
    }
  }

  return true;
}


bool XFoil::xicalc() {
  //-------------------------------------------------------------
  //     sets bl arc length array on each airfoil side and wake
  //-------------------------------------------------------------

  
    for (int ibl = 0; ibl <= boundaryLayerWorkflow.lattice.trailingEdgeIndex.top; ++ibl) {
      boundaryLayerWorkflow.lattice.xssi.top[ibl] = sst - foil.foil_shape.spline_length[boundaryLayerWorkflow.lattice.stationToPanel.get(1)[ibl]];
    }
  
    for (int ibl = 0; ibl <= boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom; ++ibl) {
      boundaryLayerWorkflow.lattice.xssi.bottom[ibl] = foil.foil_shape.spline_length[boundaryLayerWorkflow.lattice.stationToPanel.get(2)[ibl]] - sst;
    }

    // Wake: start from TE, duplicate TE value at first wake station
    boundaryLayerWorkflow.lattice.xssi.bottom[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom + 1] = boundaryLayerWorkflow.lattice.xssi.bottom[boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom];
    for (int ibl = boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom + 2; ibl < boundaryLayerWorkflow.lattice.stationCount.bottom; ++ibl) {
      boundaryLayerWorkflow.lattice.xssi.bottom[ibl] = boundaryLayerWorkflow.lattice.xssi.bottom[ibl - 1] +
                          (foil.wake_shape.points.col(boundaryLayerWorkflow.lattice.stationToPanel.get(2)[ibl]) -
                           foil.wake_shape.points.col(boundaryLayerWorkflow.lattice.stationToPanel.get(2)[ibl - 1]))
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
    for (int iw0 = 0; iw0 < nw; iw0++)
      wgap[iw0] = 0.0;
  }

  else {
    //----- set te flap (wake gap) array (0-based: iw0=0..nw-1)
    for (int iw0 = 0; iw0 < nw; iw0++) {
      const int te_bot_0b = boundaryLayerWorkflow.lattice.trailingEdgeIndex.bottom; // 0-based TE for array indexing
      const double zn = 1.0 - (boundaryLayerWorkflow.lattice.xssi.bottom[te_bot_0b + (iw0 + 1)] -
                               boundaryLayerWorkflow.lattice.xssi.bottom[te_bot_0b]) /
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

  if (boundaryLayerWorkflow.lattice.transitionLocation.get(is) >= 1.0) {
    return boundaryLayerWorkflow.lattice.xssi.get(is)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is)];
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
    str = foil.edge.sle + (foil.foil_shape.spline_length[0] - foil.edge.sle) * boundaryLayerWorkflow.lattice.transitionLocation.top;
  } else {
    str = foil.edge.sle + (foil.foil_shape.spline_length[foil.foil_shape.n - 1] - foil.edge.sle) * boundaryLayerWorkflow.lattice.transitionLocation.bottom;
  }
  str = spline::sinvrt(str, boundaryLayerWorkflow.lattice.transitionLocation.get(is), w1, w3, foil.foil_shape.spline_length.head(point_count), point_count);
  xiforc = std::min((str - sst), boundaryLayerWorkflow.lattice.xssi.get(is)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is)]);
  if (xiforc < 0.0) {
    ss << " ***  stagnation point is past trip on side " << is << "\n";
    writeString(ss.str());

    xiforc = boundaryLayerWorkflow.lattice.xssi.get(is)[boundaryLayerWorkflow.lattice.trailingEdgeIndex.get(is)];
  }

  return xiforc;
}
