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


/** -----------------------------------------------------------
 *     sets  bl location -> panel location  pointer array ipan
 * -----------------------------------------------------------*/
bool XFoil::iblpan() {
  std::stringstream ss;
  const int point_count = foil.foil_shape.n;

  //-- top surface first
  // store ipan with 0-based BL station index, and set vti at 0-based
  for (int i = 0; i <= i_stagnation; i++) {
    ipan.top[i] = i_stagnation - i; // panel index
    vti.top[i] = 1.0;
  }
  // store TE as 0-based logical index
  iblte.top = i_stagnation;
  // nbl exclusive upper bound (0-based TE -> 1-based TE station + 1)
  nbl.top = iblte.top + 2;

  //-- bottom surface next
  // Bottom side: station 0 just after stagnation on bottom
  for (int index = 0; index <= point_count - i_stagnation; ++index) {
    ipan.bottom[index] = i_stagnation + 1 + index;
    vti.bottom[index] = -1.0;
  }

  //-- wake
  iblte.bottom = point_count - i_stagnation - 2; // logical 0-based TE

  for (int iw = 0; iw < nw; iw++) {
    int i = point_count + iw; // panel index in wake
    int index = iblte.bottom + iw + 2; // 1-based BL station for wake (bottom)
    ipan.bottom[index - 1] = i;        // ipan is 0-based in BL station
    vti.bottom[index - 1] = -1.0;
  }

  nbl.bottom = iblte.bottom + nw + 2;

  //-- upper wake pointers (for plotting only)
  for (int iw = 0; iw < nw; iw++) {
    // copy wake panel pointer from bottom to top (for plotting)
    ipan.top[iblte.top + iw + 1] = ipan.bottom[iblte.bottom + iw + 1]; // both sides are 0-based indices, hence -1 vs vti
    vti.top[iblte.top + iw + 1] = 1.0;
  }
  int iblmax = std::max(iblte.top, iblte.bottom) + nw + 2;
  if (iblmax > IVX) {
    ss << "iblpan :  ***  bl array overflow\n";
    ss << "Increase IVX to at least " << iblmax << "\n";
    writeString(ss.str());
    return false;
  }

  lipan = true;
  return true;
}


/** ---------------------------------------------
 *     sets the bl newton system line number
 *     corresponding to each bl station.
 * --------------------------------------------- */
bool XFoil::iblsys() {
  int iv = 0;
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      iv++;
      isys.get(is)[ibl] = iv;
    }
  }

  nsys = iv;
  if (nsys > 2 * IVX) {
    writeString("*** iblsys: bl system array overflow. ***");
    return false;
  }

  return true;
}

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
    int itrold = itran.get(is);

    tran = false;
    turb = false;
    itran.get(is) = iblte.get(is);

    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      MixedModeStationContext ctx = prepareMixedModeStation(is, ibl, itrold, ami);
      bool converged = performMixedModeNewtonIteration(is, ibl, itrold, ctx, deps,
                                                       senswt, sens, sennew, ami);
      if (!converged) {
        handleMixedModeNonConvergence(is, ibl, ctx, ami);
      }

      sens = sennew;

      if (ibl < itran.get(is)) {
        ctau.get(is)[ibl] = ctx.ami;
      } else {
        ctau.get(is)[ibl] = ctx.cti;
      }
      thet.get(is)[ibl] = ctx.thi;
      dstr.get(is)[ibl] = ctx.dsi;
      uedg.get(is)[ibl] = ctx.uei;
      mass.get(is)[ibl] = ctx.dsi * ctx.uei;
      ctq.get(is)[ibl] = blData2.cqz.scalar;

      {
        blData updatedCurrent =
            blprv(state.current(), ctx.xsi, ctx.ami, ctx.cti, ctx.thi,
                  ctx.dsi, ctx.dswaki, ctx.uei);
        state.current() = updatedCurrent;
      }
      blkin(state);

      stepbl(state);

      if (tran || ibl == iblte.get(is)) {
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

bool XFoil::mrchdu() { return mrchdu(boundaryLayerState, boundaryLayerLattice); }

XFoil::MixedModeStationContext XFoil::prepareMixedModeStation(int side, int ibl,
                                                              int itrold,
                                                              double& ami) {
  MixedModeStationContext ctx;

  ctx.simi = (ibl == 0);
  ctx.wake = ibl > iblte.get(side);
  ctx.xsi = xssi.get(side)[ibl];
  ctx.uei = uedg.get(side)[ibl];
  ctx.thi = thet.get(side)[ibl];
  ctx.dsi = dstr.get(side)[ibl];

  if (ibl < itrold) {
    ami = ctau.get(side)[ibl];
    ctx.cti = 0.03;
  } else {
    ctx.cti = ctau.get(side)[ibl];
    if (ctx.cti <= 0.0) {
      ctx.cti = 0.03;
    }
  }
  ctx.ami = ami;

  if (ctx.wake) {
    int iw = ibl - iblte.get(side);
    ctx.dswaki = wgap[iw - 1];
  } else {
    ctx.dswaki = 0.0;
  }

  double thickness_limit = (ibl <= iblte.get(side)) ? 1.02 : 1.00005;
  ctx.dsi = std::max(ctx.dsi - ctx.dswaki, thickness_limit * ctx.thi) + ctx.dswaki;

  simi = ctx.simi;
  wake = ctx.wake;

  return ctx;
}

bool XFoil::isStartOfWake(int side, int stationIndex) const {
  return stationIndex == iblte.get(side) + 1;
}

void XFoil::updateSystemMatricesForStation(int side, int stationIndex,
                                           MixedModeStationContext& ctx) {
  if (isStartOfWake(side, stationIndex)) {
    ctx.tte = thet.get(1)[iblte.top] + thet.get(2)[iblte.bottom];
    ctx.dte =
        dstr.get(1)[iblte.top] + dstr.get(2)[iblte.bottom] + foil.edge.ante;
    ctx.cte =
        (ctau.get(1)[iblte.top] * thet.get(1)[iblte.top] +
         ctau.get(2)[iblte.bottom] * thet.get(2)[iblte.bottom]) /
        ctx.tte;
    tesys(ctx.cte, ctx.tte, ctx.dte);
  } else {
    blsys(boundaryLayerState, boundaryLayerLattice);
  }
}

void XFoil::initializeFirstIterationState(int side, int stationIndex,
                                          int previousTransition,
                                          MixedModeStationContext& ctx,
                                          double& ueref, double& hkref,
                                          double& ami) {
  ueref = blData2.param.uz;
  hkref = blData2.hkz.scalar;

  const bool inLaminarWindow =
      stationIndex < itran.get(side) && stationIndex >= previousTransition;
  if (inLaminarWindow) {
    double uem;
    double dsm;
    double thm;
    if (stationIndex > 0) {
      uem = uedg.get(side)[stationIndex - 1];
      dsm = dstr.get(side)[stationIndex - 1];
      thm = thet.get(side)[stationIndex - 1];
    } else {
      uem = uedg.get(side)[stationIndex];
      dsm = dstr.get(side)[stationIndex];
      thm = thet.get(side)[stationIndex];
    }
    const double uem_sq = uem * uem;
    const double msq =
        uem_sq * hstinv / (gm1bl * (1.0 - 0.5 * uem_sq * hstinv));
    const auto hkin_result = boundary_layer::hkin(dsm / thm, msq);
    hkref = hkin_result.hk;
  }

  if (stationIndex < previousTransition) {
    if (tran) {
      ctau.get(side)[stationIndex] = 0.03;
    }
    if (turb) {
      const double prev =
          (stationIndex >= 1) ? ctau.get(side)[stationIndex - 1]
                              : ctau.get(side)[stationIndex];
      ctau.get(side)[stationIndex] = prev;
    }
    if (tran || turb) {
      ctx.cti = ctau.get(side)[stationIndex - 1];
      blData2.param.sz = ctx.cti;
    }
  }
}

void XFoil::configureSimilarityRow(double ueref) {
  blc.a2(3, 0) = 0.0;
  blc.a2(3, 1) = 0.0;
  blc.a2(3, 2) = 0.0;
  blc.a2(3, 3) = blData2.param.uz_uei;
  blc.rhs[3] = ueref - blData2.param.uz;
}

void XFoil::configureViscousRow(double hkref, double ueref, double senswt,
                                bool resetSensitivity, bool averageSensitivity,
                                double& sens, double& sennew) {
  blc.a2(3, 0) = 0.0;
  blc.a2(3, 1) = blData2.hkz.t();
  blc.a2(3, 2) = blData2.hkz.d();
  blc.a2(3, 3) = blData2.hkz.u() * blData2.param.uz_uei;
  blc.rhs[3] = 1.0;

  const double delta_sen =
      blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs)[3];

  sennew = senswt * delta_sen * hkref / ueref;
  if (resetSensitivity) {
    sens = sennew;
  } else if (averageSensitivity) {
    sens = 0.5 * (sens + sennew);
  }

  blc.a2(3, 1) = blData2.hkz.t() * hkref;
  blc.a2(3, 2) = blData2.hkz.d() * hkref;
  blc.a2(3, 3) =
      (blData2.hkz.u() * hkref + sens / ueref) * blData2.param.uz_uei;
  blc.rhs[3] = -(hkref * hkref) * (blData2.hkz.scalar / hkref - 1.0) -
               sens * (blData2.param.uz / ueref - 1.0);
}

bool XFoil::applyMixedModeNewtonStep(int side, int stationIndex, double deps,
                                     double& ami,
                                     MixedModeStationContext& ctx) {
  blc.rhs = blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);

  ctx.dmax =
      std::max(std::fabs(blc.rhs[1] / ctx.thi), std::fabs(blc.rhs[2] / ctx.dsi));
  if (stationIndex >= itran.get(side)) {
    ctx.dmax = std::max(ctx.dmax, std::fabs(blc.rhs[0] / (10.0 * ctx.cti)));
  }

  rlx = 1.0;
  if (ctx.dmax > 0.3) {
    rlx = 0.3 / ctx.dmax;
  }

  if (stationIndex < itran.get(side)) {
    ami += rlx * blc.rhs[0];
    ctx.ami = ami;
  }
  if (stationIndex >= itran.get(side)) {
    ctx.cti += rlx * blc.rhs[0];
  }
  ctx.thi += rlx * blc.rhs[1];
  ctx.dsi += rlx * blc.rhs[2];
  ctx.uei += rlx * blc.rhs[3];

  if (stationIndex >= itran.get(side)) {
    ctx.cti = std::clamp(ctx.cti, 0.0000001, 0.30);
  }

  const double hklim = (stationIndex <= iblte.get(side)) ? 1.02 : 1.00005;
  const double uei_sq = ctx.uei * ctx.uei;
  const double msq = uei_sq * hstinv /
                     (gm1bl * (1.0 - 0.5 * uei_sq * hstinv));
  double dsw = ctx.dsi - ctx.dswaki;
  dslim(dsw, ctx.thi, msq, hklim);
  ctx.dsi = dsw + ctx.dswaki;

  return ctx.dmax <= deps;
}

void XFoil::checkTransitionIfNeeded(int side, int ibl, bool skipCheck,
                                    int laminarAdvance, double& ami) {
  if (skipCheck || turb) {
    return;
  }

  trchek();
  ami = blData2.param.amplz;
  if (tran) {
    itran.get(side) = ibl;
  } else {
    itran.get(side) = ibl + laminarAdvance;
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
      blData updatedCurrent = blprv(boundaryLayerState.current(), ctx.xsi, ami,
                                    ctx.cti, ctx.thi, ctx.dsi, ctx.dswaki,
                                    ctx.uei);
      boundaryLayerState.current() = updatedCurrent;
    }
    blkin(boundaryLayerState);

    checkTransitionIfNeeded(side, ibl, ctx.simi, 1, ami);

    const bool startOfWake = isStartOfWake(side, ibl);
    updateSystemMatricesForStation(side, ibl, ctx);

    if (itbl == 1) {
      initializeFirstIterationState(side, ibl, itrold, ctx, ueref, hkref, ami);
    }

    if (ctx.simi || startOfWake) {
      configureSimilarityRow(ueref);
    } else {
      const bool resetSensitivity = (itbl <= 5);
      const bool averageSensitivity = (itbl > 5 && itbl <= 15);
      configureViscousRow(hkref, ueref, senswt, resetSensitivity,
                          averageSensitivity, sens, sennew);
    }

    if (applyMixedModeNewtonStep(side, ibl, deps, ami, ctx)) {
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
    if (ibl <= iblte.get(side)) {
      ctx.thi = thet.get(side)[ibl - 1] *
                sqrt(xssi.get(side)[ibl] / xssi.get(side)[ibl - 1]);
      ctx.dsi = dstr.get(side)[ibl - 1] *
                sqrt(xssi.get(side)[ibl] / xssi.get(side)[ibl - 1]);
      ctx.uei = uedg.get(side)[ibl - 1];
    } else {
      if (ibl == iblte.get(side) + 1) {
        ctx.cti = ctx.cte;
        ctx.thi = ctx.tte;
        ctx.dsi = ctx.dte;
        ctx.uei = uedg.get(side)[ibl - 1];
      } else {
        ctx.thi = thet.get(side)[ibl - 1];
        double ratlen = (xssi.get(side)[ibl] - xssi.get(side)[ibl - 1]) /
                        (10.0 * dstr.get(side)[ibl - 1]);
        ctx.dsi = (dstr.get(side)[ibl - 1] + ctx.thi * ratlen) /
                  (1.0 + ratlen);
        ctx.uei = uedg.get(side)[ibl - 1];
      }
    }

    if (ibl == itran.get(side)) {
      ctx.cti = 0.05;
    }
    if (ibl > itran.get(side)) {
      ctx.cti = ctau.get(side)[ibl - 1];
    }
  }

  {
    blData updatedCurrent = blprv(boundaryLayerState.current(), ctx.xsi, ami,
                                  ctx.cti, ctx.thi, ctx.dsi, ctx.dswaki,
                                  ctx.uei);
    boundaryLayerState.current() = updatedCurrent;
  }
  blkin(boundaryLayerState);

  checkTransitionIfNeeded(side, ibl, ctx.simi, 2, ami);

  if (ibl < itran.get(side)) {
    blData2 = blvar(blData2, FlowRegimeEnum::Laminar);
    blmid(boundaryLayerState, FlowRegimeEnum::Laminar);
  }
  if (ibl >= itran.get(side)) {
    blData2 = blvar(blData2, FlowRegimeEnum::Turbulent);
    blmid(boundaryLayerState, FlowRegimeEnum::Turbulent);
  }
  if (ctx.wake) {
    blData2 = blvar(blData2, FlowRegimeEnum::Wake);
    blmid(boundaryLayerState, FlowRegimeEnum::Wake);
  }

  ctx.ami = ami;
}

double XFoil::calcHtarg(int ibl, int is, bool wake) {
  if (ibl < itran.get(is)) {
    return blData1.hkz.scalar +
           0.03 * (blData2.param.xz - blData1.param.xz) / blData1.param.tz;
  } else if (ibl == itran.get(is)) {
    return blData1.hkz.scalar +
           (0.03 * (xt - blData1.param.xz) - 0.15 * (blData2.param.xz - xt)) /
               blData1.param.tz;
  } else if (wake) {
    const double cst =
        0.03 * (blData2.param.xz - blData1.param.xz) / blData1.param.tz;
    auto euler = [](double hk2, double hk1, double cst) {
      return hk2 - (hk2 + cst * pow(hk2 - 1, 3) - hk1) /
                       (1 + 3 * cst * pow(hk2 - 1, 2));
    };
    blData2.hkz.scalar = blData1.hkz.scalar;
    for (int i = 0; i < 3; i++) {
      blData2.hkz.scalar = euler(blData2.hkz.scalar, blData1.hkz.scalar, cst);
    }
    return blData2.hkz.scalar;
  } else {
    return blData1.hkz.scalar -
           0.15 * (blData2.param.xz - blData1.param.xz) / blData1.param.tz;
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
    xsi = xssi.get(is)[0];
    uei = uedg.get(is)[0];

    ucon = uei / xsi;
    tsq = 0.45 / (ucon * 6.0 * reybl);
    thi = sqrt(tsq);
    dsi = 2.2 * thi;
    ami = 0.0;

    //---- initialize ctau for first turbulent station
    cti = 0.03;

    tran = false;
    turb = false;
    itran.get(is) = iblte.get(is);

    //---- march downstream
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      //int ibl = ibl - 1;
      int ibm = ibl + 1;
      int iw = ibl - iblte.get(is);
      simi = (ibl == 0);
      wake = ibl > iblte.get(is);

      //------ prescribed quantities (xssi now 0-based)
      xsi = xssi.get(is)[ibl];
      uei = uedg.get(is)[ibl];

      if (wake) {
        iw = ibl - iblte.get(is);
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
              blprv(boundaryLayerState.current(), xsi, ami, cti, thi, dsi,
                    dswaki, uei);
          boundaryLayerState.current() = updatedCurrent;
        }
        blkin(boundaryLayerState);

        //-------- check for transition and set appropriate flags and things
        if ((!simi) && (!turb)) {
          trchek();
          ami = blData2.param.amplz;

          //--------- fixed bug   md 7 jun 99
          if (tran) {
            itran.get(is) = ibl;
            if (cti <= 0.0) {
              cti = 0.03;
              blData2.param.sz = cti;
            }
          } else
            itran.get(is) = ibl + 2;
        }

        if (ibl == iblte.get(is) + 1) {
          tte = thet.get(1)[iblte.top] +
                thet.get(2)[iblte.bottom];
          dte = dstr.get(1)[iblte.top] +
                dstr.get(2)[iblte.bottom] + foil.edge.ante;
          cte = (ctau.get(1)[iblte.top] *
                     thet.get(1)[iblte.top] +
                 ctau.get(2)[iblte.bottom] *
                     thet.get(2)[iblte.bottom]) /
                tte;
          tesys(cte, tte, dte);
        } else
          blsys(boundaryLayerState, boundaryLayerLattice);

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
          if (ibl < itran.get(is))
            dmax = std::max(dmax, fabs(blc.rhs[0] / 10.0));
          if (ibl >= itran.get(is))
            dmax = std::max(dmax, fabs(blc.rhs[0] / cti));

          rlx = 1.0;
          if (dmax > 0.3)
            rlx = 0.3 / dmax;
          //--------- see if direct mode is not applicable
          if (ibl != iblte.get(is) + 1) {
            //---------- calculate resulting kinematic shape parameter hk
            msq =
                uei * uei * hstinv / (gm1bl * (1.0 - 0.5 * uei * uei * hstinv));
            htest = (dsi + rlx * blc.rhs[2]) / (thi + rlx * blc.rhs[1]);
            boundary_layer::KineticShapeParameterResult hkin_result =
                boundary_layer::hkin(htest, msq);
            hktest = hkin_result.hk;

            //---------- decide whether to do direct or inverse problem based on
            // hk
            if (ibl < itran.get(is))
              hmax = hlmax;
            if (ibl >= itran.get(is))
              hmax = htmax;
            direct = (hktest < hmax);
          }
          if (direct) {
            //---------- update as usual
            if (ibl >= itran.get(is))
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
          blc.a2(3, 1) = blData2.hkz.t();
          blc.a2(3, 2) = blData2.hkz.d();
          blc.a2(3, 3) = blData2.hkz.u();
          blc.rhs[3] = htarg - blData2.hkz.scalar;
          blc.rhs = blc.a2.block(0, 0, 4, 4).fullPivLu().solve(blc.rhs);

          dmax = std::max(fabs(blc.rhs[1] / thi), fabs(blc.rhs[2] / dsi));
          if (ibl >= itran.get(is))
            dmax = std::max(dmax, fabs(blc.rhs[0] / cti));
          rlx = 1.0;
          if (dmax > 0.3)
            rlx = 0.3 / dmax;
          //--------- update variables
          if (ibl >= itran.get(is))
            cti = cti + rlx * blc.rhs[0];
          thi = thi + rlx * blc.rhs[1];
          dsi = dsi + rlx * blc.rhs[2];
          uei = uei + rlx * blc.rhs[3];
        }
        //-------- eliminate absurd transients

        if (ibl >= itran.get(is)) {
          cti = std::min(cti, 0.30);
          cti = std::max(cti, 0.0000001);
        }
        if (ibl <= iblte.get(is))
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
            if (ibl <= iblte.get(is)) {
              thi = thet.get(is)[ibl - 1] *
                    sqrt(xssi.get(is)[ibl] / xssi.get(is)[ibl - 1]);
              dsi = dstr.get(is)[ibl - 1] *
                    sqrt(xssi.get(is)[ibl] / xssi.get(is)[ibl - 1]);
            } else {
              if (ibl == iblte.get(is) + 1) {
                cti = cte;
                thi = tte;
                dsi = dte;
              } else {
                thi = thet.get(is)[ibl - 1];
                ratlen = (xssi.get(is)[ibl] - xssi.get(is)[ibl - 1]) /
                         (10.0 * dstr.get(is)[ibl - 1]);
                dsi = (dstr.get(is)[ibl - 1] +
                       thi * ratlen) /
                      (1.0 + ratlen);
              }
            }
            if (ibl == itran.get(is))
              cti = 0.05;
            if (ibl > itran.get(is))
              cti = ctau.get(is)[ibl - 1];

            uei = uedg.get(is)[ibl];

            if (ibl < nbl.get(is) - 1)
              uei = 0.5 *
                    (uedg.get(is)[ibl - 1] +
                     uedg.get(is)[ibl + 1]);
          }
        }
        // 109
        {
          blData updatedCurrent =
              blprv(boundaryLayerState.current(), xsi, ami, cti, thi, dsi,
                    dswaki, uei);
          boundaryLayerState.current() = updatedCurrent;
        }
        blkin(boundaryLayerState);
        //------- check for transition and set appropriate flags and things
        if ((!simi) && (!turb)) {
          trchek();
          ami = blData2.param.amplz;
          if (tran)
            itran.get(is) = ibl;
          if (!tran)
            itran.get(is) = ibl + 2;
        }
        //------- set all other extrapolated values for current station
        if (ibl < itran.get(is))
          blData2 = blvar(blData2, FlowRegimeEnum::Laminar);
        if (ibl >= itran.get(is))
          blData2 = blvar(blData2, FlowRegimeEnum::Turbulent);
        if (wake)
          blData2 = blvar(blData2, FlowRegimeEnum::Wake);
        if (ibl < itran.get(is))
          blmid(boundaryLayerState, FlowRegimeEnum::Laminar);
        if (ibl >= itran.get(is))
          blmid(boundaryLayerState, FlowRegimeEnum::Turbulent);
        if (wake)
          blmid(boundaryLayerState, FlowRegimeEnum::Wake);
      }
      //------ store primary variables
      if (ibl < itran.get(is))
        ctau.get(is)[ibl] = ami;
      if (ibl >= itran.get(is))
        ctau.get(is)[ibl] = cti;
      thet.get(is)[ibl] = thi;
      dstr.get(is)[ibl] = dsi;
      uedg.get(is)[ibl] = uei;
      mass.get(is)[ibl] = dsi * uei;
      ctq.get(is)[ibl] = blData2.cqz.scalar;

      //------ set "1" variables to "2" variables for next streamwise station
      {
        blData updatedCurrent =
            blprv(boundaryLayerState.current(), xsi, ami, cti, thi, dsi,
                  dswaki, uei);
        boundaryLayerState.current() = updatedCurrent;
      }
      blkin(boundaryLayerState);

      stepbl(boundaryLayerState);

      //------ turbulent intervals will follow transition interval or te
      if (tran || ibl == iblte.get(is)) {
        turb = true;
      }

      tran = false;

      if (ibl == iblte.get(is)) {
        thi = thet.get(1)[iblte.top] +
              thet.get(2)[iblte.bottom];
        dsi = dstr.get(1)[iblte.top] +
              dstr.get(2)[iblte.bottom] + foil.edge.ante;
      }
    }
  }
  return true;
}

SetblInputView XFoil::makeSetblInputView() const {
  return SetblInputView{lblini,
                        uedg,
                        ctau,
                        thet,
                        dstr,
                        mass,
                        ctq,
                        itran};
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
                         uedg,
                         ctau,
                         thet,
                         dstr,
                         mass,
                         ctq,
                         itran,
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
  jvte1 = isys.top[iblte.top];
  jvte2 = isys.bottom[iblte.bottom];

  dule1 = output.uedg.top[0] - usav.top[0];
  dule2 = output.uedg.bottom[0] - usav.bottom[0];

  //---- set le and te ue sensitivities wrt all m values
  const auto le_te_sensitivities = computeLeTeSensitivities(
      ipan.get(1)[0], ipan.get(2)[0], ipan.get(1)[iblte.top],
      ipan.get(2)[iblte.bottom]);
  ule1_m = le_te_sensitivities.ule1_m;
  ule2_m = le_te_sensitivities.ule2_m;
  ute1_m = le_te_sensitivities.ute1_m;
  ute2_m = le_te_sensitivities.ute2_m;

  ule1_a = uinv_a.get(1)[0];
  ule2_a = uinv_a.get(2)[0];

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
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      
      int iv = isys.get(is)[ibl];

      output.simi = (ibl == 0);
      output.wake = (ibl > iblte.get(is));
      output.tran = (ibl == output.itran.get(is));
      output.turb = (ibl > output.itran.get(is));

      //---- set primary variables for current station
      xsi = xssi.get(is)[ibl];
      if (ibl < output.itran.get(is))
        ami = output.ctau.get(is)[ibl];
      else
        cti = output.ctau.get(is)[ibl];
      uei = output.uedg.get(is)[ibl];
      thi = output.thet.get(is)[ibl];
      mdi = output.mass.get(is)[ibl];

      dsi = mdi / uei;

      if (output.wake) {
        int iw = ibl - iblte.get(is);
        dswaki = wgap[iw - 1];
      } else
        dswaki = 0.0;

      //---- set derivatives of dsi (= d2)
      d2_m2 = 1.0 / uei;
      d2_u2 = -dsi / uei;

      for (int js = 1; js <= 2; js++) {
        for (int jbl = 0; jbl < nbl.get(js) - 1; ++jbl) {
          int jv = isys.get(js)[jbl];
          u2_m[jv] = -vti.get(is)[ibl] * vti.get(js)[jbl] *
                     dij(ipan.get(is)[ibl], ipan.get(js)[jbl]);
          d2_m[jv] = d2_u2 * u2_m[jv];
        }
      }
      d2_m[iv] = d2_m[iv] + d2_m2;

      u2_a = uinv_a.get(is)[ibl];
      d2_a = d2_u2 * u2_a;

  //---- "forced" changes due to mismatch between edge velocities and
  // usav=uinv+dij*output.mass
      due2 = output.uedg.get(is)[ibl] - usav.get(is)[ibl];
      dds2 = d2_u2 * due2;

      {
        blData updatedCurrent =
            blprv(boundaryLayerState.current(), xsi, ami, cti, thi, dsi,
                  dswaki, uei);
        boundaryLayerState.current() = updatedCurrent;
      } // cti
      blkin(boundaryLayerState);

      //---- check for transition and set output.tran, xt, etc. if found
      if (output.tran) {
        trchek();
        ami = blData2.param.amplz;
      }

      if (ibl == output.itran.get(is) && !output.tran) {
        // TRACE("setbl: xtr???  n1=%d n2=%d: \n", ampl1, ampl2);

        ss << "setbl: xtr???  n1=" << blData1.param.amplz
           << " n2=" << blData2.param.amplz << ":\n";
        writeString(ss.str());
        ss.str("");
      }

      //---- assemble 10x4 linearized system for dctau, dth, dds, due, dxi
      //	   at the previous "1" station and the current "2" station

      if (ibl == iblte.get(is) + 1) {
        //----- define quantities at start of output.wake, adding te base thickness to
        // dstar
        tte = output.thet.get(1)[iblte.top] +
              output.thet.get(2)[iblte.bottom];
        dte = output.dstr.get(1)[iblte.top] +
              output.dstr.get(2)[iblte.bottom] + foil.edge.ante;
        cte = (output.ctau.get(1)[iblte.top] *
                   output.thet.get(1)[iblte.top] +
               output.ctau.get(2)[iblte.bottom] *
                   output.thet.get(2)[iblte.bottom]) /
               tte;
        tesys(cte, tte, dte);

        tte_tte1 = 1.0;
        tte_tte2 = 1.0;
        dte_mte1 = 1.0 / output.uedg.top[iblte.top];
        dte_ute1 = -output.dstr.get(1)[iblte.top] /
                    output.uedg.top[iblte.top];
        dte_mte2 = 1.0 / output.uedg.bottom[iblte.bottom];
        dte_ute2 = -output.dstr.get(2)[iblte.bottom] /
                    output.uedg.bottom[iblte.bottom];
        cte_cte1 = output.thet.get(1)[iblte.top] / tte;
        cte_cte2 = output.thet.get(2)[iblte.bottom] / tte;
        cte_tte1 = (output.ctau.get(1)[iblte.top] - cte) / tte;
        cte_tte2 = (output.ctau.get(2)[iblte.bottom] - cte) / tte;

        //----- re-define d1 sensitivities wrt m since d1 depends on both te ds
        // values
      for (int js = 1; js <= 2; js++) {
        for (int jbl = 0; jbl < nbl.get(js) - 1; ++jbl) {
            int jv = isys.get(js)[jbl];
            d1_m[jv] = dte_ute1 * ute1_m[jv] + dte_ute2 * ute2_m[jv];
          }
        }
        d1_m[jvte1] = d1_m[jvte1] + dte_mte1;
        d1_m[jvte2] = d1_m[jvte2] + dte_mte2;

        //----- "forced" changes from  output.uedg --- usav=uinv+dij*output.mass	mismatch
        due1 = 0.0;
        dds1 =
            dte_ute1 *
                (output.uedg.top[iblte.top] -
                 usav.top[iblte.top]) +
            dte_ute2 *
                (output.uedg.bottom[iblte.bottom] -
                 usav.bottom[iblte.bottom]);
      } else {
        blsys(boundaryLayerState, boundaryLayerLattice);
      }

      //---- save wall shear and equil. max shear coefficient for plotting
      // output
      output.ctq.get(is)[ibl] = blData2.cqz.scalar;

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

      if (ibl == iblte.get(is) + 1) {
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

      if (ibl == iblte.get(is)) {
        //----- set "2" variables at te to output.wake correlations for next station

        output.turb = true;
        output.wake = true;
        blData2 = blvar(blData2, FlowRegimeEnum::Wake);
        blmid(boundaryLayerState, FlowRegimeEnum::Wake);
      }
      u1_m = u2_m;
      d1_m = d2_m;

      u1_a = u2_a;
      d1_a = d2_a;

      due1 = due2;
      dds1 = dds2;

      //---- set bl variables for next station
      //			for (icom=1; icom<= ncom;icom++)
      // com1[icom] = com2[icom];
      stepbl(boundaryLayerState);

      //---- next streamwise station
    }

    //---- next airfoil side
  }

  return output;
}


bool XFoil::stepbl(BoundaryLayerState& state) {
  state.previous() = state.current();
  return true;
}


bool XFoil::stfind() {
  //-----------------------------------------
  //     locates stagnation point arc length
  //     location sst and panel index ist.
  //-----------------------------------------

  int i;
  bool bFound = false;
  const int point_count = foil.foil_shape.n;
  for (i = 0; i < point_count - 1; i++) {
    if (surface_vortex(0, i) >= 0.0 && surface_vortex(0, i + 1) < 0.0) {
      bFound = true;
      break;
    }
  }

  if (!bFound) {
    writeString("stfind: Stagnation point not found. Continuing ...\n");
    i = point_count / 2;
  }

  i_stagnation = i;
  const double dgam = surface_vortex(0, i + 1) - surface_vortex(0, i);
  const double ds = foil.foil_shape.spline_length[i + 1] - foil.foil_shape.spline_length[i];

  //---- evaluate so as to minimize roundoff for very small gam[i] or gam[i+1]
  if (surface_vortex(0, i) < -surface_vortex(0, i + 1))
    sst = foil.foil_shape.spline_length[i] - ds * (surface_vortex(0, i) / dgam);
  else
    sst = foil.foil_shape.spline_length[i + 1] - ds * (surface_vortex(0, i + 1) / dgam);

  //---- tweak stagnation point if it falls right on a node (very unlikely)
  if (sst <= foil.foil_shape.spline_length[i])
    sst = foil.foil_shape.spline_length[i] + 0.0000001;
  if (sst >= foil.foil_shape.spline_length[i + 1])
    sst = foil.foil_shape.spline_length[i + 1] - 0.0000001;

  sst_go = (sst - foil.foil_shape.spline_length[i + 1]) / dgam;
  sst_gp = (foil.foil_shape.spline_length[i] - sst) / dgam;

  return true;
}


bool XFoil::stmove() {
  //--------------------------------------------------
  //    moves stagnation point location to new panel.
  //---------------------------------------------------
  
  //-- locate new stagnation point arc length sst from gam distribution
  int istold = i_stagnation;
  stfind();

  if (istold == i_stagnation) {
    //--- recalculate new arc length array
    xicalc();
  } else {
    //--- set new bl position -> panel position  pointers
    iblpan();

    //--- set new inviscid bl edge velocity uinv from qinv
    uicalc();

    //--- recalculate new arc length array
    xicalc();

    //--- set  bl position -> system line  pointers
    iblsys();

    if (i_stagnation > istold) {
      //---- increase in number of points on top side (is=1)
      int idif = i_stagnation - istold;

      itran.top += idif;
      itran.bottom -= idif;

      //---- move top side bl variables downstream
      for (int ibl = nbl.top - 2; ibl >= idif; ibl--) {
        ctau.top[ibl] = ctau.top[ibl - idif];
        thet.top[ibl] = thet.top[ibl - idif];
        dstr.top[ibl] = dstr.top[ibl - idif];
        uedg.top[ibl] = uedg.top[ibl - idif];
      }

      //---- set bl variables between old and new stagnation point
      const double dudx =
          uedg.top[idif] / xssi.top[idif];
      for (int ibl = idif; ibl >= 1; ibl--) {
        ctau.top[ibl - 1] = ctau.top[idif];
        thet.top[ibl - 1] = thet.top[idif];
        dstr.top[ibl - 1] = dstr.top[idif];
        uedg.top[ibl - 1] = dudx * xssi.top[ibl - 1];
      }

      //---- move bottom side bl variables upstream
      for (int ibl = 0; ibl < nbl.bottom - 1; ibl++) {
        ctau.bottom[ibl] = ctau.bottom[ibl + idif];
        thet.bottom[ibl] = thet.bottom[ibl + idif];
        dstr.bottom[ibl] = dstr.bottom[ibl + idif];
        uedg.bottom[ibl] = uedg.bottom[ibl + idif];
      }
    } else {
      //---- increase in number of points on bottom side (is=2)
      int idif = istold - i_stagnation;

      itran.top = itran.top - idif;
      itran.bottom = itran.bottom + idif;

      //---- move bottom side bl variables downstream
      for (int ibl = nbl.bottom - 1; ibl >= idif + 1; ibl--) {
        ctau.bottom[ibl - 1] = ctau.bottom[(ibl - idif) - 1];
        thet.bottom[ibl - 1] = thet.bottom[(ibl - idif) - 1];
        dstr.bottom[ibl - 1] = dstr.bottom[(ibl - idif) - 1];
        uedg.bottom[ibl - 1] = uedg.bottom[(ibl - idif) - 1];
      }

      //---- set bl variables between old and new stagnation point
      const double dudx =
          uedg.bottom[idif] / xssi.bottom[idif];
      for (int ibl = idif; ibl >= 1; ibl--) {
        ctau.bottom[ibl - 1] = ctau.bottom[idif];
        thet.bottom[ibl - 1] = thet.bottom[idif];
        dstr.bottom[ibl - 1] = dstr.bottom[idif];
        uedg.bottom[ibl - 1] = dudx * xssi.bottom[ibl - 1];
      }

      //---- move top side bl variables upstream
      for (int ibl = 0; ibl < nbl.top - 1; ibl++) {
        ctau.top[ibl] = ctau.top[ibl + idif];
        thet.top[ibl] = thet.top[ibl + idif];
        dstr.top[ibl] = dstr.top[ibl + idif];
        uedg.top[ibl] = uedg.top[ibl + idif];
      }
    }
  }

  //-- set new mass array since ue has been tweaked
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      mass.get(is)[ibl] = dstr.get(is)[ibl] * uedg.get(is)[ibl];
    }
  }

  return true;
}


// trailing-edge calculations handled by Edge::updateFromFoilShape()

bool XFoil::tesys(double cte, double tte, double dte) {
  //--------------------------------------------------------
  //	   sets up "dummy" bl system between airfoil te point
  //	   and first wake point infinitesimally behind te.
  //--------------------------------------------------------

  blc.clear();

  blData2 = blvar(blData2, FlowRegimeEnum::Wake);

  blc.a1(0, 0) = -1.0;
  blc.a2(0, 0) = 1.0;
  blc.rhs[0] = cte - blData2.param.sz;

  blc.a1(1, 1) = -1.0;
  blc.a2(1, 1) = 1.0;
  blc.rhs[1] = tte - blData2.param.tz;

  blc.a1(2, 2) = -1.0;
  blc.a2(2, 2) = 1.0;
  blc.rhs[2] = dte - blData2.param.dz - blData2.param.dwz;

  return true;
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
      BoundaryLayerUtil::axset(blData1.hkz.scalar, blData1.param.tz, blData1.rtz.scalar,
            blData1.param.amplz, blData2.hkz.scalar, blData2.param.tz,
            blData2.rtz.scalar, blData2.param.amplz, amcrit);

  //---- set initial guess for iterate n2 (ampl2) at x2
  blData2.param.amplz = blData1.param.amplz +
                        ax_result.ax * (blData2.param.xz - blData1.param.xz);
  //---- solve implicit system for amplification ampl2
  auto iterateAmplification = [&]() -> bool {
    for (int itam = 0; itam < 30; itam++) {
      //---- define weighting factors wf1,wf2 for defining "t" quantities
      if (blData2.param.amplz <= amcrit) {
        //------ there is no transition yet,  "t" is the same as "2"
        amplt = blData2.param.amplz;
        amplt_a2 = 1.0;
        sfa = 1.0;
        sfa_a1 = 0.0;
        sfa_a2 = 0.0;
      } else {
        //------ there is transition in x1..x2, "t" is set from n1, n2
        amplt = amcrit;
        amplt_a2 = 0.0;
        sfa = (amplt - blData1.param.amplz) /
              (blData2.param.amplz - blData1.param.amplz);
        sfa_a1 = (sfa - 1.0) / (blData2.param.amplz - blData1.param.amplz);
        sfa_a2 = (-sfa) / (blData2.param.amplz - blData1.param.amplz);
      }

      if (xiforc < blData2.param.xz) {
        sfx = (xiforc - blData1.param.xz) /
              (blData2.param.xz - blData1.param.xz);
        sfx_x1 = (sfx - 1.0) / (blData2.param.xz - blData1.param.xz);
        sfx_x2 = (-sfx) / (blData2.param.xz - blData1.param.xz);
        sfx_xf = 1.0 / (blData2.param.xz - blData1.param.xz);
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
      xt = blData1.param.xz * (1 - wf) + blData2.param.xz * wf;
      tt = blData1.param.tz * (1 - wf) + blData2.param.tz * wf;
      dt = blData1.param.dz * (1 - wf) + blData2.param.dz * wf;
      ut = blData1.param.uz * (1 - wf) + blData2.param.uz * wf;

      xt_a2 = (blData2.param.xz - blData1.param.xz) * wf_a2;
      tt_a2 = (blData2.param.tz - blData1.param.tz) * wf_a2;
      dt_a2 = (blData2.param.dz - blData1.param.dz) * wf_a2;
      ut_a2 = (blData2.param.uz - blData1.param.uz) * wf_a2;

      //---- temporarily set "2" variables from "t" for blkin
      blData2.param.xz = xt;
      blData2.param.tz = tt;
      blData2.param.dz = dt;
      blData2.param.uz = ut;

      //---- calculate laminar secondary "t" variables hkt, rtt
      blkin(boundaryLayerState);

      blData::blVector hkt = blData2.hkz;
      blData::blVector rtt = blData2.rtz;

      //---- restore clobbered "2" variables, except for ampl2
      amsave = blData2.param.amplz;

      restoreblData(2);

      blData2.param.amplz = amsave;

      //---- calculate amplification rate ax over current x1-xt interval
      ax_result = BoundaryLayerUtil::axset(blData1.hkz.scalar, blData1.param.tz,
                        blData1.rtz.scalar, blData1.param.amplz, hkt.scalar, tt, rtt.scalar,
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
      res = blData2.param.amplz - blData1.param.amplz -
            ax_result.ax * (blData2.param.xz - blData1.param.xz);
      res_a2 = 1.0 - ax_result.ax_a2 * (blData2.param.xz - blData1.param.xz);

      da2 = -res / res_a2;

      rlx = 1.0;
      dxt = xt_a2 * da2;

      if (rlx * fabs(dxt / (blData2.param.xz - blData1.param.xz)) > 0.05) {
        rlx = 0.05 * fabs((blData2.param.xz - blData1.param.xz) / dxt);
      }

      if (rlx * fabs(da2) > 1.0) {
        rlx = 1.0 * fabs(1.0 / da2);
      }

      //---- check if converged
      if (fabs(da2) < daeps) {
        return true;
      }

      if ((blData2.param.amplz > amcrit &&
           blData2.param.amplz + rlx * da2 < amcrit) ||
          (blData2.param.amplz < amcrit &&
           blData2.param.amplz + rlx * da2 > amcrit)) {
        //------ limited newton step so ampl2 doesn't step across amcrit either
        // way
        blData2.param.amplz = amcrit;
      } else {
        //------ regular newton step
        blData2.param.amplz = blData2.param.amplz + rlx * da2;
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
  trfree = (blData2.param.amplz >= amcrit);
  trforc = (xiforc > blData1.param.xz) && (xiforc <= blData2.param.xz);

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

  xt_a1 = (blData2.param.xz - blData1.param.xz) * wf_a1;
  tt_a1 = (blData2.param.tz - blData1.param.tz) * wf_a1;
  dt_a1 = (blData2.param.dz - blData1.param.dz) * wf_a1;
  ut_a1 = (blData2.param.uz - blData1.param.uz) * wf_a1;

  xt_x1 += (blData2.param.xz - blData1.param.xz) * wf_x1;
  tt_x1 = (blData2.param.tz - blData1.param.tz) * wf_x1;
  dt_x1 = (blData2.param.dz - blData1.param.dz) * wf_x1;
  ut_x1 = (blData2.param.uz - blData1.param.uz) * wf_x1;

  xt_x2 += (blData2.param.xz - blData1.param.xz) * wf_x2;
  tt_x2 = (blData2.param.tz - blData1.param.tz) * wf_x2;
  dt_x2 = (blData2.param.dz - blData1.param.dz) * wf_x2;
  ut_x2 = (blData2.param.uz - blData1.param.uz) * wf_x2;

  xt_xf = (blData2.param.xz - blData1.param.xz) * wf_xf;

  //---- at this point, ax = ax( hk1, t1, rt1, a1, hkt, tt, rtt, at )
  blData::blVector hkt = blData2.hkz;
  blData::blVector rtt = blData2.rtz;

  //---- set sensitivities of ax( t1 d1 u1 a1 t2 d2 u2 a2 ms re )
  double ax_t1 = ax_result.ax_hk1 * blData1.hkz.t() + ax_result.ax_t1 +
                 ax_result.ax_rt1 * blData1.rtz.t() +
                 (ax_result.ax_hk2 * hkt.t() + ax_result.ax_t2 +
                  ax_result.ax_rt2 * rtt.t()) *
                     tt_t1;
  double ax_d1 =
      ax_result.ax_hk1 * blData1.hkz.d() + (ax_result.ax_hk2 * hkt.d()) * dt_d1;
  double ax_u1 =
      ax_result.ax_hk1 * blData1.hkz.u() + ax_result.ax_rt1 * blData1.rtz.u() +
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
                 ax_result.ax_hk1 * blData1.hkz.ms() +
                 ax_result.ax_rt1 * blData1.rtz.ms();
  double ax_re =
      ax_result.ax_rt2 * rtt.re() + ax_result.ax_rt1 * blData1.rtz.re();

  //---- set sensitivities of residual res
  z_ax = -(blData2.param.xz - blData1.param.xz);

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
  double wf2 = (xt - blData1.param.xz) / (blData2.param.xz - blData1.param.xz);
  double wf2_xt = 1.0 / (blData2.param.xz - blData1.param.xz);

  double wf2_a1 = wf2_xt * xt_a1;
  double wf2_x1 =
      wf2_xt * xt_x1 + (wf2 - 1.0) / (blData2.param.xz - blData1.param.xz);
  double wf2_x2 = wf2_xt * xt_x2 - wf2 / (blData2.param.xz - blData1.param.xz);
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
  tt = blData1.param.tz * wf1 + blData2.param.tz * wf2;
  tt_a1 = blData1.param.tz * wf1_a1 + blData2.param.tz * wf2_a1;
  tt_x1 = blData1.param.tz * wf1_x1 + blData2.param.tz * wf2_x1;
  tt_x2 = blData1.param.tz * wf1_x2 + blData2.param.tz * wf2_x2;
  tt_t1 = blData1.param.tz * wf1_t1 + blData2.param.tz * wf2_t1 + wf1;
  tt_t2 = blData1.param.tz * wf1_t2 + blData2.param.tz * wf2_t2 + wf2;
  tt_d1 = blData1.param.tz * wf1_d1 + blData2.param.tz * wf2_d1;
  tt_d2 = blData1.param.tz * wf1_d2 + blData2.param.tz * wf2_d2;
  tt_u1 = blData1.param.tz * wf1_u1 + blData2.param.tz * wf2_u1;
  tt_u2 = blData1.param.tz * wf1_u2 + blData2.param.tz * wf2_u2;
  tt_ms = blData1.param.tz * wf1_ms + blData2.param.tz * wf2_ms;
  tt_re = blData1.param.tz * wf1_re + blData2.param.tz * wf2_re;
  tt_xf = blData1.param.tz * wf1_xf + blData2.param.tz * wf2_xf;

  dt = blData1.param.dz * wf1 + blData2.param.dz * wf2;
  dt_a1 = blData1.param.dz * wf1_a1 + blData2.param.dz * wf2_a1;
  dt_x1 = blData1.param.dz * wf1_x1 + blData2.param.dz * wf2_x1;
  dt_x2 = blData1.param.dz * wf1_x2 + blData2.param.dz * wf2_x2;
  dt_t1 = blData1.param.dz * wf1_t1 + blData2.param.dz * wf2_t1;
  dt_t2 = blData1.param.dz * wf1_t2 + blData2.param.dz * wf2_t2;
  dt_d1 = blData1.param.dz * wf1_d1 + blData2.param.dz * wf2_d1 + wf1;
  dt_d2 = blData1.param.dz * wf1_d2 + blData2.param.dz * wf2_d2 + wf2;
  dt_u1 = blData1.param.dz * wf1_u1 + blData2.param.dz * wf2_u1;
  dt_u2 = blData1.param.dz * wf1_u2 + blData2.param.dz * wf2_u2;
  dt_ms = blData1.param.dz * wf1_ms + blData2.param.dz * wf2_ms;
  dt_re = blData1.param.dz * wf1_re + blData2.param.dz * wf2_re;
  dt_xf = blData1.param.dz * wf1_xf + blData2.param.dz * wf2_xf;

  ut = blData1.param.uz * wf1 + blData2.param.uz * wf2;
  ut_a1 = blData1.param.uz * wf1_a1 + blData2.param.uz * wf2_a1;
  ut_x1 = blData1.param.uz * wf1_x1 + blData2.param.uz * wf2_x1;
  ut_x2 = blData1.param.uz * wf1_x2 + blData2.param.uz * wf2_x2;
  ut_t1 = blData1.param.uz * wf1_t1 + blData2.param.uz * wf2_t1;
  ut_t2 = blData1.param.uz * wf1_t2 + blData2.param.uz * wf2_t2;
  ut_d1 = blData1.param.uz * wf1_d1 + blData2.param.uz * wf2_d1;
  ut_d2 = blData1.param.uz * wf1_d2 + blData2.param.uz * wf2_d2;
  ut_u1 = blData1.param.uz * wf1_u1 + blData2.param.uz * wf2_u1 + wf1;
  ut_u2 = blData1.param.uz * wf1_u2 + blData2.param.uz * wf2_u2 + wf2;
  ut_ms = blData1.param.uz * wf1_ms + blData2.param.uz * wf2_ms;
  ut_re = blData1.param.uz * wf1_re + blData2.param.uz * wf2_re;
  ut_xf = blData1.param.uz * wf1_xf + blData2.param.uz * wf2_xf;

  //---- set primary "t" variables at xt  (really placed into "2" variables)
  blData2.param.xz = xt;
  blData2.param.tz = tt;
  blData2.param.dz = dt;
  blData2.param.uz = ut;

  blData2.param.amplz = amcrit;
  blData2.param.sz = 0.0;

  //---- calculate laminar secondary "t" variables
  blkin(boundaryLayerState);
  blData2 = blvar(blData2, FlowRegimeEnum::Laminar);

  //---- calculate x1-xt midpoint cfm value
  SkinFrictionCoefficients laminarSkinFriction =
      blmid(boundaryLayerState, FlowRegimeEnum::Laminar);

  //=    at this point, all "2" variables are really "t" variables at xt

  //---- set up newton system for dam, dth, dds, due, dxi  at  x1 and xt
  blc = blDiffSolver.solve(FlowRegimeEnum::Laminar, boundaryLayerState, laminarSkinFriction, amcrit);

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
  blData2 = blvar(blData2, FlowRegimeEnum::Turbulent);

  //---- set initial shear coefficient value st at transition point
  //-    ( note that cq2, cq2_t2, etc. are really "cqt", "cqt_tt", etc.)

  ctr = 1.8 * exp(-3.3 / (blData2.hkz.scalar - 1.0));
  ctr_hk2 = ctr * 3.3 / (blData2.hkz.scalar - 1.0) / (blData2.hkz.scalar - 1.0);

  st = ctr * blData2.cqz.scalar;
  st_tt =
      ctr * blData2.cqz.t() + blData2.cqz.scalar * ctr_hk2 * blData2.hkz.t();
  st_dt =
      ctr * blData2.cqz.d() + blData2.cqz.scalar * ctr_hk2 * blData2.hkz.d();
  st_ut =
      ctr * blData2.cqz.u() + blData2.cqz.scalar * ctr_hk2 * blData2.hkz.u();
  st_ms =
      ctr * blData2.cqz.ms() + blData2.cqz.scalar * ctr_hk2 * blData2.hkz.ms();
  st_re = ctr * blData2.cqz.re();

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

  blData2.param.amplz = 0.0;
  blData2.param.sz = st;

  //---- recalculate turbulent secondary "t" variables using proper cti
  blData2 = blvar(blData2, FlowRegimeEnum::Turbulent);

  stepbl(boundaryLayerState);
  restoreblData(2);

  //---- calculate xt-x2 midpoint cfm value
  SkinFrictionCoefficients turbulentSkinFriction =
      blmid(boundaryLayerState, FlowRegimeEnum::Turbulent);

  //---- set up newton system for dct, dth, dds, due, dxi  at  xt and x2
  blc = blDiffSolver.solve(FlowRegimeEnum::Turbulent, boundaryLayerState, turbulentSkinFriction, amcrit);

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
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      double dui = 0.0;
      for (int js = 1; js <= 2; js++) {
        for (int jbl = 0; jbl < nbl.get(js) - 1; ++jbl) {
          double ue_m = -vti.get(is)[ibl] * vti.get(js)[jbl] *
                        dij(ipan.get(is)[ibl],
                            ipan.get(js)[jbl]);
          dui += ue_m * mass.get(js)[jbl];
        }
      }
      uedg.get(is)[ibl] = uinv.get(is)[ibl] + dui;
    }
  }
  return true;
}


bool XFoil::uicalc() {
  //--------------------------------------------------------------
  //     sets inviscid ue from panel inviscid tangential velocity
  //--------------------------------------------------------------
  for (int is = 1; is <= 2; is++) {
    uinv.get(is)[0] = 0.0;
    uinv_a.get(is)[0] = 0.0;
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      int i = ipan.get(is)[ibl];
      uinv.get(is)[ibl] = vti.get(is)[ibl] * qinv[i];
      uinv_a.get(is)[ibl] = vti.get(is)[ibl] * qinv_a[i];
    }
  }

  return true;
}


bool XFoil::xicalc() {
  //-------------------------------------------------------------
  //     sets bl arc length array on each airfoil side and wake
  //-------------------------------------------------------------

  
    for (int ibl = 0; ibl <= iblte.top; ++ibl) {
      xssi.top[ibl] = sst - foil.foil_shape.spline_length[ipan.get(1)[ibl]];
    }
  
    for (int ibl = 0; ibl <= iblte.bottom; ++ibl) {
      xssi.bottom[ibl] = foil.foil_shape.spline_length[ipan.get(2)[ibl]] - sst;
    }

    // Wake: start from TE, duplicate TE value at first wake station
    xssi.bottom[iblte.bottom + 1] = xssi.bottom[iblte.bottom];
    for (int ibl = iblte.bottom + 2; ibl < nbl.bottom; ++ibl) {
      xssi.bottom[ibl] = xssi.bottom[ibl - 1] +
                          (foil.wake_shape.points.col(ipan.get(2)[ibl]) -
                           foil.wake_shape.points.col(ipan.get(2)[ibl - 1]))
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
      const int te_bot_0b = iblte.bottom; // 0-based TE for array indexing
      const double zn = 1.0 - (xssi.bottom[te_bot_0b + (iw0 + 1)] -
                               xssi.bottom[te_bot_0b]) /
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

  if (xstrip.get(is) >= 1.0) {
    return xssi.get(is)[iblte.get(is)];
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
    str = foil.edge.sle + (foil.foil_shape.spline_length[0] - foil.edge.sle) * xstrip.top;
  } else {
    str = foil.edge.sle + (foil.foil_shape.spline_length[foil.foil_shape.n - 1] - foil.edge.sle) * xstrip.bottom;
  }
  str = spline::sinvrt(str, xstrip.get(is), w1, w3, foil.foil_shape.spline_length.head(point_count), point_count);
  xiforc = std::min((str - sst), xssi.get(is)[iblte.get(is)]);
  if (xiforc < 0.0) {
    ss << " ***  stagnation point is past trip on side " << is << "\n";
    writeString(ss.str());

    xiforc = xssi.get(is)[iblte.get(is)];
  }

  return xiforc;
}
