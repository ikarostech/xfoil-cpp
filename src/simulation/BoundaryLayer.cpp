#include "BoundaryLayer.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

#include "XFoil.h"

using BoundaryContext = BoundaryLayerWorkflow::MixedModeStationContext;

bool BoundaryLayerWorkflow::isStartOfWake(const XFoil& xfoil, int side,
                                          int stationIndex) {
  return stationIndex == lattice.trailingEdgeIndex.get(side) + 1;
}

void BoundaryLayerWorkflow::updateSystemMatricesForStation(
    XFoil& xfoil, int side, int stationIndex, BoundaryContext& ctx) {
  if (isStartOfWake(xfoil, side, stationIndex)) {
    ctx.tte = lattice.thet.get(1)[lattice.trailingEdgeIndex.top] +
              lattice.thet.get(2)[lattice.trailingEdgeIndex.bottom];
    ctx.dte = lattice.dstr.get(1)[lattice.trailingEdgeIndex.top] +
              lattice.dstr.get(2)[lattice.trailingEdgeIndex.bottom] + xfoil.foil.edge.ante;
    ctx.cte =
        (lattice.ctau.get(1)[lattice.trailingEdgeIndex.top] * lattice.thet.get(1)[lattice.trailingEdgeIndex.top] +
         lattice.ctau.get(2)[lattice.trailingEdgeIndex.bottom] *
             lattice.thet.get(2)[lattice.trailingEdgeIndex.bottom]) /
        ctx.tte;
    tesys(xfoil, ctx.cte, ctx.tte, ctx.dte);
  } else {
    xfoil.blsys(xfoil.boundaryLayerState, lattice);
  }
}

void BoundaryLayerWorkflow::initializeFirstIterationState(
    XFoil& xfoil, int side, int stationIndex, int previousTransition,
    BoundaryContext& ctx, double& ueref, double& hkref, double& ami) {
  ueref = xfoil.blData2.param.uz;
  hkref = xfoil.blData2.hkz.scalar;

  const bool inLaminarWindow =
      stationIndex < lattice.transitionIndex.get(side) && stationIndex >= previousTransition;
  if (inLaminarWindow) {
    double uem;
    double dsm;
    double thm;
    if (stationIndex > 0) {
      uem = lattice.uedg.get(side)[stationIndex - 1];
      dsm = lattice.dstr.get(side)[stationIndex - 1];
      thm = lattice.thet.get(side)[stationIndex - 1];
    } else {
      uem = lattice.uedg.get(side)[stationIndex];
      dsm = lattice.dstr.get(side)[stationIndex];
      thm = lattice.thet.get(side)[stationIndex];
    }
    const double uem_sq = uem * uem;
    const double msq =
        uem_sq * xfoil.hstinv /
        (xfoil.gm1bl * (1.0 - 0.5 * uem_sq * xfoil.hstinv));
    const auto hkin_result =
        boundary_layer::hkin(dsm / thm, msq);
    hkref = hkin_result.hk;
  }

  if (stationIndex < previousTransition) {
    if (xfoil.tran) {
      lattice.ctau.get(side)[stationIndex] = 0.03;
    }
    if (xfoil.turb) {
      const double prev =
          (stationIndex >= 1) ? lattice.ctau.get(side)[stationIndex - 1]
                              : lattice.ctau.get(side)[stationIndex];
      lattice.ctau.get(side)[stationIndex] = prev;
    }
    if (xfoil.tran || xfoil.turb) {
      ctx.cti = lattice.ctau.get(side)[stationIndex - 1];
      xfoil.blData2.param.sz = ctx.cti;
    }
  }
}

void BoundaryLayerWorkflow::configureSimilarityRow(XFoil& xfoil,
                                                   double ueref) {
  xfoil.blc.a2(3, 0) = 0.0;
  xfoil.blc.a2(3, 1) = 0.0;
  xfoil.blc.a2(3, 2) = 0.0;
  xfoil.blc.a2(3, 3) = xfoil.blData2.param.uz_uei;
  xfoil.blc.rhs[3] = ueref - xfoil.blData2.param.uz;
}

void BoundaryLayerWorkflow::configureViscousRow(XFoil& xfoil, double hkref,
                                                double ueref, double senswt,
                                                bool resetSensitivity,
                                                bool averageSensitivity,
                                                double& sens, double& sennew) {
  xfoil.blc.a2(3, 0) = 0.0;
  xfoil.blc.a2(3, 1) = xfoil.blData2.hkz.t();
  xfoil.blc.a2(3, 2) = xfoil.blData2.hkz.d();
  xfoil.blc.a2(3, 3) = xfoil.blData2.hkz.u() * xfoil.blData2.param.uz_uei;
  xfoil.blc.rhs[3] = 1.0;

  const double delta_sen =
      xfoil.blc.a2.block(0, 0, 4, 4).fullPivLu().solve(xfoil.blc.rhs)[3];

  sennew = senswt * delta_sen * hkref / ueref;
  if (resetSensitivity) {
    sens = sennew;
  } else if (averageSensitivity) {
    sens = 0.5 * (sens + sennew);
  }

  xfoil.blc.a2(3, 1) = xfoil.blData2.hkz.t() * hkref;
  xfoil.blc.a2(3, 2) = xfoil.blData2.hkz.d() * hkref;
  xfoil.blc.a2(3, 3) =
      (xfoil.blData2.hkz.u() * hkref + sens / ueref) * xfoil.blData2.param.uz_uei;
  xfoil.blc.rhs[3] =
      -(hkref * hkref) * (xfoil.blData2.hkz.scalar / hkref - 1.0) -
      sens * (xfoil.blData2.param.uz / ueref - 1.0);
}

bool BoundaryLayerWorkflow::applyMixedModeNewtonStep(
    XFoil& xfoil, int side, int stationIndex, double deps, double& ami,
    BoundaryContext& ctx) {
  xfoil.blc.rhs =
      xfoil.blc.a2.block(0, 0, 4, 4).fullPivLu().solve(xfoil.blc.rhs);

  ctx.dmax = std::max(std::fabs(xfoil.blc.rhs[1] / ctx.thi),
                      std::fabs(xfoil.blc.rhs[2] / ctx.dsi));
  if (stationIndex >= lattice.transitionIndex.get(side)) {
    ctx.dmax = std::max(ctx.dmax,
                        std::fabs(xfoil.blc.rhs[0] / (10.0 * ctx.cti)));
  }

  xfoil.rlx = 1.0;
  if (ctx.dmax > 0.3) {
    xfoil.rlx = 0.3 / ctx.dmax;
  }

  if (stationIndex < lattice.transitionIndex.get(side)) {
    ami += xfoil.rlx * xfoil.blc.rhs[0];
    ctx.ami = ami;
  }
  if (stationIndex >= lattice.transitionIndex.get(side)) {
    ctx.cti += xfoil.rlx * xfoil.blc.rhs[0];
  }
  ctx.thi += xfoil.rlx * xfoil.blc.rhs[1];
  ctx.dsi += xfoil.rlx * xfoil.blc.rhs[2];
  ctx.uei += xfoil.rlx * xfoil.blc.rhs[3];

  if (stationIndex >= lattice.transitionIndex.get(side)) {
    ctx.cti = std::clamp(ctx.cti, 0.0000001, 0.30);
  }

  const double hklim =
      (stationIndex <= lattice.trailingEdgeIndex.get(side)) ? 1.02 : 1.00005;
  const double uei_sq = ctx.uei * ctx.uei;
  const double msq = uei_sq * xfoil.hstinv /
                     (xfoil.gm1bl * (1.0 - 0.5 * uei_sq * xfoil.hstinv));
  double dsw = ctx.dsi - ctx.dswaki;
  xfoil.dslim(dsw, ctx.thi, msq, hklim);
  ctx.dsi = dsw + ctx.dswaki;

  return ctx.dmax <= deps;
}

bool BoundaryLayerWorkflow::iblpan(XFoil& xfoil) {
  std::stringstream ss;
  const int point_count = xfoil.foil.foil_shape.n;

  for (int i = 0; i <= xfoil.i_stagnation; i++) {
    lattice.stationToPanel.top[i] = xfoil.i_stagnation - i;
    lattice.vti.top[i] = 1.0;
  }

  lattice.trailingEdgeIndex.top = xfoil.i_stagnation;
  lattice.stationCount.top = lattice.trailingEdgeIndex.top + 2;

  for (int index = 0; index <= point_count - xfoil.i_stagnation; ++index) {
    lattice.stationToPanel.bottom[index] = xfoil.i_stagnation + 1 + index;
    lattice.vti.bottom[index] = -1.0;
  }

  lattice.trailingEdgeIndex.bottom = point_count - xfoil.i_stagnation - 2;

  for (int iw = 0; iw < xfoil.nw; iw++) {
    const int panel = point_count + iw;
    const int index = lattice.trailingEdgeIndex.bottom + iw + 2;
    lattice.stationToPanel.bottom[index - 1] = panel;
    lattice.vti.bottom[index - 1] = -1.0;
  }

  lattice.stationCount.bottom = lattice.trailingEdgeIndex.bottom + xfoil.nw + 2;

  for (int iw = 0; iw < xfoil.nw; iw++) {
    lattice.stationToPanel.top[lattice.trailingEdgeIndex.top + iw + 1] =
        lattice.stationToPanel.bottom[lattice.trailingEdgeIndex.bottom + iw + 1];
    lattice.vti.top[lattice.trailingEdgeIndex.top + iw + 1] = 1.0;
  }

  const int iblmax = std::max(lattice.trailingEdgeIndex.top, lattice.trailingEdgeIndex.bottom) +
                     xfoil.nw + 2;
  if (iblmax > IVX) {
    ss << "iblpan :  ***  bl array overflow\n";
    ss << "Increase IVX to at least " << iblmax << "\n";
    xfoil.writeString(ss.str());
    return false;
  }

  xfoil.lipan = true;
  return true;
}

bool BoundaryLayerWorkflow::iblsys(XFoil& xfoil) {
  int iv = 0;
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 0; ibl < lattice.stationCount.get(is) - 1; ++ibl) {
      ++iv;
      lattice.stationToSystem.get(is)[ibl] = iv;
    }
  }

  xfoil.nsys = iv;
  if (xfoil.nsys > 2 * IVX) {
    xfoil.writeString("*** iblsys: bl system array overflow. ***");
    return false;
  }

  return true;
}

bool BoundaryLayerWorkflow::stfind(XFoil& xfoil) {
  int stagnation_index = 0;
  bool found = false;
  const int point_count = xfoil.foil.foil_shape.n;
  for (int i = 0; i < point_count - 1; ++i) {
    if (xfoil.surface_vortex(0, i) >= 0.0 &&
        xfoil.surface_vortex(0, i + 1) < 0.0) {
      stagnation_index = i;
      found = true;
      break;
    }
  }

  if (!found) {
    xfoil.writeString("stfind: Stagnation point not found. Continuing ...\n");
    stagnation_index = point_count / 2;
  }

  xfoil.i_stagnation = stagnation_index;
  const double dgam = xfoil.surface_vortex(0, stagnation_index + 1) -
                      xfoil.surface_vortex(0, stagnation_index);
  const double ds = xfoil.foil.foil_shape.spline_length[stagnation_index + 1] -
                    xfoil.foil.foil_shape.spline_length[stagnation_index];

  if (xfoil.surface_vortex(0, stagnation_index) <
      -xfoil.surface_vortex(0, stagnation_index + 1)) {
    xfoil.sst = xfoil.foil.foil_shape.spline_length[stagnation_index] -
                ds * (xfoil.surface_vortex(0, stagnation_index) / dgam);
  } else {
    xfoil.sst =
        xfoil.foil.foil_shape.spline_length[stagnation_index + 1] -
        ds * (xfoil.surface_vortex(0, stagnation_index + 1) / dgam);
  }

  if (xfoil.sst <= xfoil.foil.foil_shape.spline_length[stagnation_index])
    xfoil.sst =
        xfoil.foil.foil_shape.spline_length[stagnation_index] + 0.0000001;
  if (xfoil.sst >= xfoil.foil.foil_shape.spline_length[stagnation_index + 1])
    xfoil.sst = xfoil.foil.foil_shape.spline_length[stagnation_index + 1] -
                0.0000001;

  xfoil.sst_go =
      (xfoil.sst - xfoil.foil.foil_shape.spline_length[stagnation_index + 1]) /
      dgam;
  xfoil.sst_gp =
      (xfoil.foil.foil_shape.spline_length[stagnation_index] - xfoil.sst) /
      dgam;

  return true;
}

bool BoundaryLayerWorkflow::stmove(XFoil& xfoil) {
  const int previous = xfoil.i_stagnation;
  stfind(xfoil);

  if (previous == xfoil.i_stagnation) {
    xfoil.xicalc();
  } else {
    iblpan(xfoil);
    xfoil.uicalc();
    xfoil.xicalc();
    iblsys(xfoil);

    if (xfoil.i_stagnation > previous) {
      const int delta = xfoil.i_stagnation - previous;

      lattice.transitionIndex.top += delta;
      lattice.transitionIndex.bottom -= delta;

      for (int ibl = lattice.stationCount.top - 2; ibl >= delta; --ibl) {
        lattice.ctau.top[ibl] = lattice.ctau.top[ibl - delta];
        lattice.thet.top[ibl] = lattice.thet.top[ibl - delta];
        lattice.dstr.top[ibl] = lattice.dstr.top[ibl - delta];
        lattice.uedg.top[ibl] = lattice.uedg.top[ibl - delta];
      }

      const double dudx =
          lattice.uedg.top[delta] / lattice.xssi.top[delta];
      for (int ibl = delta; ibl >= 1; --ibl) {
        lattice.ctau.top[ibl - 1] = lattice.ctau.top[delta];
        lattice.thet.top[ibl - 1] = lattice.thet.top[delta];
        lattice.dstr.top[ibl - 1] = lattice.dstr.top[delta];
        lattice.uedg.top[ibl - 1] = dudx * lattice.xssi.top[ibl - 1];
      }

      for (int ibl = 0; ibl < lattice.stationCount.bottom - 1; ++ibl) {
        lattice.ctau.bottom[ibl] = lattice.ctau.bottom[ibl + delta];
        lattice.thet.bottom[ibl] = lattice.thet.bottom[ibl + delta];
        lattice.dstr.bottom[ibl] = lattice.dstr.bottom[ibl + delta];
        lattice.uedg.bottom[ibl] = lattice.uedg.bottom[ibl + delta];
      }
    } else {
      const int delta = previous - xfoil.i_stagnation;

      lattice.transitionIndex.top -= delta;
      lattice.transitionIndex.bottom += delta;

      for (int ibl = lattice.stationCount.bottom - 1; ibl >= delta + 1; --ibl) {
        lattice.ctau.bottom[ibl - 1] =
            lattice.ctau.bottom[(ibl - delta) - 1];
        lattice.thet.bottom[ibl - 1] =
            lattice.thet.bottom[(ibl - delta) - 1];
        lattice.dstr.bottom[ibl - 1] =
            lattice.dstr.bottom[(ibl - delta) - 1];
        lattice.uedg.bottom[ibl - 1] =
            lattice.uedg.bottom[(ibl - delta) - 1];
      }

      const double dudx =
          lattice.uedg.bottom[delta] / lattice.xssi.bottom[delta];
      for (int ibl = delta; ibl >= 1; --ibl) {
        lattice.ctau.bottom[ibl - 1] = lattice.ctau.bottom[delta];
        lattice.thet.bottom[ibl - 1] = lattice.thet.bottom[delta];
        lattice.dstr.bottom[ibl - 1] = lattice.dstr.bottom[delta];
        lattice.uedg.bottom[ibl - 1] = dudx * lattice.xssi.bottom[ibl - 1];
      }

      for (int ibl = 0; ibl < lattice.stationCount.top - 1; ++ibl) {
        lattice.ctau.top[ibl] = lattice.ctau.top[ibl + delta];
        lattice.thet.top[ibl] = lattice.thet.top[ibl + delta];
        lattice.dstr.top[ibl] = lattice.dstr.top[ibl + delta];
        lattice.uedg.top[ibl] = lattice.uedg.top[ibl + delta];
      }
    }
  }

  for (int is = 1; is <= 2; ++is) {
    for (int ibl = 0; ibl < lattice.stationCount.get(is) - 1; ++ibl) {
      lattice.mass.get(is)[ibl] =
          lattice.dstr.get(is)[ibl] * lattice.uedg.get(is)[ibl];
    }
  }

  return true;
}

bool BoundaryLayerWorkflow::tesys(XFoil& xfoil, double cte, double tte, double dte) {
  xfoil.blc.clear();

  xfoil.blData2 = xfoil.blvar(xfoil.blData2, FlowRegimeEnum::Wake);

  xfoil.blc.a1(0, 0) = -1.0;
  xfoil.blc.a2(0, 0) = 1.0;
  xfoil.blc.rhs[0] = cte - xfoil.blData2.param.sz;

  xfoil.blc.a1(1, 1) = -1.0;
  xfoil.blc.a2(1, 1) = 1.0;
  xfoil.blc.rhs[1] = tte - xfoil.blData2.param.tz;

  xfoil.blc.a1(2, 2) = -1.0;
  xfoil.blc.a2(2, 2) = 1.0;
  xfoil.blc.rhs[2] = dte - xfoil.blData2.param.dz - xfoil.blData2.param.dwz;

  return true;
}
