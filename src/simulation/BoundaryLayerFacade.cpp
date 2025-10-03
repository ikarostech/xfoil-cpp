#include "simulation/BoundaryLayerFacade.hpp"

#include <utility>

#include "simulation/XFoil.h"

BoundaryLayerFacade::BoundaryLayerFacade(XFoil &owner) : xfoil_(owner) {}

void BoundaryLayerFacade::reset() {
  xfoil_.setBLInitialized(false);
  xfoil_.lvconv = false;
}

bool BoundaryLayerFacade::ensureInitialized() {
  if (xfoil_.isBLInitialized()) {
    return true;
  }
  if (!buildNewtonSystem()) {
    return false;
  }
  xfoil_.setBLInitialized(true);
  return true;
}

bool BoundaryLayerFacade::iterateOnce() {
  constexpr double kConvergenceTolerance = 0.0001;

  if (!buildNewtonSystem()) {
    return false;
  }

  if (!xfoil_.blsolve()) {
    return false;
  }

  if (!xfoil_.update()) {
    return false;
  }

  if (xfoil_.lalfa) {
    xfoil_.minf_cl = xfoil_.getActualMach(xfoil_.cl, xfoil_.mach_type);
    xfoil_.reinf_cl = xfoil_.getActualReynolds(xfoil_.cl, xfoil_.reynolds_type);
    if (!xfoil_.comset()) {
      return false;
    }
  } else {
    auto qiset_result = xfoil_.qiset();
    xfoil_.qinv = std::move(qiset_result.qinv);
    xfoil_.qinv_a = std::move(qiset_result.qinv_a);
    if (!xfoil_.uicalc()) {
      return false;
    }
  }

  xfoil_.qvis = xfoil_.qvfue();
  xfoil_.surface_vortex = xfoil_.gamqv();
  if (!xfoil_.stmove()) {
    return false;
  }

  const auto cl_result = xfoil_.clcalc(xfoil_.cmref);
  xfoil_.applyClComputation(cl_result);
  xfoil_.cd = xfoil_.cdcalc();

  if (xfoil_.rmsbl < kConvergenceTolerance) {
    xfoil_.lvconv = true;
    xfoil_.avisc = xfoil_.alfa;
    xfoil_.mvisc = xfoil_.minf;
    xfoil_.writeString("----------CONVERGED----------\n\n");
  }

  return true;
}

bool BoundaryLayerFacade::buildNewtonSystem() {
  return xfoil_.setbl();
}

bool BoundaryLayerFacade::initializeWithCurrentUe() {
  return xfoil_.mrchue();
}

bool BoundaryLayerFacade::marchDisplacement() {
  return xfoil_.mrchdu();
}

bool BoundaryLayerFacade::hasConverged() const {
  return xfoil_.lvconv;
}

