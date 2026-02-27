#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <memory>
#include <sstream>

#include "XFoil.h"

namespace {
constexpr int kMaxDatPoints = 604;
constexpr int kIterationLimit = 100;
constexpr int kReynolds = 100000;
constexpr double kPi = 3.14159265358979323846;
}  // namespace

struct IterationContext {
  int iterationLimit = kIterationLimit;
  bool autoInitializeBoundaryLayer = true;
  int iterationCount = 0;
  bool encounteredErrors = false;
};

bool iterate(XFoil& xfoil, IterationContext& context) {
  if (!xfoil.viscal()) {
    xfoil.invalidateConvergedSolution();
    std::cout
        << "CpCalc: local speed too large\nCompressibility corrections invalid"
        << std::endl;
    return false;
  }

  while (context.iterationCount < context.iterationLimit &&
         !xfoil.hasConvergedSolution()) {
    if (xfoil.ViscousIter()) {
      context.iterationCount++;
    } else
      context.iterationCount = context.iterationLimit;
  }

  {
    const auto viscal_end = xfoil.ViscalEnd();
    xfoil.cpi = viscal_end.inviscidCp;
    xfoil.cpv = viscal_end.viscousCp;
  }

  if (context.iterationCount >= context.iterationLimit &&
      !xfoil.hasConvergedSolution()) {
    if (context.autoInitializeBoundaryLayer) {
      xfoil.setBLInitialized(false);
      xfoil.invalidatePanelMap();
    }
    return true;
  }
  if (!xfoil.hasConvergedSolution()) {
    context.encounteredErrors = true;
    return false;
  }
  // converged at last
  return true;
}

int loadDatFile(const std::string& filename, double x[kMaxDatPoints],
                double y[kMaxDatPoints]) {
  std::ifstream fs(filename);
  if (!fs) {
    std::cout << "Failed to open dat file" << std::endl;
    return -1;
  }

  std::string line;
  std::getline(fs, line);
  std::cout << "Foil name : " << line << std::endl;
  int cnt = 0;
  while (std::getline(fs, line)) {
    if (cnt >= kMaxDatPoints) {
      break;
    }
    std::istringstream iss(line);
    double x_coord = 0.0;
    double y_coord = 0.0;
    if (!(iss >> x_coord >> y_coord)) {
      continue;
    }
    x[cnt] = x_coord;
    y[cnt] = y_coord;
    ++cnt;
  }
  return cnt;
}

int main() {
  double x[kMaxDatPoints], y[kMaxDatPoints];
  int n = 0;

  const auto start = std::chrono::steady_clock::now();

  // auto naca = new XFoil();
  // naca->naca4(4412, 100 / 2);
  // for (int i = 0; i < naca->nb; i++) {
  //   x[i] = naca->xb[i + 1];
  //   y[i] = naca->yb[i + 1];
  // }
  // n = naca->nb;

  n = loadDatFile("sample/CLARK_Y.dat", x, y);
  if (n == -1) return 1;

  auto foil = std::make_unique<XFoil>();

  if (!foil->initXFoilGeometry(n, x, y)) {
    std::cout << "Initialization error!" << std::endl;
    return 1;
  }
  foil->initXFoilAnalysis(kReynolds, 0, 0.0, 9.0, 1.0, 1.0,
                          XFoil::ReynoldsType::CONSTANT,
                          XFoil::MachType::CONSTANT, true);

  IterationContext iterationContext{};

  for (double alpha = 0; alpha < 15; alpha += 0.5) {
    iterationContext.iterationCount = 0;
    iterationContext.encounteredErrors = false;

    foil->setBLInitialized(false);
    foil->invalidatePanelMap();

    foil->setAlpha(alpha * kPi / 180.0);
    foil->analysisState().controlByAlpha = true;
    foil->setQInf(1.0);
    std::cout << "alpha : " << alpha << std::endl;

    if (!foil->specal()) {
      std::cout << "Invalid Analysis Settings" << std::endl;
      return 1;
    }
    foil->invalidateWakeGeometry();
    foil->invalidateConvergedSolution();

    bool finished = false;
    while (!finished) {
      finished = iterate(*foil, iterationContext);
    }


    if (foil->hasConvergedSolution()) {
      std::cout << "  converged after " << iterationContext.iterationCount << " iterations"
                << std::endl;
      std::cout << "  cl : " << foil->Cl() << ", cd : " << foil->Cd()
                << ", cm : " << foil->Cm() << ", xcp : " << foil->getXcp()
                << std::endl;
    } else {
      std::cout << "  unconverged" << std::endl;
    }
  }

  const auto end = std::chrono::steady_clock::now();
  // 経過時間を計算
  const auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  std::cout << "Elapsed time: " << duration.count() << " milliseconds." << std::endl;

  return 0;
}
