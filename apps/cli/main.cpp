#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>

#include "XFoil.h"

struct IterationContext {
  int iterationLimit = 100;
  bool autoInitializeBoundaryLayer = true;
  int iterationCount = 0;
  bool encounteredErrors = false;
};

bool iterate(XFoil *xfoil, IterationContext &context) {
  if (!xfoil->viscal()) {
    xfoil->lvconv = false;
    std::cout
        << "CpCalc: local speed too large\nCompressibility corrections invalid"
        << std::endl;
    return false;
  }

  while (context.iterationCount < context.iterationLimit &&
         !xfoil->lvconv /*&& !s_bCancel*/) {
    if (xfoil->ViscousIter()) {
      // if (m_x0 && m_y0) {
      //  m_x0->append((double)m_Iterations);
      //  m_y0->append(xfoil->rmsbl);
      //}
      // if (m_x1 && m_y1) {
      //  m_x1->append((double)m_Iterations);
      //  m_y1->append(xfoil->rmxbl);
      //}
      context.iterationCount++;
    } else
      context.iterationCount = context.iterationLimit;
  }

  {
    const auto viscal_end = xfoil->ViscalEnd();
    xfoil->cpi = viscal_end.inviscidCp;
    xfoil->cpv = viscal_end.viscousCp;
  }

  if (context.iterationCount >= context.iterationLimit && !xfoil->lvconv) {
    if (context.autoInitializeBoundaryLayer) {
      xfoil->setBLInitialized(false);
      xfoil->lipan = false;
    }
    return true;
  }
  if (!xfoil->lvconv) {
    context.encounteredErrors = true;
    return false;
  } else {
    // converged at last
    return true;
  }
  return false;
}

int loadDatFile(std::string filename, double x[604], double y[604]) {
  std::ifstream fs(filename);
  if (!fs) {
    std::cout << "Failed to open dat file" << std::endl;
    return -1;
  }

  std::string line;
  std::getline(fs, line);
  std::cout << "Foil name : " << line << std::endl;
  int cnt = 0;
  while (!fs.eof()) {
    std::getline(fs, line);

    line.erase(0, line.find_first_not_of(" \t"));
    int endOfX = line.find_first_of(" \t");
    if (endOfX == -1) continue;

    std::string sx = line.substr(0, endOfX);
    std::string sy = line.substr(endOfX);

    x[cnt] = atof(sx.c_str());
    y[cnt] = atof(sy.c_str());
    cnt++;
  }
  return cnt;
}

int main() {
  double x[604], y[604];
  int n = 0;

  std::stringstream ss;
  auto start = std::chrono::steady_clock::now();

  // auto naca = new XFoil();
  // naca->naca4(4412, 100 / 2);
  // for (int i = 0; i < naca->nb; i++) {
  //   x[i] = naca->xb[i + 1];
  //   y[i] = naca->yb[i + 1];
  // }
  // n = naca->nb;

  n = loadDatFile("sample/CLARK_Y.dat", x, y);
  if (n == -1) return 1;

  XFoil *foil = new XFoil();

  if (!foil->initXFoilGeometry(n, x, y)) {
    std::cout << "Initialization error!" << std::endl;
    return 1;
  }
  foil->initXFoilAnalysis(100000, 0, 0.0, 9.0, 1.0, 1.0, XFoil::ReynoldsType::CONSTANT, XFoil::MachType::CONSTANT, true, ss);

  IterationContext iterationContext{};

  for (double alpha = 0; alpha < 15; alpha += 0.5) {
    iterationContext.iterationCount = 0;
    iterationContext.encounteredErrors = false;

    foil->setBLInitialized(false);
    foil->lipan = false;

    foil->setAlpha(alpha * 3.14159 / 180);
    foil->analysisState().controlByAlpha = true;
    foil->setQInf(1.0);
    std::cout << "alpha : " << alpha << std::endl;

    if (!foil->specal()) {
      std::cout << "Invalid Analysis Settings" << std::endl;
      return 1;
    }
    foil->lwake = false;
    foil->lvconv = false;

    while (!iterate(foil, iterationContext))
      ;

    // std::cout << ss.str() << std::endl;

    if (foil->lvconv) {
      std::cout << "  converged after " << iterationContext.iterationCount << " iterations"
                << std::endl;
      std::cout << "  cl : " << foil->cl << ", cd : " << foil->cd
                << ", cm : " << foil->cm << ", xcp : " << foil->getXcp()
                << std::endl;
    } else {
      std::cout << "  unconverged" << std::endl;
    }
  }

  auto end = std::chrono::steady_clock::now();
  // 経過時間を計算
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  std::cout << "Elapsed time: " << duration.count() << " milliseconds." << std::endl;

  return 0;
}
