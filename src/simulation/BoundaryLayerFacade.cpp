#include "simulation/BoundaryLayerFacade.hpp"

#include <sstream>
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
  XFoil& xfoil = xfoil_;
#define getActualMach(...) xfoil.getActualMach(__VA_ARGS__)
#define getActualReynolds(...) xfoil.getActualReynolds(__VA_ARGS__)
#define comset() xfoil.comset()
#define ueset(...) xfoil.ueset(__VA_ARGS__)
#define swapEdgeVelocities(...) xfoil.swapEdgeVelocities(__VA_ARGS__)
#define computeLeTeSensitivities(...) xfoil.computeLeTeSensitivities(__VA_ARGS__)
#define writeString(...) xfoil.writeString(__VA_ARGS__)
#define clearDerivativeVectors(...) xfoil.clearDerivativeVectors(__VA_ARGS__)
#define xifset(...) xfoil.xifset(__VA_ARGS__)
#define dij(...) xfoil.dij(__VA_ARGS__)
#define blprv(...) xfoil.blprv(__VA_ARGS__)
#define blkin(...) xfoil.blkin(__VA_ARGS__)
#define trchek(...) xfoil.trchek(__VA_ARGS__)
#define tesys(...) xfoil.tesys(__VA_ARGS__)
#define blsys() xfoil.blsys()
#define blvar(...) xfoil.blvar(__VA_ARGS__)
#define restoreblData(...) xfoil.restoreblData(__VA_ARGS__)
#define mrchue() xfoil.mrchue()
#define mrchdu() xfoil.mrchdu()
#define stepbl() xfoil.stepbl()
#define blmid(...) xfoil.blmid(__VA_ARGS__)
  auto& lalfa = xfoil.lalfa;
  auto& cl = xfoil.cl;
  auto& clspec = xfoil.clspec;
  auto& mach_type = xfoil.mach_type;
  auto& reynolds_type = xfoil.reynolds_type;
  auto& minf = xfoil.minf;
  auto& gm1bl = xfoil.gm1bl;
  auto& gamm1 = xfoil.gamm1;
  auto& qinfbl = xfoil.qinfbl;
  auto& qinf = xfoil.qinf;
  auto& tklam = xfoil.tklam;
  auto& tkl_msq = xfoil.tkl_msq;
  auto& tkbl = xfoil.tkbl;
  auto& tkbl_ms = xfoil.tkbl_ms;
  auto& rstbl = xfoil.rstbl;
  auto& rstbl_ms = xfoil.rstbl_ms;
  auto& hstinv = xfoil.hstinv;
  auto& hstinv_ms = xfoil.hstinv_ms;
  auto& reinf = xfoil.reinf;
  auto& hvrat = xfoil.hvrat;
  auto& reybl = xfoil.reybl;
  auto& reybl_re = xfoil.reybl_re;
  auto& reybl_ms = xfoil.reybl_ms;
  auto& amcrit = xfoil.amcrit;
  auto& acrit = xfoil.acrit;
  auto& lblini = xfoil.lblini;
  auto& uedg = xfoil.uedg;
  auto& isys = xfoil.isys;
  auto& iblte = xfoil.iblte;
  auto& uinv_a = xfoil.uinv_a;
  auto& uinv = xfoil.uinv;
  auto& xiforc = xfoil.xiforc;
  auto& tran = xfoil.tran;
  auto& turb = xfoil.turb;
  auto& simi = xfoil.simi;
  auto& wake = xfoil.wake;
  auto& xssi = xfoil.xssi;
  auto& ctau = xfoil.ctau;
  auto& itran = xfoil.itran;
  auto& mass = xfoil.mass;
  auto& vti = xfoil.vti;
  auto& blData2 = xfoil.blData2;
  auto& blData1 = xfoil.blData1;
  auto& thet = xfoil.thet;
  auto& dstr = xfoil.dstr;
  auto& ante = xfoil.ante;
  auto& ctq = xfoil.ctq;
  auto& sst_go = xfoil.sst_go;
  auto& sst_gp = xfoil.sst_gp;
  auto& nbl = xfoil.nbl;
  auto& vm = xfoil.vm;
  auto& va = xfoil.va;
  auto& vb = xfoil.vb;
  auto& vdel = xfoil.vdel;
  auto& blc = xfoil.blc;
  auto& wgap = xfoil.wgap;
  auto& sharp = xfoil.sharp;
  auto& foil = xfoil.foil;
  auto& qinv = xfoil.qinv;
  auto& qinv_a = xfoil.qinv_a;
  auto& ipan = xfoil.ipan;
  auto& nsys = xfoil.nsys;
  auto& vz = xfoil.vz;
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
  gm1bl = gamm1;

  //---- set parameters for compressibility correction
  qinfbl = qinf;
  tkbl = tklam;
  tkbl_ms = tkl_msq;

  //---- stagnation density and 1/enthalpy
  rstbl = pow((1.0 + 0.5 * gm1bl * minf * minf), (1.0 / gm1bl));
  rstbl_ms = 0.5 * rstbl / (1.0 + 0.5 * gm1bl * minf * minf);
  hstinv = gm1bl * (minf / qinfbl) * (minf / qinfbl) /
           (1.0 + 0.5 * gm1bl * minf * minf);
  hstinv_ms = gm1bl * (1.0 / qinfbl) * (1.0 / qinfbl) /
                  (1.0 + 0.5 * gm1bl * minf * minf) -
              0.5 * gm1bl * hstinv / (1.0 + 0.5 * gm1bl * minf * minf);

  //---- set reynolds number based on freestream density, velocity, viscosity
  herat = 1.0 - 0.5 * qinfbl * qinfbl * hstinv;
  herat_ms = -0.5 * qinfbl * qinfbl * hstinv_ms;

  reybl = reinf * sqrt(herat * herat * herat) * (1.0 + hvrat) / (herat + hvrat);
  reybl_re = sqrt(herat * herat * herat) * (1.0 + hvrat) / (herat + hvrat);
  reybl_ms = reybl * (1.5 / herat - 1.0 / (herat + hvrat)) * herat_ms;

  amcrit = acrit;

  if (!lblini) {
    //----- initialize bl by marching with ue (fudge at separation)
    // TRACE(" initializing bl ...\n");
    writeString("   Initializing bl ...\n");

  initializeWithCurrentUe();
    lblini = true;
  }

  //---- march bl with current ue and ds to establish transition
  marchDisplacement();

  XFoil::SidePair<VectorXd> usav = uedg;

  ueset();
  const auto swapped_edge_velocities = swapEdgeVelocities(usav);
  usav = swapped_edge_velocities.swappedUsav;
  uedg = swapped_edge_velocities.restoredUedg;
jvte1 = isys.top[iblte.top];
jvte2 = isys.bottom[iblte.bottom];

  dule1 = uedg.top[0] - usav.top[0];
  dule2 = uedg.bottom[0] - usav.bottom[0];

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

  //*** go over each boundary layer/wake
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
    xiforc = xifset(is);

    tran = false;
    turb = false;

    //**** sweep downstream setting up bl equation linearizations
    for (int ibl = 0; ibl < nbl.get(is) - 1; ++ibl) {
      
      int iv = isys.get(is)[ibl];

      simi = (ibl == 0);
      wake = (ibl > iblte.get(is));
      tran = (ibl == itran.get(is));
      turb = (ibl > itran.get(is));

      //---- set primary variables for current station
      xsi = xssi.get(is)[ibl];
      if (ibl < itran.get(is))
        ami = ctau.get(is)[ibl];
      else
        cti = ctau.get(is)[ibl];
      uei = uedg.get(is)[ibl];
      thi = thet.get(is)[ibl];
      mdi = mass.get(is)[ibl];

      dsi = mdi / uei;

      if (wake) {
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

      //---- "forced" changes due to mismatch between uedg and
      // usav=uinv+dij*mass
      due2 = uedg.get(is)[ibl] - usav.get(is)[ibl];
      dds2 = d2_u2 * due2;

      blprv(xsi, ami, cti, thi, dsi, dswaki, uei); // cti
      blkin();

      //---- check for transition and set tran, xt, etc. if found
      if (tran) {
        trchek();
        ami = blData2.param.amplz;
      }

      if (ibl == itran.get(is) && !tran) {
        // TRACE("setbl: xtr???  n1=%d n2=%d: \n", ampl1, ampl2);

        ss << "setbl: xtr???  n1=" << blData1.param.amplz
           << " n2=" << blData2.param.amplz << ":\n";
        writeString(ss.str());
        ss.str("");
      }

      //---- assemble 10x4 linearized system for dctau, dth, dds, due, dxi
      //	   at the previous "1" station and the current "2" station

      if (ibl == iblte.get(is) + 1) {
        //----- define quantities at start of wake, adding te base thickness to
        // dstar
        tte = thet.get(1)[iblte.top] +
              thet.get(2)[iblte.bottom];
        dte = dstr.get(1)[iblte.top] +
              dstr.get(2)[iblte.bottom] + ante;
        cte = (ctau.get(1)[iblte.top] *
                   thet.get(1)[iblte.top] +
               ctau.get(2)[iblte.bottom] *
                   thet.get(2)[iblte.bottom]) /
               tte;
        tesys(cte, tte, dte);

        tte_tte1 = 1.0;
        tte_tte2 = 1.0;
        dte_mte1 = 1.0 / uedg.top[iblte.top];
        dte_ute1 = -dstr.get(1)[iblte.top] /
                    uedg.top[iblte.top];
        dte_mte2 = 1.0 / uedg.bottom[iblte.bottom];
        dte_ute2 = -dstr.get(2)[iblte.bottom] /
                    uedg.bottom[iblte.bottom];
        cte_cte1 = thet.get(1)[iblte.top] / tte;
        cte_cte2 = thet.get(2)[iblte.bottom] / tte;
        cte_tte1 = (ctau.get(1)[iblte.top] - cte) / tte;
        cte_tte2 = (ctau.get(2)[iblte.bottom] - cte) / tte;

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

        //----- "forced" changes from  uedg --- usav=uinv+dij*mass	mismatch
        due1 = 0.0;
        dds1 =
            dte_ute1 *
                (uedg.top[iblte.top] -
                 usav.top[iblte.top]) +
            dte_ute2 *
                (uedg.bottom[iblte.bottom] -
                 usav.bottom[iblte.bottom]);
      } else {
        blsys();
      }

      //---- save wall shear and equil. max shear coefficient for plotting
      // output
      ctq.get(is)[ibl] = blData2.cqz.scalar;

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
        vm[0][jv][iv] = blc.a1(0, 2) * d1_m[jv] + blc.a1(0, 3) * u1_m[jv] +
                        blc.a2(0, 2) * d2_m[jv] + blc.a2(0, 3) * u2_m[jv] +
                        (blc.a1(0, 4) + blc.a2(0, 4) + blc.d_xi[0]) *
                            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
      }

      vb[iv](0, 0) = blc.a1(0, 0);
      vb[iv](0, 1) = blc.a1(0, 1);

      va[iv](0, 0) = blc.a2(0, 0);
      va[iv](0, 1) = blc.a2(0, 1);

      if (lalfa)
        vdel[iv](0, 1) = blc.d_re[0] * re_clmr + blc.d_msq[0] * msq_clmr;
      else
        vdel[iv](0, 1) = (blc.a1(0, 3) * u1_a + blc.a1(0, 2) * d1_a) +
                         (blc.a2(0, 3) * u2_a + blc.a2(0, 2) * d2_a) +
                         (blc.a1(0, 4) + blc.a2(0, 4) + blc.d_xi[0]) *
                             (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

      vdel[iv](0, 0) = blc.rhs[0] + (blc.a1(0, 3) * due1 + blc.a1(0, 2) * dds1) +
                       (blc.a2(0, 3) * due2 + blc.a2(0, 2) * dds2) +
                       (blc.a1(0, 4) + blc.a2(0, 4) + blc.d_xi[0]) *
                           (xi_ule1 * dule1 + xi_ule2 * dule2);

      for (int jv = 1; jv <= nsys; jv++) {
        vm[1][jv][iv] = blc.a1(1, 2) * d1_m[jv] + blc.a1(1, 3) * u1_m[jv] +
                        blc.a2(1, 2) * d2_m[jv] + blc.a2(1, 3) * u2_m[jv] +
                        (blc.a1(1, 4) + blc.a2(1, 4) + blc.d_xi[1]) *
                            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
      }
      vb[iv](1, 0) = blc.a1(1, 0);
      vb[iv](1, 1) = blc.a1(1, 1);

      va[iv](1, 0) = blc.a2(1, 0);
      va[iv](1, 1) = blc.a2(1, 1);

      if (lalfa)
        vdel[iv](1, 1) = blc.d_re[1] * re_clmr + blc.d_msq[1] * msq_clmr;
      else
        vdel[iv](1, 1) = (blc.a1(1, 3) * u1_a + blc.a1(1, 2) * d1_a) +
                         (blc.a2(1, 3) * u2_a + blc.a2(1, 2) * d2_a) +
                         (blc.a1(1, 4) + blc.a2(1, 4) + blc.d_xi[1]) *
                             (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

      vdel[iv](1, 0) = blc.rhs[1] + (blc.a1(1, 3) * due1 + blc.a1(1, 2) * dds1) +
                       (blc.a2(1, 3) * due2 + blc.a2(1, 2) * dds2) +
                       (blc.a1(1, 4) + blc.a2(1, 4) + blc.d_xi[1]) *
                           (xi_ule1 * dule1 + xi_ule2 * dule2);

      // memory overlap problem
      for (int jv = 1; jv <= nsys; jv++) {
        vm[2][jv][iv] = blc.a1(2, 2) * d1_m[jv] + blc.a1(2, 3) * u1_m[jv] +
                        blc.a2(2, 2) * d2_m[jv] + blc.a2(2, 3) * u2_m[jv] +
                        (blc.a1(2, 4) + blc.a2(2, 4) + blc.d_xi[2]) *
                            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
      }

      vb[iv](2, 0) = blc.a1(2, 0);
      vb[iv](2, 1) = blc.a1(2, 1);

      va[iv](2, 0) = blc.a2(2, 0);
      va[iv](2, 1) = blc.a2(2, 1);

      if (lalfa)
        vdel[iv](2, 1) = blc.d_re[2] * re_clmr + blc.d_msq[2] * msq_clmr;
      else
        vdel[iv](2, 1) = (blc.a1(2, 3) * u1_a + blc.a1(2, 2) * d1_a) +
                         (blc.a2(2, 3) * u2_a + blc.a2(2, 2) * d2_a) +
                         (blc.a1(2, 4) + blc.a2(2, 4) + blc.d_xi[2]) *
                             (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

      vdel[iv](2, 0) = blc.rhs[2] + (blc.a1(2, 3) * due1 + blc.a1(2, 2) * dds1) +
                       (blc.a2(2, 3) * due2 + blc.a2(2, 2) * dds2) +
                       (blc.a1(2, 4) + blc.a2(2, 4) + blc.d_xi[2]) *
                           (xi_ule1 * dule1 + xi_ule2 * dule2);

      if (ibl == iblte.get(is) + 1) {
        //----- redefine coefficients for tte, dte, etc
        vz[0][0] = blc.a1(0, 0) * cte_cte1;
        vz[0][1] = blc.a1(0, 0) * cte_tte1 + blc.a1(0, 1) * tte_tte1;
        vb[iv](0, 0) = blc.a1(0, 0) * cte_cte2;
        vb[iv](0, 1) = blc.a1(0, 0) * cte_tte2 + blc.a1(0, 1) * tte_tte2;

        vz[1][0] = blc.a1(1, 0) * cte_cte1;
        vz[1][1] = blc.a1(1, 0) * cte_tte1 + blc.a1(1, 1) * tte_tte1;
        vb[iv](1, 0) = blc.a1(1, 0) * cte_cte2;
        vb[iv](1, 1) = blc.a1(1, 0) * cte_tte2 + blc.a1(1, 1) * tte_tte2;

        vz[2][0] = blc.a1(2, 0) * cte_cte1;
        vz[2][1] = blc.a1(2, 0) * cte_tte1 + blc.a1(2, 1) * tte_tte1;
        vb[iv](2, 0) = blc.a1(2, 0) * cte_cte2;
        vb[iv](2, 1) = blc.a1(2, 0) * cte_tte2 + blc.a1(2, 1) * tte_tte2;
      }

      //---- turbulent intervals will follow if currently at transition interval
      if (tran) {
        turb = true;

        //------ save transition location
        itran.get(is) = ibl;
      }

      tran = false;

      if (ibl == iblte.get(is)) {
        //----- set "2" variables at te to wake correlations for next station

        turb = true;
        wake = true;
        blvar(blData2, FlowRegimeEnum::Wake);
        blmid(FlowRegimeEnum::Wake);
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
      stepbl();

      //---- next streamwise station
    }

    //---- next airfoil side
  }

  return true;

#undef getActualMach
#undef getActualReynolds
#undef comset
#undef ueset
#undef swapEdgeVelocities
#undef computeLeTeSensitivities
#undef writeString
#undef clearDerivativeVectors
#undef xifset
#undef dij
#undef blprv
#undef blkin
#undef trchek
#undef tesys
#undef blsys
#undef blvar
#undef restoreblData
#undef mrchue
#undef mrchdu
#undef stepbl
#undef blmid
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
