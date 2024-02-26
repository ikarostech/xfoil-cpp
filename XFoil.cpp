/****************************************************************************

        XFoil Class
        Copyright (C) 2000 Mark Drela
        Copyright (C) 2003 Andre Deperrois techwinder@gmail.com

        This program is free software; you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation; either version 2 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program; if not, write to the Free Software
        Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

*****************************************************************************/

#include "XFoil.h"
#include <cstring>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/StdVector"
using namespace Eigen;
// determinant
double cross2(const Eigen::Vector2d& a, const Eigen::Vector2d& b)
{
  return a[0]*b[1] - a[1]*b[0];
}

#define PI 3.141592654

bool XFoil::s_bCancel = false;
bool XFoil::s_bFullReport = false;
double XFoil::vaccel = 0.01;
const int INDEX_START_WITH = 1;



XFoil::XFoil() {
  m_pOutStream = NULL;
  //------ primary dimensioning limit parameters

  //------ derived dimensioning limit parameters
  //	nax=800;//number of points in stored polar
  //	npx=8;//number of polars and reference polars
  //	nfx=128;// number of points in one reference polar
  //	ncom = 73;

  // imx   number of complex mapping coefficients  cn

  sccon = 5.6;
  gacon = 6.70;
  gbcon = 0.75;
  gbc0 = 0.60;
  gbc1 = 0.40;
  gccon = 18.0;
  dlcon = 0.9;
  ctcon = 0.01485111754659538130244;  //(ctcon = 0.5/(gacon**2 * gbcon))

  // fortran seems to initializes variables to 0
  mvisc = 0.0;

  // initialize transition parameters until user changes them
  acrit = 9.0;
  xstrip[1] = 1.0;
  xstrip[2] = 1.0;

  //---- initialize freestream mach number to zero
  mach_type = MachType::CONSTANT;
  minf1 = 0.0;

  //---- drop tolerance for bl system solver
  vaccel = 0.01;
  //---- default viscous parameters
  reynolds_type = ReynoldsType::CONSTANT;
  reinf1 = 0.0;

  initialize();
}

XFoil::~XFoil() {}

/** ---------------------------------------------------
 *      variable initialization/default routine.
 * --------------------------------------------------- */
bool XFoil::initialize() {

  hopi = 0.50 / PI;
  qopi = 0.25 / PI;
  dtor = PI / 180.0;

  n = 0;  // so that current airfoil is not initialized

  memset(aijpiv, 0, sizeof(aijpiv));
  memset(apanel, 0, sizeof(apanel));
  memset(blsav, 0, sizeof(blsav));
  memset(aij, 0, sizeof(aij));
  memset(bij, 0, sizeof(bij));
  memset(cij, 0, sizeof(cij));
  memset(cpi, 0, sizeof(cpi));
  memset(cpv, 0, sizeof(cpv));
  memset(ctau, 0, sizeof(ctau));
  memset(ctq, 0, sizeof(ctq));
  memset(dij, 0, sizeof(dij));
  memset(dqdg, 0, sizeof(dqdg));
  memset(dqdm, 0, sizeof(dqdm));
  memset(dstr, 0, sizeof(dstr));
  memset(dzdg, 0, sizeof(dzdg));
  memset(dzdm, 0, sizeof(dzdm));
  memset(dzdn, 0, sizeof(dzdn));
  memset(iblte, 0, sizeof(iblte));
  memset(ipan, 0, sizeof(ipan));
  memset(isys, 0, sizeof(isys));
  memset(itran, 0, sizeof(itran));
  memset(mass, 0, sizeof(mass));
  memset(nbl, 0, sizeof(nbl));
  nx = VectorXd::Zero(IZX);
  ny = VectorXd::Zero(IZX);
  memset(gamu, 0, sizeof(gamu));
  memset(gam, 0, sizeof(gam));
  memset(gam_a, 0, sizeof(gam_a));
  memset(qf0, 0, sizeof(qf0));
  memset(qf1, 0, sizeof(qf1));
  memset(qf2, 0, sizeof(qf2));
  memset(qf3, 0, sizeof(qf3));
  memset(qinv, 0, sizeof(qinv));
  memset(qinvu, 0, sizeof(qinvu));
  memset(qinv_a, 0, sizeof(qinv_a));
  memset(qvis, 0, sizeof(qvis));
  spline_length.resize(IZX);
  memset(sig, 0, sizeof(sig));
  snew = VectorXd::Zero(4 * IBX);
  memset(sig, 0, sizeof(sig));
  memset(thet, 0, sizeof(thet));
  memset(uedg, 0, sizeof(uedg));
  memset(uinv, 0, sizeof(uinv));
  memset(vti, 0, sizeof(vti));
  points.resize(IZX, 2);
  xbp.resize(IBX, 0);
  dpoints_ds.resize(IZX, 2);
  
  memset(xssi, 0, sizeof(xssi));
  
  ybp.resize(IBX, 0);
  
  memset(wgap, 0, sizeof(wgap));
  memset(va, 0, sizeof(va));
  memset(vb, 0, sizeof(vb));
  memset(vdel, 0, sizeof(vdel));
  memset(vm, 0, sizeof(vm));
  memset(vs1, 0, sizeof(vs1));
  memset(vs2, 0, sizeof(vs2));
  memset(vsrez, 0, sizeof(vsrez));
  memset(vsr, 0, sizeof(vsr));
  memset(vsm, 0, sizeof(vsm));
  memset(vsx, 0, sizeof(vsx));
  memset(vz, 0, sizeof(vz));
  w1 = VectorXd::Zero(6 * IQX);
  w2 = VectorXd::Zero(6 * IQX);
  w3 = VectorXd::Zero(6 * IQX);
  w4 = VectorXd::Zero(6 * IQX);

  memset(ctau, 0, sizeof(ctau));
  memset(ctq, 0, sizeof(ctq));
  memset(dstr, 0, sizeof(dstr));
  memset(thet, 0, sizeof(thet));
  memset(uedg, 0, sizeof(uedg));
  memset(itran, 0, sizeof(itran));
  
  // mdes
  memset(qgamm, 0, sizeof(qgamm));

  //---- default cp/cv (air)
  gamma = 1.4;
  gamm1 = gamma - 1.0;

  //---- set unity freestream speed
  qinf = 1.0;

  psio = 0.0;

  cl = 0.0;
  cm = 0.0;
  cd = 0.0;

  sigte = 0.0;
  gamte = 0.0;

  awake = 0.0;
  avisc = 0.0;

  lgamu = false;
  lvisc = false;
  lwake = false;
  lblini = false;
  lipan = false;
  lqaij = false;
  ladij = false;
  lwdij = false;
  lvconv = false;

  sharp = false;
  lalfa = false;
  
  trforc = false;
  simi = false;
  tran = false;
  turb = false;
  wake = false;
  trfree = false;
  tforce[0] = false;
  tforce[1] = false;
  tforce[2] = false;

  //---- circle plane array size (largest 2  + 1 that will fit array size)
  double ann = log(double((2 * IQX) - 1)) / log(2.0);
  int nn = int(ann + 0.00001);
  int tmp = 1;

  for (int l = 0; l < nn; l++) {
    tmp = 2 * tmp;
  }
  nc1 = tmp + 1;
  //	nc1 = (int)pow(2,nn) + 1;
  if (nc1 > ICX) {
    tmp = 1;
    for (int l = 0; l < nn - 1; l++) {
      tmp = 2 * tmp;
    }
    nc1 = tmp + 1;
    //		nc1 = pow(2,(nn-1)) + 1; //257 instead of ICX in original source
    // code
  }

  //---- default cm reference location
  xcmref = 0.25;
  ycmref = 0.0;

  xoctr[1] = 1.0;
  xoctr[2] = 1.0;
  yoctr[1] = 0.0;
  yoctr[2] = 0.0;
  waklen = 1.0;

  // added techwinder : no wake yet
  nw = 0;

  // added techwinder : no flap yet
  hmom = 0.0;

  // added techwinder : fortran initializes to 0
  imxbl = 0;
  ismxbl = 0;
  ist = 0;

  dwte = 0.0;
  qinfbl = 0.0;
  tkbl = 0.0;
  tkbl_ms = 0.0;
  rstbl = 0.0;
  rstbl_ms = 0.0;
  hstinv = 0.0;
  hstinv_ms = 0.0;
  reybl = 0.0;
  reybl_ms = 0.0;
  reybl_re = 0.0;
  gambl = 0.0;
  gm1bl = 0.0;
  hvrat = 0.0;
  bule = 0.0;
  xiforc = 0.0;
  amcrit = 0.0;

  alfa = 0.0;
  amax = 0.0;
  adeg = 0.0;
  rmxbl = 0.0;
  rmsbl = 0.0;
  rlx = 0.0;
  ante = 0.0;
  cpmn = 0.0;
  clspec = 0.0;
  minf = 0.0;
  reinf = 0.0;
  minf_cl = 0.0;
  reinf_cl = 0.0;

  sle = 0.0;

  chord = 0.0;
  cl_alf = 0.0;
  cl_msq = 0.0;
  cosa = 0.0;
  sina = 0.0;
  tklam = 0.0;
  tkl_msq = 0.0;
  cpstar = 0.0;
  qstar = 0.0;
  cpmni = 0.0;
  cpmnv = 0.0;
  xcpmni = 0.0;
  xcpmnv = 0.0;
  sst = 0.0;
  sst_go = 0.0;
  sst_gp = 0.0;
  dste = 0.0;
  aste = 0.0;

  cfm = 0.0;
  cfm_ms = 0.0;
  cfm_re = 0.0;
  cfm_u1 = 0.0;
  cfm_t1 = 0.0;
  cfm_d1 = 0.0;
  cfm_u2 = 0.0;
  cfm_t2 = 0.0;
  cfm_d2 = 0.0;
  xt = 0.0;
  xt_a1 = 0.0;
  xt_ms = 0.0;
  xt_re = 0.0;
  xt_xf = 0.0;
  xt_x1 = 0.0;
  xt_t1 = 0.0;
  xt_d1 = 0.0;
  xt_u1 = 0.0;
  xt_x2 = 0.0;
  xt_t2 = 0.0;
  xt_d2 = 0.0;
  xt_u2 = 0.0;

  //---- drop tolerance for bl system solver
  vaccel = 0.01;

  //---- set minf, reinf, based on current cl-dependence
  mrcl(1.0, minf_cl, reinf_cl);

  //---- set various compressibility parameters from minf
  comset();

  return true;
}

bool XFoil::abcopy(MatrixX2d copyFrom) {

  if (n != copyFrom.rows() - 1) lblini = false;

  n = copyFrom.rows() - 1;
  for (int i = 1; i <= n; i++) {
    points.row(i).x() = copyFrom.row(i).x();
    points.row(i).y() = copyFrom.row(i).y();
  }

  //---- strip out doubled points
  int r = 1;

  while (r < n) {
    r++;
    if (points.row(r - 1).x() == points.row(r).x() && points.row(r - 1).y() == points.row(r).y()) {
      for (int j = r; j <= n - 1; j++) {
        points.row(j).x() = points.row(j + 1).x();
        points.row(j).y() = points.row(j + 1).y();
      }
      n = n - 1;
    }
  }

  spline_length.segment(1, spline_length.size() - 1) = spline::scalc(points.middleRows(1, points.rows() - 1), n, spline_length.size() - 1);
  dpoints_ds.col(0) = spline::splind(points.col(0), spline_length, n);
  dpoints_ds.col(1) = spline::splind(points.col(1), spline_length, n);
  ncalc(points, spline_length, n, nx.data(), ny.data());
  lefind(sle, points, dpoints_ds, spline_length, n);
  point_le.x() = spline::seval(sle, points.col(0), dpoints_ds.col(0), spline_length, n);
  point_le.y() = spline::seval(sle, points.col(1), dpoints_ds.col(1), spline_length, n);
  point_te.x() = 0.5 * (points.row(1).x() + points.row(n).x());
  point_te.y() = 0.5 * (points.row(1).y() + points.row(n).y());
  chord = (point_le - point_te).norm();
  tecalc();
  apcalc();

  lgamu = false;
  lwake = false;
  lqaij = false;
  ladij = false;
  lwdij = false;
  lipan = false;
  lvconv = false;

  return true;
}

double XFoil::aint(double number) {
  if (number >= 0)
    return (double)(int(number));
  else
    return (double)(-int(-number));
}

bool XFoil::apcalc() {
  double sx, sy;

  //---- set angles of airfoil panels
  for (int i = 1; i <= n - 1; i++) {
    sx = points.row(i + 1).x() - points.row(i).x();
    sy = points.row(i + 1).y() - points.row(i).y();
    if (sx == 0.0 && sy == 0.0)
      apanel[i] = atan2(-ny[i], -nx[i]);
    else
      apanel[i] = atan2(sx, -sy);
  }

  //---- TE panel
  if (sharp)
    apanel[n] = PI;
  else {
    sx = points.row(1).x() - points.row(n).x();
    sy = points.row(1).y() - points.row(n).y();
    apanel[n] = atan2(-sx, sy) + PI;
  }

  return true;
}

/** -------------------------------------------------------------
 *     atan2 function with branch cut checking.
 *
 *     increments position angle of point x,y from some previous
 *     value thold due to a change in position, ensuring that the
 *     position change does not cross the atan2 branch cut
 *     (which is in the -x direction).  for example:
 *
 *       atanc( -1.0 , -1.0 , 0.75*PI )  returns  1.25*PI , whereas
 *       atan2( -1.0 , -1.0 )            returns  -.75*PI .
 *
 *     typically, atanc is used to fill an array of angles:
 *
 *        theta(1) = atan2( y(1) , x(1) )
 *        do i=2, n
 *          theta[i] = atanc( y[i] , x[i] , theta(i-1) )
 *        end do
 *
 *     this will prevent the angle array theta(i) from jumping by
 *     +/- 2 pi when the path x(i),y(i) crosses the negative x axis.
 *
 *     input:
 *       x,y     point position coordinates
 *       thold   position angle of nearby point
 *
 *     output:
 *       atanc   position angle of x,y
 * -------------------------------------------------------------- */
double XFoil::atanc(double y, double x, double thold) {
  double tpi, thnew, dthet, dtcorr;
  tpi = 6.2831853071795864769;

  //---- set new position angle, ignoring branch cut in atan2 function for now

  thnew = atan2(y, x);
  dthet = thnew - thold;

  //---- angle change cannot exceed +/- pi, so get rid of any multiples of 2 pi
  dtcorr = dthet - tpi * int((dthet + sign(PI, dthet)) / tpi);

  //---- set correct new angle
  return thold + dtcorr;
}

/** ----------------------------------------------------------
 *      returns average amplification ax over interval 1..2
 * ----------------------------------------------------------- */
bool XFoil::axset(double hk1, double t1, double rt1, double a1, double hk2,
                  double t2, double rt2, double a2, double acrit, double &ax,
                  double &ax_hk1, double &ax_t1, double &ax_rt1, double &ax_a1,
                  double &ax_hk2, double &ax_t2, double &ax_rt2,
                  double &ax_a2) {
  //
  //==========================
  //---- 2nd-order
  double ax1 = 0.0, ax2 = 0.0, ax1_hk1 = 0.0, ax1_t1 = 0.0, ax1_rt1 = 0.0;
  double ax2_hk2 = 0.0, ax2_t2 = 0.0, ax2_rt2 = 0.0, axsq = 0.0;
  double axa = 0.0, axa_ax1 = 0.0, axa_ax2 = 0.0;
  double exn = 0.0, exn_a1 = 0.0, exn_a2 = 0.0, dax = 0.0, dax_a1 = 0.0,
         dax_a2 = 0.0, dax_t1 = 0.0, dax_t2 = 0.0;
  double f_arg = 0.0;  // ex arg

  dampl(hk1, t1, rt1, ax1, ax1_hk1, ax1_t1, ax1_rt1);
  dampl(hk2, t2, rt2, ax2, ax2_hk2, ax2_t2, ax2_rt2);

  //---- rms-average version (seems a little better on coarse grids)
  axsq = 0.5 * (ax1 * ax1 + ax2 * ax2);
  if (axsq <= 0.0) {
    axa = 0.0;
    axa_ax1 = 0.0;
    axa_ax2 = 0.0;
  } else {
    axa = sqrt(axsq);
    axa_ax1 = 0.5 * ax1 / axa;
    axa_ax2 = 0.5 * ax2 / axa;
  }

  //----- small additional term to ensure  dn/dx > 0  near  n = ncrit
  f_arg = std::min(20.0 * (acrit - 0.5 * (a1 + a2)), 20.0);
  if (f_arg <= 0.0) {
    exn = 1.0;
    exn_a1 = 0.0;
    exn_a2 = 0.0;
  } else {
    exn = exp(-f_arg);
    exn_a1 = 20.0 * 0.5 * exn;
    exn_a2 = 20.0 * 0.5 * exn;
  }

  dax = exn * 0.002 / (t1 + t2);
  dax_a1 = exn_a1 * 0.002 / (t1 + t2);
  dax_a2 = exn_a2 * 0.002 / (t1 + t2);
  dax_t1 = -dax / (t1 + t2);
  dax_t2 = -dax / (t1 + t2);

  //==========================

  ax = axa + dax;
  ax_hk1 = axa_ax1 * ax1_hk1;
  ax_t1 = axa_ax1 * ax1_t1 + dax_t1;
  ax_rt1 = axa_ax1 * ax1_rt1;
  ax_a1 = dax_a1;

  ax_hk2 = axa_ax2 * ax2_hk2;
  ax_t2 = axa_ax2 * ax2_t2 + dax_t2;
  ax_rt2 = axa_ax2 * ax2_rt2;
  ax_a2 = dax_a2;

  return true;
}

/** -----------------------------------------------------------
 *     sets up the newton system coefficients and residuals
 *
 *         ityp = 0 :  similarity station
 *         ityp = 1 :  laminar interval
 *         ityp = 2 :  turbulent interval
 *         ityp = 3 :  wake interval
 *
 *      this routine knows nothing about a transition interval,
 *      which is taken care of by trdif.
 * ------------------------------------------------------------ */
bool XFoil::bldif(int ityp) {

  double hupwt, hdcon, hl, hd_hk1, hd_hk2, hlsq, ehh;
  double upw, upw_hl, upw_hd, upw_hk1, upw_hk2, upw_u1, upw_t1, upw_d1;
  double upw_u2, upw_t2, upw_d2, upw_ms;
  double dxi, slog, scc, scc_usa;
  double ax, ax_t1, ax_rt1, ax_a1, ax_hk2, ax_t2, ax_rt2, ax_a2;
  double rezc, z_ax, hr, hr_hka;
  double hl_hk1, hl_hk2, ax_hk1, sa, cqa, cfa, hka, usa, rta, dea, da, ald;
  double gcc, hkc, hkc_hka, rezt, rezh;
  double btmp, hwa, ha, ma, xa, ta, xlog, ulog, tlog, hlog, ddlog;
  double z_cfx, z_ha, z_hwa, z_ma, z_xl, z_tl, z_cfm, z_t1, z_t2;
  double z_dix, z_hca, z_hl, z_hs1, z_hs2, z_di1, z_di2;
  double z_cfa, z_hka, z_da, z_sl, z_ul, z_dxi, z_usa, z_cqa, z_sa, z_dea;
  double z_upw, z_de1, z_de2, z_us1, z_us2, z_d1, z_d2, z_u1, z_u2, z_x1, z_x2;
  double z_s1, z_s2, z_cq1, z_cq2, z_cf1, z_cf2, z_hk1, z_hk2;
  double cfx, cfx_xa, cfx_ta, cfx_x1, cfx_x2, cfx_t1, cfx_t2, cfx_cf1, cfx_cf2,
      cfx_cfm;
  double xot1, xot2, hca, hsa, dix, dix_upw, cfx_upw;
  double uq, uq_hka, uq_cfa, uq_da;
  double f_arg;  // ex arg
                 //	double scc_us1, scc_us2;

  if (ityp == 0) {
    //----- similarity logarithmic differences  (prescribed)
    xlog = 1.0;
    ulog = bule;
    tlog = 0.5 * (1.0 - bule);
    hlog = 0.0;
    ddlog = 0.0;
  } else {
    //----- usual logarithmic differences
    xlog = log(blData2.xz / blData1.xz);
    ulog = log(blData2.uz / blData1.uz);
    tlog = log(blData2.tz / blData1.tz);
    hlog = log(blData2.hsz / blData1.hsz);
    ddlog = 1.0;
  }

  for (int k = 1; k <= 4; k++) {
    vsrez[k] = 0.0;
    vsm[k] = 0.0;
    vsr[k] = 0.0;
    vsx[k] = 0.0;
    for (int l = 1; l <= 5; l++) {
      vs1[k][l] = 0.0;
      vs2[k][l] = 0.0;
    }
  }

  //---- set triggering constant for local upwinding
  hupwt = 1.0;

  hdcon = 5.0 * hupwt / blData2.hkz / blData2.hkz;
  hd_hk1 = 0.0;
  hd_hk2 = -hdcon * 2.0 / blData2.hkz;

  //---- use less upwinding in the wake
  if (ityp == 3) {
    hdcon = hupwt / blData2.hkz / blData2.hkz;
    hd_hk1 = 0.0;
    hd_hk2 = -hdcon * 2.0 / blData2.hkz;
  }
  //
  //---- local upwinding is based on local change in  log(hk-1)
  //-    (mainly kicks in at transition)
  f_arg = fabs((blData2.hkz - 1.0) / (blData1.hkz - 1.0));
  hl = log(f_arg);
  hl_hk1 = -1.0 / (blData1.hkz - 1.0);
  hl_hk2 = 1.0 / (blData2.hkz - 1.0);

  hlsq = std::min(hl * hl, 15.0);
  ehh = exp(-hlsq * hdcon);
  upw = 1.0 - 0.5 * ehh;
  upw_hl = ehh * hl * hdcon;
  upw_hd = 0.5 * ehh * hlsq;

  upw_hk1 = upw_hl * hl_hk1 + upw_hd * hd_hk1;
  upw_hk2 = upw_hl * hl_hk2 + upw_hd * hd_hk2;

  upw_u1 = upw_hk1 * blData1.hkz_uz;
  upw_t1 = upw_hk1 * blData1.hkz_tz;
  upw_d1 = upw_hk1 * blData1.hkz_dz;
  upw_u2 = upw_hk2 * blData2.hkz_uz;
  upw_t2 = upw_hk2 * blData2.hkz_tz;
  upw_d2 = upw_hk2 * blData2.hkz_dz;
  upw_ms = upw_hk1 * blData1.hkz_ms + upw_hk2 * blData2.hkz_ms;

  if (ityp == 0) {
    //***** le point -->  set zero amplification factor
    vs2[1][1] = 1.0;
    vsr[1] = 0.0;
    vsrez[1] = -blData2.amplz;
  } else {
    if (ityp == 1) {
      //***** laminar part -->  set amplification equation
      //----- set average amplification ax over interval x1..x2

      axset(blData1.hkz, blData1.tz, blData1.rtz, blData1.amplz, blData2.hkz, blData2.tz, blData2.rtz, blData2.amplz, amcrit, ax,
            ax_hk1, ax_t1, ax_rt1, ax_a1, ax_hk2, ax_t2, ax_rt2, ax_a2);

      rezc = blData2.amplz - blData1.amplz - ax * (blData2.xz - blData1.xz);
      z_ax = -(blData2.xz - blData1.xz);

      vs1[1][1] = z_ax * ax_a1 - 1.0;
      vs1[1][2] = z_ax * (ax_hk1 * blData1.hkz_tz + ax_t1 + ax_rt1 * blData1.rtz_tz);
      vs1[1][3] = z_ax * (ax_hk1 * blData1.hkz_dz);
      vs1[1][4] = z_ax * (ax_hk1 * blData1.hkz_uz + ax_rt1 * blData1.rtz_uz);
      vs1[1][5] = ax;
      vs2[1][1] = z_ax * ax_a2 + 1.0;
      vs2[1][2] = z_ax * (ax_hk2 * blData2.hkz_tz + ax_t2 + ax_rt2 * blData2.rtz_tz);
      vs2[1][3] = z_ax * (ax_hk2 * blData2.hkz_dz);
      vs2[1][4] = z_ax * (ax_hk2 * blData2.hkz_uz + ax_rt2 * blData2.rtz_uz);
      vs2[1][5] = -ax;
      vsm[1] = z_ax * (ax_hk1 * blData1.hkz_ms + ax_rt1 * blData1.rtz_ms + ax_hk2 * blData2.hkz_ms +
                       ax_rt2 * blData2.rtz_ms);
      vsr[1] = z_ax * (ax_rt1 * blData1.rtz_re + ax_rt2 * blData2.rtz_re);
      vsx[1] = 0.0;
      vsrez[1] = -rezc;
    } else {
      //***** turbulent part -->  set shear lag equation

      sa = (1.0 - upw) * blData1.sz + upw * blData2.sz;
      cqa = (1.0 - upw) * blData1.cqz + upw * blData2.cqz;
      cfa = (1.0 - upw) * blData1.cfz + upw * blData2.cfz;
      hka = (1.0 - upw) * blData1.hkz + upw * blData2.hkz;

      usa = 0.5 * (blData1.usz + blData2.usz);
      rta = 0.5 * (blData1.rtz + blData2.rtz);
      dea = 0.5 * (blData1.dez + blData2.dez);
      da = 0.5 * (blData1.dz + blData2.dz);

      if (ityp == 3)
        ald = dlcon;  //------ increased dissipation length in wake (decrease
                      // its reciprocal)
      else
        ald = 1.0;

      //----- set and linearize  equilibrium 1/ue due/dx   ...  new  12 oct 94
      if (ityp == 2) {
        gcc = gccon;
        hkc = hka - 1.0 - gcc / rta;
        hkc_hka = 1.0;
        if (hkc < 0.01) {
          hkc = 0.01;
          hkc_hka = 0.0;
        }
      } else {        
        hkc = hka - 1.0;
        hkc_hka = 1.0;
      }

      hr = hkc / (gacon * ald * hka);
      hr_hka = hkc_hka / (gacon * ald * hka) - hr / hka;

      uq = (0.5 * cfa - hr * hr) / (gbcon * da);
      uq_hka = -2.0 * hr * hr_hka / (gbcon * da);
      uq_cfa = 0.5 / (gbcon * da);
      uq_da = -uq / da;

      scc = sccon * 1.333 / (1.0 + usa);
      scc_usa = -scc / (1.0 + usa);

      slog = log(blData2.sz / blData1.sz);
      dxi = blData2.xz - blData1.xz;

      rezc = scc * (cqa - sa * ald) * dxi - dea * 2.0 * slog +
             dea * 2.0 * (uq * dxi - ulog);

      z_cfa = dea * 2.0 * uq_cfa * dxi;
      z_hka = dea * 2.0 * uq_hka * dxi;
      z_da = dea * 2.0 * uq_da * dxi;
      z_sl = -dea * 2.0;
      z_ul = -dea * 2.0;
      z_dxi = scc * (cqa - sa * ald) + dea * 2.0 * uq;
      z_usa = scc_usa * (cqa - sa * ald) * dxi;
      z_cqa = scc * dxi;
      z_sa = -scc * dxi * ald;
      z_dea = 2.0 * (uq * dxi - ulog - slog);

      z_upw = z_cqa * (blData2.cqz - blData1.cqz) + z_sa * (blData2.sz - blData1.sz) + z_cfa * (blData2.cfz - blData1.cfz) +
              z_hka * (blData2.hkz - blData1.hkz);
      z_de1 = 0.5 * z_dea;
      z_de2 = 0.5 * z_dea;
      z_us1 = 0.5 * z_usa;
      z_us2 = 0.5 * z_usa;
      z_d1 = 0.5 * z_da;
      z_d2 = 0.5 * z_da;
      z_u1 = -z_ul / blData1.uz;
      z_u2 = z_ul / blData2.uz;
      z_x1 = -z_dxi;
      z_x2 = z_dxi;
      z_s1 = (1.0 - upw) * z_sa - z_sl / blData1.sz;
      z_s2 = upw * z_sa + z_sl / blData2.sz;
      z_cq1 = (1.0 - upw) * z_cqa;
      z_cq2 = upw * z_cqa;
      z_cf1 = (1.0 - upw) * z_cfa;
      z_cf2 = upw * z_cfa;
      z_hk1 = (1.0 - upw) * z_hka;
      z_hk2 = upw * z_hka;

      vs1[1][1] = z_s1;
      vs1[1][2] = z_upw * upw_t1 + z_de1 * blData1.dez_tz + z_us1 * blData1.usz_tz;
      vs1[1][3] = z_d1 + z_upw * upw_d1 + z_de1 * blData1.dez_dz + z_us1 * blData1.usz_dz;
      vs1[1][4] = z_u1 + z_upw * upw_u1 + z_de1 * blData1.dez_uz + z_us1 * blData1.usz_uz;
      vs1[1][5] = z_x1;
      vs2[1][1] = z_s2;
      vs2[1][2] = z_upw * upw_t2 + z_de2 * blData2.dez_tz + z_us2 * blData2.usz_tz;
      vs2[1][3] = z_d2 + z_upw * upw_d2 + z_de2 * blData2.dez_dz + z_us2 * blData2.usz_dz;
      vs2[1][4] = z_u2 + z_upw * upw_u2 + z_de2 * blData2.dez_uz + z_us2 * blData2.usz_uz;
      vs2[1][5] = z_x2;
      vsm[1] = z_upw * upw_ms + z_de1 * blData1.dez_ms + z_us1 * blData1.usz_ms +
               z_de2 * blData2.dez_ms + z_us2 * blData2.usz_ms;

      vs1[1][2] = vs1[1][2] + z_cq1 * blData1.cqz_tz + z_cf1 * blData1.cfz_tz + z_hk1 * blData1.hkz_tz;
      vs1[1][3] = vs1[1][3] + z_cq1 * blData1.cqz_dz + z_cf1 * blData1.cfz_dz + z_hk1 * blData1.hkz_dz;
      vs1[1][4] = vs1[1][4] + z_cq1 * blData1.cqz_uz + z_cf1 * blData1.cfz_uz + z_hk1 * blData1.hkz_uz;

      vs2[1][2] = vs2[1][2] + z_cq2 * blData2.cqz_tz + z_cf2 * blData2.cfz_tz + z_hk2 * blData2.hkz_tz;
      vs2[1][3] = vs2[1][3] + z_cq2 * blData2.cqz_dz + z_cf2 * blData2.cfz_dz + z_hk2 * blData2.hkz_dz;
      vs2[1][4] = vs2[1][4] + z_cq2 * blData2.cqz_uz + z_cf2 * blData2.cfz_uz + z_hk2 * blData2.hkz_uz;

      vsm[1] = vsm[1] + z_cq1 * blData1.cqz_ms + z_cf1 * blData1.cfz_ms + z_hk1 * blData1.hkz_ms +
               z_cq2 * blData2.cqz_ms + z_cf2 * blData2.cfz_ms + z_hk2 * blData2.hkz_ms;
      vsr[1] =
          z_cq1 * blData1.cqz_re + z_cf1 * blData1.cfz_re + z_cq2 * blData2.cqz_re + z_cf2 * blData2.cfz_re;
      vsx[1] = 0.0;
      vsrez[1] = -rezc;
    }
  }  // endif

  //**** set up momentum equation
  ha = 0.5 * (blData1.hz + blData2.hz);
  ma = 0.5 * (blData1.mz + blData2.mz);
  xa = 0.5 * (blData1.xz + blData2.xz);
  ta = 0.5 * (blData1.tz + blData2.tz);
  hwa = 0.5 * (blData1.dwz / blData1.tz + blData2.dwz / blData2.tz);

  //---- set cf term, using central value cfm for better accuracy in drag
  cfx = 0.50 * cfm * xa / ta + 0.25 * (blData1.cfz * blData1.xz / blData1.tz + blData2.cfz * blData2.xz / blData2.tz);
  cfx_xa = 0.50 * cfm / ta;
  cfx_ta = -.50 * cfm * xa / ta / ta;

  cfx_x1 = 0.25 * blData1.cfz / blData1.tz + cfx_xa * 0.5;
  cfx_x2 = 0.25 * blData2.cfz / blData2.tz + cfx_xa * 0.5;
  cfx_t1 = -.25 * blData1.cfz * blData1.xz / blData1.tz / blData1.tz + cfx_ta * 0.5;
  cfx_t2 = -.25 * blData2.cfz * blData2.xz / blData2.tz / blData2.tz + cfx_ta * 0.5;
  cfx_cf1 = 0.25 * blData1.xz / blData1.tz;
  cfx_cf2 = 0.25 * blData2.xz / blData2.tz;
  cfx_cfm = 0.50 * xa / ta;

  btmp = ha + 2.0 - ma + hwa;

  rezt = tlog + btmp * ulog - xlog * 0.5 * cfx;
  z_cfx = -xlog * 0.5;
  z_ha = ulog;
  z_hwa = ulog;
  z_ma = -ulog;
  z_xl = -ddlog * 0.5 * cfx;
  z_ul = ddlog * btmp;
  z_tl = ddlog;

  z_cfm = z_cfx * cfx_cfm;
  z_cf1 = z_cfx * cfx_cf1;
  z_cf2 = z_cfx * cfx_cf2;

  z_t1 =
      -z_tl / blData1.tz + z_cfx * cfx_t1 + z_hwa * 0.5 * (-blData1.dwz / blData1.tz / blData1.tz);
  z_t2 =
      z_tl / blData2.tz + z_cfx * cfx_t2 + z_hwa * 0.5 * (-blData2.dwz / blData2.tz / blData2.tz);
  z_x1 = -z_xl / blData1.xz + z_cfx * cfx_x1;
  z_x2 = z_xl / blData2.xz + z_cfx * cfx_x2;
  z_u1 = -z_ul / blData1.uz;
  z_u2 = z_ul / blData2.uz;

  vs1[2][2] = 0.5 * z_ha * blData1.hz_tz + z_cfm * cfm_t1 + z_cf1 * blData1.cfz_tz + z_t1;
  vs1[2][3] = 0.5 * z_ha * blData1.hz_dz + z_cfm * cfm_d1 + z_cf1 * blData1.cfz_dz;
  vs1[2][4] = 0.5 * z_ma * blData1.mz_uz + z_cfm * cfm_u1 + z_cf1 * blData1.cfz_uz + z_u1;
  vs1[2][5] = z_x1;
  vs2[2][2] = 0.5 * z_ha * blData2.hz_tz + z_cfm * cfm_t2 + z_cf2 * blData2.cfz_tz + z_t2;
  vs2[2][3] = 0.5 * z_ha * blData2.hz_dz + z_cfm * cfm_d2 + z_cf2 * blData2.cfz_dz;
  vs2[2][4] = 0.5 * z_ma * blData2.mz_uz + z_cfm * cfm_u2 + z_cf2 * blData2.cfz_uz + z_u2;
  vs2[2][5] = z_x2;

  vsm[2] = 0.5 * z_ma * blData1.mz_ms + z_cfm * cfm_ms + z_cf1 * blData1.cfz_ms +
           0.5 * z_ma * blData2.mz_ms + z_cf2 * blData2.cfz_ms;
  vsr[2] = z_cfm * cfm_re + z_cf1 * blData1.cfz_re + z_cf2 * blData2.cfz_re;
  vsx[2] = 0.0;
  vsrez[2] = -rezt;

  //**** set up shape parameter equation

  xot1 = blData1.xz / blData1.tz;
  xot2 = blData2.xz / blData2.tz;

  ha = 0.5 * (blData1.hz + blData2.hz);
  hsa = 0.5 * (blData1.hsz + blData2.hsz);
  hca = 0.5 * (blData1.hcz + blData2.hcz);
  hwa = 0.5 * (blData1.dwz / blData1.tz + blData2.dwz / blData2.tz);

  dix = (1.0 - upw) * blData1.diz * xot1 + upw * blData2.diz * xot2;
  cfx = (1.0 - upw) * blData1.cfz * xot1 + upw * blData2.cfz * xot2;
  dix_upw = blData2.diz * xot2 - blData1.diz * xot1;
  cfx_upw = blData2.cfz * xot2 - blData1.cfz * xot1;

  btmp = 2.0 * hca / hsa + 1.0 - ha - hwa;

  rezh = hlog + btmp * ulog + xlog * (0.5 * cfx - dix);
  z_cfx = xlog * 0.5;
  z_dix = -xlog;
  z_hca = 2.0 * ulog / hsa;
  z_ha = -ulog;
  z_hwa = -ulog;
  z_xl = ddlog * (0.5 * cfx - dix);
  z_ul = ddlog * btmp;
  z_hl = ddlog;

  z_upw = z_cfx * cfx_upw + z_dix * dix_upw;

  z_hs1 = -hca * ulog / hsa / hsa - z_hl / blData1.hsz;
  z_hs2 = -hca * ulog / hsa / hsa + z_hl / blData2.hsz;

  z_cf1 = (1.0 - upw) * z_cfx * xot1;
  z_cf2 = upw * z_cfx * xot2;
  z_di1 = (1.0 - upw) * z_dix * xot1;
  z_di2 = upw * z_dix * xot2;

  z_t1 = (1.0 - upw) * (z_cfx * blData1.cfz + z_dix * blData1.diz) * (-xot1 / blData1.tz);
  z_t2 = upw * (z_cfx * blData2.cfz+ z_dix * blData2.diz) * (-xot2 / blData2.tz);
  z_x1 = (1.0 - upw) * (z_cfx * blData1.cfz + z_dix * blData1.diz) / blData1.tz - z_xl / blData1.xz;
  z_x2 = upw * (z_cfx * blData2.cfz + z_dix * blData2.diz) / blData2.tz + z_xl / blData2.xz;
  z_u1 = -z_ul / blData1.uz;
  z_u2 = z_ul / blData2.uz;

  z_t1 = z_t1 + z_hwa * 0.5 * (-blData1.dwz / blData1.tz / blData1.tz);
  z_t2 = z_t2 + z_hwa * 0.5 * (-blData2.dwz / blData2.tz / blData2.tz);

  vs1[3][1] = z_di1 * blData1.diz_sz;
  vs1[3][2] = z_hs1 * blData1.hsz_tz + z_cf1 * blData1.cfz_tz + z_di1 * blData1.diz_tz + z_t1;
  vs1[3][3] = z_hs1 * blData1.hsz_dz + z_cf1 * blData1.cfz_dz + z_di1 * blData1.diz_dz;
  vs1[3][4] = z_hs1 * blData1.hsz_uz + z_cf1 * blData1.cfz_uz + z_di1 * blData1.diz_uz + z_u1;
  vs1[3][5] = z_x1;
  vs2[3][1] = z_di2 * blData2.diz_sz;
  vs2[3][2] = z_hs2 * blData2.hsz_tz + z_cf2 * blData2.cfz_tz + z_di2 * blData2.diz_tz + z_t2;
  vs2[3][3] = z_hs2 * blData2.hsz_dz + z_cf2 * blData2.cfz_dz + z_di2 * blData2.diz_dz;
  vs2[3][4] = z_hs2 * blData2.hsz_uz + z_cf2 * blData2.cfz_uz + z_di2 * blData2.diz_uz + z_u2;
  vs2[3][5] = z_x2;
  vsm[3] = z_hs1 * blData1.hsz_ms + z_cf1 * blData1.cfz_ms + z_di1 * blData1.diz_ms + z_hs2 * blData2.hsz_ms +
           z_cf2 * blData2.cfz_ms + z_di2 * blData2.diz_ms;
  vsr[3] = z_hs1 * blData1.hsz_re + z_cf1 * blData1.cfz_re + z_di1 * blData1.diz_re + z_hs2 * blData2.hsz_re +
           z_cf2 * blData2.cfz_re + z_di2 * blData2.diz_re;

  vs1[3][2] =
      vs1[3][2] + 0.5 * (z_hca * blData1.hcz_tz + z_ha * blData1.hz_tz) + z_upw * upw_t1;
  vs1[3][3] =
      vs1[3][3] + 0.5 * (z_hca * blData1.hcz_dz + z_ha * blData1.hz_dz) + z_upw * upw_d1;
  vs1[3][4] = vs1[3][4] + 0.5 * (z_hca * blData1.hcz_uz) + z_upw * upw_u1;
  vs2[3][2] =
      vs2[3][2] + 0.5 * (z_hca * blData2.hcz_tz + z_ha * blData2.hz_tz) + z_upw * upw_t2;
  vs2[3][3] =
      vs2[3][3] + 0.5 * (z_hca * blData2.hcz_dz + z_ha * blData2.hz_dz) + z_upw * upw_d2;
  vs2[3][4] = vs2[3][4] + 0.5 * (z_hca * blData2.hcz_uz) + z_upw * upw_u2;

  vsm[3] =
      vsm[3] + 0.5 * (z_hca * blData1.hcz_ms) + z_upw * upw_ms + 0.5 * (z_hca * blData2.hcz_ms);

  vsx[3] = 0.0;
  vsrez[3] = -rezh;

  return true;
}

bool XFoil::blkin() {
  //----------------------------------------------------------
  //     calculates turbulence-independent secondary "2"
  //     variables from the primary "2" variables.
  //----------------------------------------------------------
  double tr2, herat, he_u2, he_ms, v2_he, hk2_h2, hk2_m2;
  //---- set edge mach number ** 2
  blData2.mz = blData2.uz * blData2.uz * hstinv / (gm1bl * (1.0 - 0.5 * blData2.uz * blData2.uz * hstinv));
  tr2 = 1.0 + 0.5 * gm1bl * blData2.mz;
  blData2.mz_uz = 2.0 * blData2.mz * tr2 / blData2.uz;
  blData2.mz_ms = blData2.uz * blData2.uz * tr2 / (gm1bl * (1.0 - 0.5 * blData2.uz * blData2.uz * hstinv)) * hstinv_ms;

  //---- set edge density (isentropic relation)
  blData2.rz = rstbl * pow(tr2, (-1.0 / gm1bl));
  blData2.rz_uz = -blData2.rz / tr2 * 0.5 * blData2.mz_uz;
  blData2.rz_ms = -blData2.rz / tr2 * 0.5 * blData2.mz_ms + rstbl_ms * pow(tr2, (-1.0 / gm1bl));

  //---- set shape parameter
  blData2.hz = blData2.dz / blData2.tz;
  blData2.hz_dz = 1.0 / blData2.tz;
  blData2.hz_tz = -blData2.hz / blData2.tz;

  //---- set edge static/stagnation enthalpy
  herat = 1.0 - 0.5 * blData2.uz * blData2.uz * hstinv;
  he_u2 = -blData2.uz * hstinv;
  he_ms = -0.5 * blData2.uz * blData2.uz * hstinv_ms;
  //---- set molecular viscosity
  blData2.vz = sqrt(herat * herat * herat) * (1.0 + hvrat) / (herat + hvrat) / reybl;
  v2_he = blData2.vz * (1.5 / herat - 1.0 / (herat + hvrat));

  blData2.vz_uz = v2_he * he_u2;
  blData2.vz_ms = -blData2.vz / reybl * reybl_ms + v2_he * he_ms;
  blData2.vz_re = -blData2.vz / reybl * reybl_re;

  //---- set kinematic shape parameter
  hkin(blData2.hz, blData2.mz, blData2.hkz, hk2_h2, hk2_m2);

  blData2.hkz_uz = hk2_m2 * blData2.mz_uz;
  blData2.hkz_tz = hk2_h2 * blData2.hz_tz;
  blData2.hkz_dz = hk2_h2 * blData2.hz_dz;
  blData2.hkz_ms = hk2_m2 * blData2.mz_ms;

  //---- set momentum thickness reynolds number
  blData2.rtz = blData2.rz * blData2.uz * blData2.tz / blData2.vz;
  blData2.rtz_uz = blData2.rtz * (1.0 / blData2.uz + blData2.rz_uz / blData2.rz - blData2.vz_uz / blData2.vz);
  blData2.rtz_tz = blData2.rtz / blData2.tz;
  blData2.rtz_ms = blData2.rtz * (blData2.rz_ms / blData2.rz - blData2.vz_ms / blData2.vz);
  blData2.rtz_re = blData2.rtz * (-blData2.vz_re / blData2.vz);

  return true;
}

bool XFoil::blmid(int ityp) {
  //----------------------------------------------------
  //     calculates midpoint skin friction cfm
  //
  //      ityp = 1 :  laminar
  //      ityp = 2 :  turbulent
  //      ityp = 3 :  turbulent wake
  //----------------------------------------------------
  //
  double hka, rta, ma, cfm_rta, cfm_ma;
  double cfml, cfml_hka, cfml_rta, cfml_ma, cfm_hka;
  //---- set similarity variables if not defined
  if (simi) {
    blData1.hkz = blData2.hkz;
    blData1.hkz_tz = blData2.hkz_tz;
    blData1.hkz_dz = blData2.hkz_dz;
    blData1.hkz_uz = blData2.hkz_uz;
    blData1.hkz_ms = blData2.hkz_ms;
    blData1.rtz = blData2.rtz;
    blData1.rtz_tz = blData2.rtz_tz;
    blData1.rtz_uz = blData2.rtz_uz;
    blData1.rtz_ms = blData2.rtz_ms;
    blData1.rtz_re = blData2.rtz_re;
    blData1.mz = blData2.mz;
    blData1.mz_uz = blData2.mz_uz;
    blData1.mz_ms = blData2.mz_ms;
  }

  //---- define stuff for midpoint cf
  hka = 0.5 * (blData1.hkz + blData2.hkz);
  rta = 0.5 * (blData1.rtz + blData2.rtz);
  ma = 0.5 * (blData1.mz + blData2.mz);

  //---- midpoint skin friction coefficient  (zero in wake)
  if (ityp == 3) {
    cfm = 0.0;
    cfm_hka = 0.0;
    cfm_rta = 0.0;
    cfm_ma = 0.0;
    cfm_ms = 0.0;
  } else {
     
    if (ityp == 1) {
      C_f c_f = cfl(hka, rta);
      cfm = c_f.cf;
      cfm_hka = c_f.hk;
      cfm_rta = c_f.rt;
      cfm_ma = c_f.msq;
    } else {
      C_f c_fl = cfl(hka, rta);
      cfml = c_fl.cf;
      cfml_hka = c_fl.hk;
      cfml_rta = c_fl.rt;
      cfml_ma = c_fl.msq;
      C_f c_ft = cft(hka, rta, ma);
      cfm = c_ft.cf;
      cfm_hka = c_ft.hk;
      cfm_rta = c_ft.rt;
      cfm_ma = c_ft.msq;
      if (cfml > cfm) {
        cfm = cfml;
        cfm_hka = cfml_hka;
        cfm_rta = cfml_rta;
        cfm_ma = cfml_ma;
      }
    }
  }
  cfm_u1 = 0.5 * (cfm_hka * blData1.hkz_uz + cfm_ma * blData1.mz_uz + cfm_rta * blData1.rtz_uz);
  cfm_t1 = 0.5 * (cfm_hka * blData1.hkz_tz + cfm_rta * blData1.rtz_tz);
  cfm_d1 = 0.5 * (cfm_hka * blData1.hkz_dz);

  cfm_u2 = 0.5 * (cfm_hka * blData2.hkz_uz + cfm_ma * blData2.mz_uz + cfm_rta * blData2.rtz_uz);
  cfm_t2 = 0.5 * (cfm_hka * blData2.hkz_tz + cfm_rta * blData2.rtz_tz);
  cfm_d2 = 0.5 * (cfm_hka * blData2.hkz_dz);

  cfm_ms = 0.5 * (cfm_hka * blData1.hkz_ms + cfm_ma * blData1.mz_ms + cfm_rta * blData1.rtz_ms +
                  cfm_hka * blData2.hkz_ms + cfm_ma * blData2.mz_ms + cfm_rta * blData2.rtz_ms);
  cfm_re = 0.5 * (cfm_rta * blData1.rtz_re + cfm_rta * blData2.rtz_re);

  return true;
}

/** ----------------------------------------------------------
 *     set bl primary "2" variables from parameter list
 *  ---------------------------------------------------------- */
bool XFoil::blprv(double xsi, double ami, double cti, double thi, double dsi,
                  double dswaki, double uei) {
  blData2.xz = xsi;
  blData2.amplz = ami;
  blData2.sz = cti;
  blData2.tz = thi;
  blData2.dz = dsi - dswaki;
  blData2.dwz = dswaki;

  blData2.uz = uei * (1.0 - tkbl) / (1.0 - tkbl * (uei / qinfbl) * (uei / qinfbl));
  blData2.uz_uei = (1.0 + tkbl * (2.0 * blData2.uz * uei / qinfbl / qinfbl - 1.0)) /
           (1.0 - tkbl * (uei / qinfbl) * (uei / qinfbl));
  blData2.uz_ms = (blData2.uz * (uei / qinfbl) * (uei / qinfbl) - uei) * tkbl_ms /
          (1.0 - tkbl * (uei / qinfbl) * (uei / qinfbl));
  return true;
}

/** -----------------------------------------------------------------
 *      custom solver for coupled viscous-inviscid newton system:
 *
 *        a  |  |  .  |  |  .  |    d       r       s
 *        b  a  |  .  |  |  .  |    d       r       s
 *        |  b  a  .  |  |  .  |    d       r       s
 *        .  .  .  .  |  |  .  |    d   =   r - dre s
 *        |  |  |  b  a  |  .  |    d       r       s
 *        |  z  |  |  b  a  .  |    d       r       s
 *        .  .  .  .  .  .  .  |    d       r       s
 *        |  |  |  |  |  |  b  a    d       r       s
 *
 *       a, b, z  3x3  blocks containing linearized bl equation coefficients
 *       |        3x1  vectors containing mass defect influence
 *                     coefficients on ue
 *       d        3x1  unknown vectors (newton deltas for ctau][ theta][ m)
 *       r        3x1  residual vectors
 *       s        3x1  re influence vectors
 * ------------------------------------------------------------------ */
bool XFoil::blsolve() {
  
  double vtmp, vtmp3;

  int ivte1 = isys[iblte[1]][1];
  //
  for (int iv = 1; iv <= nsys; iv++) {
    //
    int ivp = iv + 1;
    //
    //====== invert va[iv] block
    //
    //------ normalize first row
    double pivot = 1.0 / va[1][1][iv];
    va[1][2][iv] *= pivot;
    for (int l = iv; l <= nsys; l++) vm[1][l][iv] *= pivot;
    vdel[1][1][iv] *= pivot;
    vdel[1][2][iv] *= pivot;
    //
    //------ eliminate lower first column in va block
    for (int k = 2; k <= 3; k++) {
      vtmp = va[k][1][iv];
      va[k][2][iv] -= vtmp * va[1][2][iv];
      for (int l = iv; l <= nsys; l++) vm[k][l][iv] -= vtmp * vm[1][l][iv];
      vdel[k][1][iv] -= vtmp * vdel[1][1][iv];
      vdel[k][2][iv] -= vtmp * vdel[1][2][iv];
    }
    //
    //------ normalize second row
    pivot = 1.0 / va[2][2][iv];
    for (int l = iv; l <= nsys; l++) vm[2][l][iv] *= pivot;
    vdel[2][1][iv] *= pivot;
    vdel[2][2][iv] *= pivot;
    //
    //------ eliminate lower second column in va block
    
    vtmp = va[3][2][iv];
    for (int l = iv; l <= nsys; l++) vm[3][l][iv] -= vtmp * vm[2][l][iv];
    vdel[3][1][iv] -= vtmp * vdel[2][1][iv];
    vdel[3][2][iv] -= vtmp * vdel[2][2][iv];

    //------ normalize third row
    pivot = 1.0 / vm[3][iv][iv];
    for (int l = ivp; l <= nsys; l++) vm[3][l][iv] *= pivot;
    vdel[3][1][iv] *= pivot;
    vdel[3][2][iv] *= pivot;
    //
    //
    //------ eliminate upper third column in va block
    double vtmp1 = vm[1][iv][iv];
    double vtmp2 = vm[2][iv][iv];
    for (int l = ivp; l <= nsys; l++) {
      vm[1][l][iv] -= vtmp1 * vm[3][l][iv];
      vm[2][l][iv] -= vtmp2 * vm[3][l][iv];
    }
    vdel[1][1][iv] -= vtmp1 * vdel[3][1][iv];
    vdel[2][1][iv] -= vtmp2 * vdel[3][1][iv];
    vdel[1][2][iv] -= vtmp1 * vdel[3][2][iv];
    vdel[2][2][iv] -= vtmp2 * vdel[3][2][iv];
    //
    //------ eliminate upper second column in va block
    vtmp = va[1][2][iv];
    for (int l = ivp; l <= nsys; l++) vm[1][l][iv] -= vtmp * vm[2][l][iv];

    vdel[1][1][iv] -= vtmp * vdel[2][1][iv];
    vdel[1][2][iv] -= vtmp * vdel[2][2][iv];
    //
    //
    if (iv != nsys) {
      //
      //====== eliminate vb(iv+1) block][ rows  1 -> 3
      for (int k = 1; k <= 3; k++) {
        vtmp1 = vb[k][1][ivp];
        vtmp2 = vb[k][2][ivp];
        vtmp3 = vm[k][iv][ivp];
        for (int l = ivp; l <= nsys; l++)
          vm[k][l][ivp] -= (vtmp1 * vm[1][l][iv] + vtmp2 * vm[2][l][iv] +
                            vtmp3 * vm[3][l][iv]);
        vdel[k][1][ivp] -= (vtmp1 * vdel[1][1][iv] + vtmp2 * vdel[2][1][iv] +
                            vtmp3 * vdel[3][1][iv]);
        vdel[k][2][ivp] -= (vtmp1 * vdel[1][2][iv] + vtmp2 * vdel[2][2][iv] +
                            vtmp3 * vdel[3][2][iv]);
      }
      //
      if (iv == ivte1) {
        //------- eliminate vz block
        int ivz = isys[iblte[2] + 1][2];
        //
        for (int k = 1; k <= 3; k++) {
          vtmp1 = vz[k][1];
          vtmp2 = vz[k][2];
          for (int l = ivp; l <= nsys; l++) {
            vm[k][l][ivz] -= (vtmp1 * vm[1][l][iv] + vtmp2 * vm[2][l][iv]);
          }
          vdel[k][1][ivz] -= (vtmp1 * vdel[1][1][iv] + vtmp2 * vdel[2][1][iv]);
          vdel[k][2][ivz] -= (vtmp1 * vdel[1][2][iv] + vtmp2 * vdel[2][2][iv]);
        }
      }
      //
      if (ivp != nsys) {
        //
        //====== eliminate lower vm column
        for (int kv = iv + 2; kv <= nsys; kv++) {
          vtmp1 = vm[1][iv][kv];
          vtmp2 = vm[2][iv][kv];
          vtmp3 = vm[3][iv][kv];
          //
          if (fabs(vtmp1) > vaccel) {
            for (int l = ivp; l <= nsys; l++) vm[1][l][kv] -= vtmp1 * vm[3][l][iv];
            vdel[1][1][kv] -= vtmp1 * vdel[3][1][iv];
            vdel[1][2][kv] -= vtmp1 * vdel[3][2][iv];
          }
          //
          if (fabs(vtmp2) > vaccel) {
            for (int l = ivp; l <= nsys; l++) vm[2][l][kv] -= vtmp2 * vm[3][l][iv];
            vdel[2][1][kv] -= vtmp2 * vdel[3][1][iv];
            vdel[2][2][kv] -= vtmp2 * vdel[3][2][iv];
          }
          //
          if (fabs(vtmp3) > vaccel) {
            for (int l = ivp; l <= nsys; l++) vm[3][l][kv] -= vtmp3 * vm[3][l][iv];
            vdel[3][1][kv] -= vtmp3 * vdel[3][1][iv];
            vdel[3][2][kv] -= vtmp3 * vdel[3][2][iv];
          }
          //
        }
      }
    }
  }  // 1000

  //
  for (int iv = nsys; iv >= 2; iv--) {
    //------ eliminate upper vm columns
    vtmp = vdel[3][1][iv];
    for (int kv = iv - 1; kv >= 1; kv--) {
      vdel[1][1][kv] -= vm[1][iv][kv] * vtmp;
      vdel[2][1][kv] -= vm[2][iv][kv] * vtmp;
      vdel[3][1][kv] -= vm[3][iv][kv] * vtmp;
    }
    vtmp = vdel[3][2][iv];
    for (int kv = iv - 1; kv >= 1; kv--) {
      vdel[1][2][kv] -= vm[1][iv][kv] * vtmp;
      vdel[2][2][kv] -= vm[2][iv][kv] * vtmp;
      vdel[3][2][kv] -= vm[3][iv][kv] * vtmp;
    }
    //
  }
  return true;
}

/** ------------------------------------------------------------------
 *
 *      sets up the bl newton system governing the current interval:
 *
 *      |       ||da1|     |       ||da2|       |     |
 *      |  vs1  ||dt1|  +  |  vs2  ||dt2|   =   |vsrez|
 *      |       ||dd1|     |       ||dd2|       |     |
 *               |du1|              |du2|
 *               |dx1|              |dx2|
 *
 *         3x5    5x1         3x5    5x1          3x1
 *
 *      the system as shown corresponds to a laminar station
 *      if tran, then  ds2  replaces  da2
 *      if turb, then  ds1, ds2  replace  da1, da2
 *
 * ------------------------------------------------------------------ */
bool XFoil::blsys() { 

  //---- calculate secondary bl variables and their sensitivities
  if (wake) {
    blvar(3);
    blmid(3);
  } else {
    if (turb || tran) {
      blvar(2);
      blmid(2);
    } else {
      blvar(1);
      blmid(1);
    }
  }

  //---- for the similarity station, "1" and "2" variables are the same
  if (simi) {
    //		for(int icom=1;icom<= ncom;icom++) com1[icom] = com2[icom];
    stepbl();
  }

  //---- set up appropriate finite difference system for current interval
  if (tran)
    trdif();
  else if (simi)
    bldif(0);
  else if (!turb)
    bldif(1);
  else if (wake)
    bldif(3);
  else
    bldif(2);

  if (simi) {
    //----- at similarity station, "1" variables are really "2" variables
    for (int k = 1; k <= 4; k++) {
      for (int l = 1; l <= 5; l++) {
        vs2[k][l] = vs1[k][l] + vs2[k][l];
        vs1[k][l] = 0.0;
      }
    }
  }

  //---- change system over into incompressible uei and mach
  for (int k = 1; k <= 4; k++) {
    //------ residual derivatives wrt compressible uec
    double res_u1 = vs1[k][4];
    double res_u2 = vs2[k][4];
    double res_ms = vsm[k];

    //------ combine with derivatives of compressible  u1,u2 = uec(uei m)
    vs1[k][4] = res_u1 * blData1.uz_uei;
    vs2[k][4] = res_u2 * blData2.uz_uei;
    vsm[k] = res_u1 * blData1.uz_ms + res_u2 * blData2.uz_ms + res_ms;
  }
  return true;
}

/** ----------------------------------------------------
 *      calculates all secondary "2" variables from
 *      the primary "2" variables x2, u2, t2, d2, s2.
 *      also calculates the sensitivities of the
 *      secondary variables wrt the primary variables.
 *
 *       ityp = 1 :  laminar
 *       ityp = 2 :  turbulent
 *       ityp = 3 :  turbulent wake
 * ---------------------------------------------------- */
bool XFoil::blvar(int ityp) {
  double hs2_hk2, hs2_rt2, hs2_m2;
  double hc2_hk2, hc2_m2, us2_hs2, us2_hk2, us2_h2;
  double hkc, hkc_hk2, hkc_rt2, hkb, usb;
  double cq2_hs2, cq2_us2, cq2_hk2;
  double cq2_rt2, cq2_h2, cf2_hk2, cf2_rt2, cf2_m2;
  double cf2l, cf2l_hk2, cf2l_rt2, cf2l_m2;
  double cf2t, cf2t_hk2;
  double cf2t_rt2, cf2t_m2, cf2t_u2, cf2t_t2;
  double cf2t_d2, cf2t_ms, cf2t_re, di2_hs2;
  double di2_us2, di2_cf2t, hmin, hm_rt2;
  double grt, fl, fl_hk2, fl_rt2, tfl;
  double dfac, df_fl, df_hk2, df_rt2;
  double di2l, di2l_hk2, di2l_rt2, de2_hk2, hdmax;

  //	double gbcon, gccon, ctcon, hkc2;//were are they initialized?
  if (ityp == 3) blData2.hkz = std::max(blData2.hkz, 1.00005);
  if (ityp != 3) blData2.hkz = std::max(blData2.hkz, 1.05000);

  //---- density thickness shape parameter     ( h** )
  hct(blData2.hkz, blData2.mz, blData2.hcz, hc2_hk2, hc2_m2);
  blData2.hcz_uz = hc2_hk2 * blData2.hkz_uz + hc2_m2 * blData2.mz_uz;
  blData2.hcz_tz = hc2_hk2 * blData2.hkz_tz;
  blData2.hcz_dz = hc2_hk2 * blData2.hkz_dz;
  blData2.hcz_ms = hc2_hk2 * blData2.hkz_ms + hc2_m2 * blData2.mz_ms;

  //---- set ke thickness shape parameter from  h - h*  correlations
  if (ityp == 1)
    hsl(blData2.hkz, blData2.hsz, hs2_hk2, hs2_rt2, hs2_m2);
  else
    hst(blData2.hkz, blData2.rtz, blData2.mz, blData2.hsz, hs2_hk2, hs2_rt2, hs2_m2);

  blData2.hsz_uz = hs2_hk2 * blData2.hkz_uz + hs2_rt2 * blData2.rtz_uz + hs2_m2 * blData2.mz_uz;
  blData2.hsz_tz = hs2_hk2 * blData2.hkz_tz + hs2_rt2 * blData2.rtz_tz;
  blData2.hsz_dz = hs2_hk2 * blData2.hkz_dz;
  blData2.hsz_ms = hs2_hk2 * blData2.hkz_ms + hs2_rt2 * blData2.rtz_ms + hs2_m2 * blData2.mz_ms;
  blData2.hsz_re = hs2_rt2 * blData2.rtz_re;

  //---- normalized slip velocity  us
  blData2.usz = 0.5 * blData2.hsz * (1.0 - (blData2.hkz - 1.0) / (gbcon * blData2.hz));
  us2_hs2 = 0.5 * (1.0 - (blData2.hkz - 1.0) / (gbcon * blData2.hz));
  us2_hk2 = 0.5 * blData2.hsz * (-1.0 / (gbcon * blData2.hz));
  us2_h2 = 0.5 * blData2.hsz * (blData2.hkz - 1.0) / (gbcon * blData2.hz * blData2.hz);

  blData2.usz_uz = us2_hs2 * blData2.hsz_uz + us2_hk2 * blData2.hkz_uz;
  blData2.usz_tz = us2_hs2 * blData2.hsz_tz + us2_hk2 * blData2.hkz_tz + us2_h2 * blData2.hz_tz;
  blData2.usz_dz = us2_hs2 * blData2.hsz_dz + us2_hk2 * blData2.hkz_dz + us2_h2 * blData2.hz_dz;
  blData2.usz_ms = us2_hs2 * blData2.hsz_ms + us2_hk2 * blData2.hkz_ms;
  blData2.usz_re = us2_hs2 * blData2.hsz_re;

  if (ityp <= 2 && blData2.usz > 0.95) {
    //       write(*,*) 'blvar: us clamped:', us2
    blData2.usz = 0.98;
    blData2.usz_uz = 0.0;
    blData2.usz_tz = 0.0;
    blData2.usz_dz = 0.0;
    blData2.usz_ms = 0.0;
    blData2.usz_re = 0.0;
  }

  if (ityp == 3 && blData2.usz > 0.99995) {
    //       write(*,*) 'blvar: wake us clamped:', us2
    blData2.usz = 0.99995;
    blData2.usz_uz = 0.0;
    blData2.usz_tz = 0.0;
    blData2.usz_dz = 0.0;
    blData2.usz_ms = 0.0;
    blData2.usz_re = 0.0;
  }

  //---- equilibrium wake layer shear coefficient (ctau)eq ** 1/2
  //   ...  new  12 oct 94
  hkc = blData2.hkz - 1.0;
  hkc_hk2 = 1.0;
  hkc_rt2 = 0.0;
  if (ityp == 2) {
    const double gcc = gccon;
    hkc = blData2.hkz - 1.0 - gcc / blData2.rtz;
    hkc_hk2 = 1.0;
    hkc_rt2 = gcc / blData2.rtz / blData2.rtz;
    if (hkc < 0.01) {
      hkc = 0.01;
      hkc_hk2 = 0.0;
      hkc_rt2 = 0.0;
    }
  }

  hkb = blData2.hkz - 1.0;
  usb = 1.0 - blData2.usz;
  blData2.cqz = sqrt(ctcon * blData2.hsz * hkb * hkc * hkc / (usb * blData2.hz * blData2.hkz * blData2.hkz));
  cq2_hs2 = ctcon * hkb * hkc * hkc / (usb * blData2.hz * blData2.hkz * blData2.hkz) * 0.5 / blData2.cqz;
  cq2_us2 =
      ctcon * blData2.hsz * hkb * hkc * hkc / (usb * blData2.hz * blData2.hkz * blData2.hkz) / usb * 0.5 / blData2.cqz;
  cq2_hk2 = ctcon * blData2.hsz * hkc * hkc / (usb * blData2.hz * blData2.hkz * blData2.hkz) * 0.5 / blData2.cqz -
            ctcon * blData2.hsz * hkb * hkc * hkc / (usb * blData2.hz * blData2.hkz * blData2.hkz * blData2.hkz) * 2.0 *
                0.5 / blData2.cqz +
            ctcon * blData2.hsz * hkb * hkc / (usb * blData2.hz * blData2.hkz * blData2.hkz) * 2.0 * 0.5 / blData2.cqz *
                hkc_hk2;
  cq2_rt2 = ctcon * blData2.hsz * hkb * hkc / (usb * blData2.hz * blData2.hkz * blData2.hkz) * 2.0 * 0.5 / blData2.cqz *
            hkc_rt2;
  cq2_h2 =
      -ctcon * blData2.hsz * hkb * hkc * hkc / (usb * blData2.hz * blData2.hkz * blData2.hkz) / blData2.hz * 0.5 / blData2.cqz;

  blData2.cqz_uz = cq2_hs2 * blData2.hsz_uz + cq2_us2 * blData2.usz_uz + cq2_hk2 * blData2.hkz_uz;
  blData2.cqz_tz = cq2_hs2 * blData2.hsz_tz + cq2_us2 * blData2.usz_tz + cq2_hk2 * blData2.hkz_tz;
  blData2.cqz_dz = cq2_hs2 * blData2.hsz_dz + cq2_us2 * blData2.usz_dz + cq2_hk2 * blData2.hkz_dz;
  blData2.cqz_ms = cq2_hs2 * blData2.hsz_ms + cq2_us2 * blData2.usz_ms + cq2_hk2 * blData2.hkz_ms;
  blData2.cqz_re = cq2_hs2 * blData2.hsz_re + cq2_us2 * blData2.usz_re;

  blData2.cqz_uz = blData2.cqz_uz + cq2_rt2 * blData2.rtz_uz;
  blData2.cqz_tz = blData2.cqz_tz + cq2_h2 * blData2.hz_tz + cq2_rt2 * blData2.rtz_tz;
  blData2.cqz_dz = blData2.cqz_dz + cq2_h2 * blData2.hz_dz;
  blData2.cqz_ms = blData2.cqz_ms + cq2_rt2 * blData2.rtz_ms;
  blData2.cqz_re = blData2.cqz_re + cq2_rt2 * blData2.rtz_re;

  //---- set skin friction coefficient
  if (ityp == 3) {
    //----- wake
    blData2.cfz = 0.0;
    cf2_hk2 = 0.0;
    cf2_rt2 = 0.0;
    cf2_m2 = 0.0;
  } else {
    if (ityp == 1) {
      //----- laminar
      C_f c_f = cfl(blData2.hkz, blData2.rtz);
      blData2.cfz = c_f.cf;
      cf2_hk2 = c_f.hk;
      cf2_rt2 = c_f.rt;
      cf2_m2 = c_f.msq;
    }
    else {
      //----- turbulent
      C_f c_fl = cfl(blData2.hkz, blData2.rtz);
      cf2l = c_fl.cf;
      cf2l_hk2 = c_fl.hk;
      cf2l_rt2 = c_fl.rt;
      cf2l_m2 = c_fl.msq;
      C_f c_ft = cft(blData2.hkz, blData2.rtz, blData2.mz);
      blData2.cfz = c_ft.cf;
      cf2_hk2 = c_ft.hk;
      cf2_rt2 = c_ft.rt;
      cf2_m2 = c_ft.msq;
      if (cf2l > blData2.cfz) {
        //------- laminar cf is greater than turbulent cf -- use laminar
        //-       (this will only occur for unreasonably small rtheta)
        blData2.cfz = cf2l;
        cf2_hk2 = cf2l_hk2;
        cf2_rt2 = cf2l_rt2;
        cf2_m2 = cf2l_m2;
      }
    }
  }

  blData2.cfz_uz = cf2_hk2 * blData2.hkz_uz + cf2_rt2 * blData2.rtz_uz + cf2_m2 * blData2.mz_uz;
  blData2.cfz_tz = cf2_hk2 * blData2.hkz_tz + cf2_rt2 * blData2.rtz_tz;
  blData2.cfz_dz = cf2_hk2 * blData2.hkz_dz;
  blData2.cfz_ms = cf2_hk2 * blData2.hkz_ms + cf2_rt2 * blData2.rtz_ms + cf2_m2 * blData2.mz_ms;
  blData2.cfz_re = cf2_rt2 * blData2.rtz_re;

  //---- dissipation function    2 cd / h*
  if (ityp == 1) {
    //----- laminar
    double di2_hk2, di2_rt2;
    dil(blData2.hkz, blData2.rtz, blData2.diz, di2_hk2, di2_rt2);

    blData2.diz_uz = di2_hk2 * blData2.hkz_uz + di2_rt2 * blData2.rtz_uz;
    blData2.diz_tz = di2_hk2 * blData2.hkz_tz + di2_rt2 * blData2.rtz_tz;
    blData2.diz_dz = di2_hk2 * blData2.hkz_dz;
    blData2.diz_sz = 0.0;
    blData2.diz_ms = di2_hk2 * blData2.hkz_ms + di2_rt2 * blData2.rtz_ms;
    blData2.diz_re = di2_rt2 * blData2.rtz_re;
  } else {
    if (ityp == 2) {
      //----- turbulent wall contribution
      C_f c_ft = cft(blData2.hkz, blData2.rtz, blData2.mz);
      cf2t = c_ft.cf;
      cf2t_hk2 = c_ft.hk;
      cf2t_rt2 = c_ft.rt;
      cf2t_m2 = c_ft.msq;

      cf2t_u2 = cf2t_hk2 * blData2.hkz_uz + cf2t_rt2 * blData2.rtz_uz + cf2t_m2 * blData2.mz_uz;
      cf2t_t2 = cf2t_hk2 * blData2.hkz_tz + cf2t_rt2 * blData2.rtz_tz;
      cf2t_d2 = cf2t_hk2 * blData2.hkz_dz;
      cf2t_ms = cf2t_hk2 * blData2.hkz_ms + cf2t_rt2 * blData2.rtz_ms + cf2t_m2 * blData2.mz_ms;
      cf2t_re = cf2t_rt2 * blData2.rtz_re;

      blData2.diz = (0.5 * cf2t * blData2.usz) * 2.0 / blData2.hsz;
      di2_hs2 = -(0.5 * cf2t * blData2.usz) * 2.0 / blData2.hsz / blData2.hsz;
      di2_us2 = (0.5 * cf2t) * 2.0 / blData2.hsz;
      di2_cf2t = (0.5 * blData2.usz) * 2.0 / blData2.hsz;

      blData2.diz_sz = 0.0;
      blData2.diz_uz = di2_hs2 * blData2.hsz_uz + di2_us2 * blData2.usz_uz + di2_cf2t * cf2t_u2;
      blData2.diz_tz = di2_hs2 * blData2.hsz_tz + di2_us2 * blData2.usz_tz + di2_cf2t * cf2t_t2;
      blData2.diz_dz = di2_hs2 * blData2.hsz_dz + di2_us2 * blData2.usz_dz + di2_cf2t * cf2t_d2;
      blData2.diz_ms = di2_hs2 * blData2.hsz_ms + di2_us2 * blData2.usz_ms + di2_cf2t * cf2t_ms;
      blData2.diz_re = di2_hs2 * blData2.hsz_re + di2_us2 * blData2.usz_re + di2_cf2t * cf2t_re;

      //----- set minimum hk for wake layer to still exist
      grt = log(blData2.rtz);
      hmin = 1.0 + 2.1 / grt;
      hm_rt2 = -(2.1 / grt / grt) / blData2.rtz;

      //----- set factor dfac for correcting wall dissipation for very low hk
      fl = (blData2.hkz - 1.0) / (hmin - 1.0);
      fl_hk2 = 1.0 / (hmin - 1.0);
      fl_rt2 = (-fl / (hmin - 1.0)) * hm_rt2;

      tfl = tanh(fl);
      dfac = 0.5 + 0.5 * tfl;
      df_fl = 0.5 * (1.0 - tfl * tfl);

      df_hk2 = df_fl * fl_hk2;
      df_rt2 = df_fl * fl_rt2;

      blData2.diz_sz = blData2.diz_sz * dfac;
      blData2.diz_uz = blData2.diz_uz * dfac + blData2.diz * (df_hk2 * blData2.hkz_uz + df_rt2 * blData2.rtz_uz);
      blData2.diz_tz = blData2.diz_tz * dfac + blData2.diz * (df_hk2 * blData2.hkz_tz + df_rt2 * blData2.rtz_tz);
      blData2.diz_dz = blData2.diz_dz * dfac + blData2.diz * (df_hk2 * blData2.hkz_dz);
      blData2.diz_ms = blData2.diz_ms * dfac + blData2.diz * (df_hk2 * blData2.hkz_ms + df_rt2 * blData2.rtz_ms);
      blData2.diz_re = blData2.diz_re * dfac + blData2.diz * (df_rt2 * blData2.rtz_re);
      blData2.diz = blData2.diz * dfac;
    } else {
      //----- zero wall contribution for wake
      blData2.diz = 0.0;
      blData2.diz_sz = 0.0;
      blData2.diz_uz = 0.0;
      blData2.diz_tz = 0.0;
      blData2.diz_dz = 0.0;
      blData2.diz_ms = 0.0;
      blData2.diz_re = 0.0;
    }
  }
  //---- add on turbulent outer layer contribution
  if (ityp != 1) {
    double dd = blData2.sz * blData2.sz * (0.995 - blData2.usz) * 2.0 / blData2.hsz;
    double dd_hs2 = -blData2.sz * blData2.sz * (0.995 - blData2.usz) * 2.0 / blData2.hsz / blData2.hsz;
    double dd_us2 = -blData2.sz * blData2.sz * 2.0 / blData2.hsz;
    const double dd_s2 = blData2.sz * 2.0 * (0.995 - blData2.usz) * 2.0 / blData2.hsz;

    blData2.diz = blData2.diz + dd;
    blData2.diz_sz = dd_s2;
    blData2.diz_uz = blData2.diz_uz + dd_hs2 * blData2.hsz_uz + dd_us2 * blData2.usz_uz;
    blData2.diz_tz = blData2.diz_tz + dd_hs2 * blData2.hsz_tz + dd_us2 * blData2.usz_tz;
    blData2.diz_dz = blData2.diz_dz + dd_hs2 * blData2.hsz_dz + dd_us2 * blData2.usz_dz;
    blData2.diz_ms = blData2.diz_ms + dd_hs2 * blData2.hsz_ms + dd_us2 * blData2.usz_ms;
    blData2.diz_re = blData2.diz_re + dd_hs2 * blData2.hsz_re + dd_us2 * blData2.usz_re;

    //----- add laminar stress contribution to outer layer cd
    dd = 0.15 * (0.995 - blData2.usz) * (0.995 - blData2.usz) / blData2.rtz * 2.0 / blData2.hsz;
    dd_us2 = -0.15 * (0.995 - blData2.usz) * 2.0 / blData2.rtz * 2.0 / blData2.hsz;
    dd_hs2 = -dd / blData2.hsz;
    const double dd_rt2 = -dd / blData2.rtz;

    blData2.diz = blData2.diz + dd;
    blData2.diz_uz = blData2.diz_uz + dd_hs2 * blData2.hsz_uz + dd_us2 * blData2.usz_uz + dd_rt2 * blData2.rtz_uz;
    blData2.diz_tz = blData2.diz_tz + dd_hs2 * blData2.hsz_tz + dd_us2 * blData2.usz_tz + dd_rt2 * blData2.rtz_tz;
    blData2.diz_dz = blData2.diz_dz + dd_hs2 * blData2.hsz_dz + dd_us2 * blData2.usz_dz;
    blData2.diz_ms = blData2.diz_ms + dd_hs2 * blData2.hsz_ms + dd_us2 * blData2.usz_ms + dd_rt2 * blData2.rtz_ms;
    blData2.diz_re = blData2.diz_re + dd_hs2 * blData2.hsz_re + dd_us2 * blData2.usz_re + dd_rt2 * blData2.rtz_re;
  }

  if (ityp == 2) {
    dil(blData2.hkz, blData2.rtz, di2l, di2l_hk2, di2l_rt2);

    if (di2l > blData2.diz) {
      //------- laminar cd is greater than turbulent cd -- use laminar
      //-       (this will only occur for unreasonably small rtheta)
      blData2.diz = di2l;
      blData2.diz_sz = 0.0;
      blData2.diz_uz = di2l_hk2 * blData2.hkz_uz + di2l_rt2 * blData2.rtz_uz;
      blData2.diz_tz = di2l_hk2 * blData2.hkz_tz + di2l_rt2 * blData2.rtz_tz;
      blData2.diz_dz = di2l_hk2 * blData2.hkz_dz;
      blData2.diz_ms = di2l_hk2 * blData2.hkz_ms + di2l_rt2 * blData2.rtz_ms;
      blData2.diz_re = di2l_rt2 * blData2.rtz_re;
    }
  }

  if (ityp == 3) {
    //------ laminar wake cd
    dilw(blData2.hkz, blData2.rtz, di2l, di2l_hk2, di2l_rt2);
    if (di2l > blData2.diz) {
      //------- laminar wake cd is greater than turbulent cd -- use laminar
      //-       (this will only occur for unreasonably small rtheta)
      blData2.diz = di2l;
      blData2.diz_sz = 0.0;
      blData2.diz_uz = di2l_hk2 * blData2.hkz_uz + di2l_rt2 * blData2.rtz_uz;
      blData2.diz_tz = di2l_hk2 * blData2.hkz_tz + di2l_rt2 * blData2.rtz_tz;
      blData2.diz_dz = di2l_hk2 * blData2.hkz_dz;
      blData2.diz_ms = di2l_hk2 * blData2.hkz_ms + di2l_rt2 * blData2.rtz_ms;
      blData2.diz_re = di2l_rt2 * blData2.rtz_re;
    }

    //----- double dissipation for the wake (two wake halves)
    blData2.diz = blData2.diz * 2.0;
    blData2.diz_sz = blData2.diz_sz * 2.0;
    blData2.diz_uz = blData2.diz_uz * 2.0;
    blData2.diz_tz = blData2.diz_tz * 2.0;
    blData2.diz_dz = blData2.diz_dz * 2.0;
    blData2.diz_ms = blData2.diz_ms * 2.0;
    blData2.diz_re = blData2.diz_re * 2.0;
  }

  //---- bl thickness (delta) from simplified green's correlation
  blData2.dez = (3.15 + 1.72 / (blData2.hkz - 1.0)) * blData2.tz + blData2.dz;
  de2_hk2 = (-1.72 / (blData2.hkz - 1.0) / (blData2.hkz - 1.0)) * blData2.tz;

  blData2.dez_uz = de2_hk2 * blData2.hkz_uz;
  blData2.dez_tz = de2_hk2 * blData2.hkz_tz + (3.15 + 1.72 / (blData2.hkz - 1.0));
  blData2.dez_dz = de2_hk2 * blData2.hkz_dz + 1.0;
  blData2.dez_ms = de2_hk2 * blData2.hkz_ms;

  hdmax = 12.0;
  if (blData2.dez > hdmax * blData2.tz) {
    blData2.dez = hdmax * blData2.tz;
    blData2.dez_uz = 0.0;
    blData2.dez_tz = hdmax;
    blData2.dez_dz = 0.0;
    blData2.dez_ms = 0.0;
  }

  return true;
}

/**
 * @brief find index and angle of max foil curv diff
 * 
 * @param x foil x parameters
 * @param y foil y parameters
 * @param n foil plot size
 * @return PairIndex index: index of max angle diff, value: angle of max angle diff[degree]
 */
double XFoil::cang(MatrixX2d points) {
  
  double max_angle = 0;
  //---- go over each point, calculating corner angle
  for (int i = 1; i < points.rows() - 1; i++) {
    Vector2d delta_former = points.row(i) - points.row(i - 1);
    Vector2d delta_later =  points.row(i) - points.row(i + 1);

    double sin = cross2(delta_later, delta_former) / delta_former.norm() / delta_later.norm();
    double delta_angle = asin(sin) * 180.0 / PI;

    max_angle = max(fabs(delta_angle), max_angle);
  }
  return max_angle;
}

bool XFoil::cdcalc() {
  double dx;

  double sa = sin(alfa);
  double ca = cos(alfa);

  if (lvisc && lblini) {
    //---- set variables at the end of the wake
    double thwake = thet[nbl[2]][2];
    double urat = uedg[nbl[2]][2] / qinf;
    double uewake = uedg[nbl[2]][2] * (1.0 - tklam) / (1.0 - tklam * urat * urat);
    double shwake = dstr[nbl[2]][2] / thet[nbl[2]][2];

    //---- extrapolate wake to downstream infinity using squire-young relation
    //      (reduces errors of the wake not being long enough)
    cd = 2.0 * thwake * pow((uewake / qinf), (0.5 * (5.0 + shwake)));
  } else {
    cd = 0.0;
  }

  return true;
}


/**
 * @brief calculate skin friction coefficient(C_f) in lamier
 * 
 * @param hk kinematic shape parameter
 * @param rt momentum-thickness reynolds number
 */
XFoil::C_f XFoil::cfl(double hk, double rt) {
  C_f c_f = C_f();
  if (hk < 5.5) {
    double tmp = pow(5.5 - hk, 3) / (hk + 1.0);
    c_f.cf = (0.0727 * tmp - 0.07) / rt;
    c_f.hk = (-0.0727 * tmp * 3.0 / (5.5 - hk) - 0.0727 * tmp / (hk + 1.0)) / rt;
  } else {
    double tmp = 1.0 - 1.0 / (hk - 4.5);
    c_f.cf = (0.015 * tmp * tmp - 0.07) / rt;
    c_f.hk = 0.015 * tmp * 2.0 / pow(hk - 4.5, 2.0) / rt;
  }
  c_f.rt = -c_f.cf / rt;
  c_f.msq = 0.0;
  return c_f;
}

XFoil::C_f XFoil::cft(double hk, double rt, double msq) {
  C_f c_f = C_f();
  double gam = 1.4;

  //---- turbulent skin friction function  ( cf )    (coles)
  double gm1 = 1.4 - 1.0;
  double fc = sqrt(1.0 + 0.5 * gm1 * msq);
  double grt = std::max(log(rt / fc), 3.0);

  double gex = -1.74 - 0.31 * hk;

  double f_arg = std::max(-1.33 * hk, -20.0);
  
  double tanh_hk = tanh(4.0 - hk / 0.875);

  double cfo = 0.3 * exp(f_arg) * pow((grt / 2.3026), gex);
  c_f.cf = (cfo + 0.00011 * (tanh_hk - 1.0)) / fc;
  c_f.hk = (-1.33 * cfo - 0.31 * log(grt / 2.3026) * cfo -
           0.00011 * (1.0 - pow(tanh_hk, 2)) / 0.875) /
          fc;
  c_f.rt = gex * cfo / (fc * grt) / rt;
  c_f.msq = gex * cfo / (fc * grt) * (-0.25 * gm1 / fc / fc) -
           0.25 * gm1 * (c_f.cf) / pow(fc, 2);

  return c_f;
}

bool XFoil::clcalc(double xref, double yref) {
  // modified techwinder : all other variables are member variables (ex fortran
  // common)
  //-----------------------------------------------------------
  //	   integrates surface pressures to get cl and cm.
  //	   integrates skin friction to get cdf.
  //	   calculates dcl/dalpha for prescribed-cl routines.
  //-----------------------------------------------------------

  // techwinder addition : calculate XCp position

  //---- moment-reference coordinates
  ////ccc	   xref = 0.25
  ////ccc	   yref = 0.

  double beta, beta_msq, bfac, bfac_msq, cginc;
  double cpi_gam, cpc_cpi;
  double cpg1, cpg1_msq, cpg1_alf;
  
  double sa = sin(alfa);
  double ca = cos(alfa);

  xcp = 0.0;

  beta = sqrt(1.0 - minf * minf);
  beta_msq = -0.5 / beta;

  bfac = 0.5 * minf * minf / (1.0 + beta);
  bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;

  cl = 0.0;
  cm = 0.0;

  cl_alf = 0.0;
  cl_msq = 0.0;

  
  cginc = 1.0 - (gam[1] / qinf) * (gam[1] / qinf);
  cpg1 = cginc / (beta + bfac * cginc);
  cpg1_msq = -cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

  cpi_gam = -2.0 * gam[1] / qinf / qinf;
  cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
  cpg1_alf = cpc_cpi * cpi_gam * gam_a[1];

  for (int i = 1; i <= n; i++) {
    int ip = i + 1;
    if (i == n) ip = 1;

    cginc = 1.0 - (gam[ip] / qinf) * (gam[ip] / qinf);
    double cpg2 = cginc / (beta + bfac * cginc);
    double cpg2_msq = -cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

    cpi_gam = -2.0 * gam[ip] / qinf / qinf;
    cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
    double cpg2_alf = cpc_cpi * cpi_gam * gam_a[ip];

    const double dx = (points.row(ip).x() - points.row(i).x()) * ca + (points.row(ip).y() - points.row(i).y()) * sa;
    const double dy = (points.row(ip).y() - points.row(i).y()) * ca - (points.row(ip).x() - points.row(i).x()) * sa;
    const double dg = cpg2 - cpg1;

    const double ax =
        (0.5 * (points.row(ip).x() + points.row(i).x()) - xref) * ca + (0.5 * (points.row(ip).y() + points.row(i).y()) - yref) * sa;
    const double ay =
        (0.5 * (points.row(ip).y() + points.row(i).y()) - yref) * ca - (0.5 * (points.row(ip).x() + points.row(i).x()) - xref) * sa;
    const double ag = 0.5 * (cpg2 + cpg1);

    const double dx_alf = -(points.row(ip).x() - points.row(i).x()) * sa + (points.row(ip).y() - points.row(i).y()) * ca;
    const double ag_alf = 0.5 * (cpg2_alf + cpg1_alf);
    const double ag_msq = 0.5 * (cpg2_msq + cpg1_msq);

    cl = cl + dx * ag;
    cm = cm - dx * (ag * ax + dg * dx / 12.0) - dy * (ag * ay + dg * dy / 12.0);

    xcp += dx * ag * (points.row(ip).x() + points.row(i).x()) / 2.0;

    cl_alf = cl_alf + dx * ag_alf + ag * dx_alf;
    cl_msq = cl_msq + dx * ag_msq;

    cpg1 = cpg2;
    cpg1_alf = cpg2_alf;
    cpg1_msq = cpg2_msq;
  }

  if (fabs(cl) > 0.0)
    xcp /= cl;
  else
    xcp = 0.0;

  return true;
}

bool XFoil::comset() {
  //---- set karman-tsien parameter tklam
  double beta, beta_msq;
  beta = sqrt(1.0 - minf * minf);
  beta_msq = -0.5 / beta;

  tklam = minf * minf / (1.0 + beta) / (1.0 + beta);
  tkl_msq =
      1.0 / (1.0 + beta) / (1.0 + beta) - 2.0 * tklam / (1.0 + beta) * beta_msq;

  //---- set sonic pressure coefficient and speed
  if (minf == 0.0) {
    cpstar = -999.0;
    qstar = 999.0;
  } else {
    cpstar = 2.0 / (gamma * minf * minf) *
             (pow(((1.0 + 0.5 * gamm1 * minf * minf) / (1.0 + 0.5 * gamm1)),
                  (gamma / gamm1)) -
              1.0);
    qstar = qinf / minf *
            sqrt((1.0 + 0.5 * gamm1 * minf * minf) / (1.0 + 0.5 * gamm1));
  }

  return true;
}

/** ---------------------------------------------
 *      sets compressible cp from speed.
 * ---------------------------------------------- */
bool XFoil::cpcalc(int n, const double q[], double qinf, double minf, double cp[]) {

  bool denneg;
  double beta, bfac;

  beta = sqrt(1.0 - minf * minf);
  bfac = 0.5 * minf * minf / (1.0 + beta);

  denneg = false;

  for (int i = 1; i <= n; i++) {
    const double cpinc = 1.0 - (q[i] / qinf) * (q[i] / qinf);
    const double den = beta + bfac * cpinc;
    cp[i] = cpinc / den;
    if (den <= 0.0) denneg = true;
  }

  if (denneg) {
    writeString(
        "CpCalc: local speed too larger\n Compressibility corrections "
        "invalid\n",
        true);
    return false;
  }

  return true;
}

void XFoil::writeString(std::string str, bool bFullReport) {
  *m_pOutStream << str;
}

/** ==============================================================
 *      amplification rate routine for envelope e^n method.
 *      reference:
 *                 drela, m., giles, m.,
 *                "viscous/inviscid analysis of transonic and
 *                 low reynolds number airfoils",
 *                 aiaa journal, oct. 1987.
 *
 *      new version.   march 1991       (latest bug fix  july 93)
 *           - m(h) correlation made more accurate up to h=20
 *           - for h > 5, non-similar profiles are used
 *             instead of falkner-skan profiles.  these
 *             non-similar profiles have smaller reverse
 *             velocities, are more representative of typical
 *             separation bubble profiles.
 * --------------------------------------------------------------
 *
 *      input :   hk     kinematic shape parameter
 *                th     momentum thickness
 *                rt     momentum-thickness reynolds number
 *
 *      output:   ax     envelope spatial amplification rate
 *                ax_(.) sensitivity of ax to parameter (.)
 *
 *
 *      usage: the log of the envelope amplitude n(x) is
 *             calculated by integrating ax (= dn/dx) with
 *             respect to the streamwise distance x.
 *                       x
 *                      /
 *               n(x) = | ax(h(x),th(x),rth(x)) dx
 *                      /
 *                       0
 *             the integration can be started from the leading
 *             edge since ax will be returned as zero when rt
 *             is below the critical rtheta.  transition occurs
 *             when n(x) reaches ncrit (ncrit= 9 is "standard").
 * ============================================================== */
bool XFoil::dampl(double hk, double th, double rt, double &ax, double &ax_hk,
                  double &ax_th, double &ax_rt) {
  double dgr = 0.08;

  double hmi = 0.0, hmi_hk = 0.0, aa = 0.0, aa_hk = 0.0, bb = 0.0, bb_hk = 0.0,
         grcrit = 0.0, grc_hk = 0.0, gr = 0.0, gr_rt = 0.0;

  hmi = 1.0 / (hk - 1.0);
  hmi_hk = -hmi * hmi;

  //---- log10(critical rth) - h   correlation for falkner-skan profiles
  aa = 2.492 * pow(hmi, 0.43);
  aa_hk = (aa / hmi) * 0.43 * hmi_hk;
  bb = tanh(14.0 * hmi - 9.24);
  bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk;
  grcrit = aa + 0.7 * (bb + 1.0);
  grc_hk = aa_hk + 0.7 * bb_hk;
  gr = log10(rt);
  gr_rt = 1.0 / (2.3025851 * rt);
  if (gr < grcrit - dgr) {
    //----- no amplification for rtheta < rcrit
    ax = 0.0;
    ax_hk = 0.0;
    ax_th = 0.0;
    ax_rt = 0.0;
  } else {
    //----- set steep cubic ramp used to turn on ax smoothly as rtheta
    //-     exceeds rcrit (previously, this was done discontinuously).
    //-     the ramp goes between  -dgr < log10(rtheta/rcrit) < dgr

    const double rnorm = (gr - (grcrit - dgr)) / (2.0 * dgr);
    double rn_hk = -grc_hk / (2.0 * dgr);
    double rn_rt = gr_rt / (2.0 * dgr);

    double rfac, rfac_hk, rfac_rt;
    if (rnorm >= 1.0) {
      rfac = 1.0;
      rfac_hk = 0.0;
      rfac_rt = 0.0;
    } else {
      rfac = 3.0 * rnorm * rnorm - 2.0 * rnorm * rnorm * rnorm;
      const double rfac_rn = 6.0 * rnorm - 6.0 * rnorm * rnorm;

      rfac_hk = rfac_rn * rn_hk;
      rfac_rt = rfac_rn * rn_rt;
    }

    //----- amplification envelope slope correlation for falkner-skan
    const double f_arg = 3.87 * hmi - 2.52;
    const double arg_hk = 3.87 * hmi_hk;

    const double ex = exp(-f_arg * f_arg);
    const double ex_hk = ex * (-2.0 * f_arg * arg_hk);

    const double dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
    const double dadr_hk = 0.028 - 0.0345 * ex_hk;

    //----- new m(h) correlation    1 march 91
    const double af = -0.05 + 2.7 * hmi - 5.5 * hmi * hmi + 3.0 * hmi * hmi * hmi;
    const double af_hmi = 2.7 - 11.0 * hmi + 9.0 * hmi * hmi;
    const double af_hk = af_hmi * hmi_hk;

    ax = (af * dadr / th) * rfac;
    ax_hk = (af_hk * dadr / th + af * dadr_hk / th) * rfac +
            (af * dadr / th) * rfac_hk;
    ax_th = -(ax) / th;
    ax_rt = (af * dadr / th) * rfac_rt;
  }

  return true;
}

/** Laminar dissipation function  ( 2 cd/h* )     (from Falkner-Skan)*/
bool XFoil::dil(double hk, double rt, double &di, double &di_hk,
                double &di_rt) {
  if (hk < 4.0) {
    di = (0.00205 * pow((4.0 - hk), 5.5) + 0.207) / rt;
    di_hk = (-.00205 * 5.5 * pow((4.0 - hk), 4.5)) / rt;
  } else {
    double hkb = hk - 4.0;
    double den = 1.0 + 0.02 * hkb * hkb;
    di = (-.0016 * hkb * hkb / den + 0.207) / rt;
    di_hk =
        (-.0016 * 2.0 * hkb * (1.0 / den - 0.02 * hkb * hkb / den / den)) / rt;
  }
  di_rt = -(di) / rt;

  return true;
}

bool XFoil::dilw(double hk, double rt, double &di, double &di_hk,
                 double &di_rt) {
  //	double msq = 0.0;
  double hs, hs_hk, hs_rt, hs_msq;

  hsl(hk, hs, hs_hk, hs_rt, hs_msq);
  //---- laminar wake dissipation function  ( 2 cd/h* )
  double rcd = 1.10 * (1.0 - 1.0 / hk) * (1.0 - 1.0 / hk) / hk;
  double rcd_hk = -1.10 * (1.0 - 1.0 / hk) * 2.0 / hk / hk / hk - rcd / hk;

  di = 2.0 * rcd / (hs * rt);
  di_hk = 2.0 * rcd_hk / (hs * rt) - ((di) / hs) * hs_hk;
  di_rt = -(di) / rt - ((di) / hs) * hs_rt;

  return true;
}

bool XFoil::dslim(double &dstr, double thet, double msq, double hklim) {
  double h, hk, hk_h, hk_m, dh;
  h = (dstr) / thet;

  hkin(h, msq, hk, hk_h, hk_m);

  dh = std::max(0.0, hklim - hk) / hk_h;
  dstr = (dstr) + dh * thet;

  return true;
}

/** ------------------------------------------------
 *     finds minimum cp on dist for cavitation work
 * ------------------------------------------------ */
bool XFoil::fcpmin() {
  xcpmni = points.row(1).x();
  xcpmnv = points.row(1).x();
  cpmni = cpi[1];
  cpmnv = cpv[1];

  for (int i = 2; i <= n + nw; i++) {
    if (cpi[i] < cpmni) {
      xcpmni = points.row(i).x();
      cpmni = cpi[i];
    }
    if (cpv[i] < cpmnv) {
      xcpmnv = points.row(i).x();
      cpmnv = cpv[i];
    }
  }

  if (lvisc)
    cpmn = cpmnv;
  else {
    cpmn = cpmni;
    cpmnv = cpmni;
    xcpmnv = xcpmni;
  }

  return true;
}

bool XFoil::gamqv() {
  for (int i = 1; i <= n; i++) {
    gam[i] = qvis[i];
    gam_a[i] = qinv_a[i];
  }

  return true;
}

/**
 *   Solves general nxn system in nn unknowns
 *   with arbitrary number (nrhs) of righthand sides.
 *   assumes system is invertible...
 *    ...if it isn't, a divide by zero will result.
 *
 *   z is the coefficient matrix...
 *     ...destroyed during solution process.
 *   r is the righthand side(s)...
 *     ...replaced by the solution vector(s).
 *
 *                              mark drela  1984
 */
bool XFoil::Gauss(int nn, double z[][6], double r[5]) {
  // dimension z(nsiz,nsiz), r(nsiz,nrhs)

  int loc;
  int np, nnpp, nt, k;

  double temp, ztmp;

  for (np = 1; np <= nn - 1; np++) {
    nnpp = np + 1;
    //------ find max pivot index nx
    int nx = np;
    for (nt = nnpp; nt <= nn; nt++) {
      if (fabs(z[nt][np]) > fabs(z[nx][np])) nx = nt;
    }

    double pivot = 1.0 / z[nx][np];

    //------ switch pivots
    z[nx][np] = z[np][np];

    //------ switch rows & normalize pivot row
    for (loc = nnpp; loc <= nn; loc++) {
      temp = z[nx][loc] * pivot;
      z[nx][loc] = z[np][loc];
      z[np][loc] = temp;
    }

    temp = r[nx] * pivot;
    r[nx] = r[np];
    r[np] = temp;

    //------ forward eliminate everything
    for (k = nnpp; k <= nn; k++) {
      ztmp = z[k][np];
      for (loc = nnpp; loc <= nn; loc++)
        z[k][loc] = z[k][loc] - ztmp * z[np][loc];
      r[k] = r[k] - ztmp * r[np];
    }
  }

  //---- solve for last row
  r[nn] = r[nn] / z[nn][nn];

  //---- back substitute everything
  for (np = nn - 1; np >= 1; np--) {
    nnpp = np + 1;
    for (k = nnpp; k <= nn; k++) r[np] = r[np] - z[np][k] * r[k];
  }

  return true;
}

bool XFoil::getxyf(MatrixX2d points, MatrixX2d dpoints_ds, VectorXd s,
                   int n, double &tops, double &bots, double xf, double &yf) {
  double topy, boty, yrel;

  tops = s[1] + (points.row(1).x() - xf);
  bots = s[n] - (points.row(n).x() - xf);
  sinvrt(tops, xf, points.col(0), dpoints_ds.col(0), s, n);
  sinvrt(bots, xf, points.col(0), dpoints_ds.col(0), s, n);
  topy = spline::seval(tops, points.col(1), dpoints_ds.col(1), s, n);
  boty = spline::seval(bots, points.col(1), dpoints_ds.col(1), s, n);

  yrel = yf;

  yf = topy * yrel + boty * (1.0 - yrel);
  return true;
}

/** --------------------------------------------------------------
 *     Calculates two surface vorticity (gamma) distributions
 *     for alpha = 0, 90  degrees.  These are superimposed
 *     in specal or speccl for specified alpha or cl.
 *-------------------------------------------------------------- */
bool XFoil::ggcalc() {

  double psi, psi_n, res;
  double bbb[IQX];
  //	double psiinf;

  cosa = cos(alfa);
  sina = sin(alfa);

  writeString("   Calculating unit vorticity distributions ...\n");

  for (int i = 1; i <= n; i++) {
    gam[i] = 0.0;
    gamu[i][1] = 0.0;
    gamu[i][2] = 0.0;
  }
  psio = 0.0;

  //---- set up matrix system for  psi = psio  on airfoil surface.
  //-    the unknowns are (dgamma)i and dpsio.
  for (int i = 1; i <= n; i++) {
    //------ calculate psi and dpsi/dgamma array for current node
    psilin(i, points.row(i).x(), points.row(i).y(), nx[i], ny[i], psi, psi_n, false, true);

    const double res1 = qinf * points.row(i).y();
    const double res2 = -qinf * points.row(i).x();

    //------ dres/dgamma
    for (int j = 1; j <= n; j++) {
      aij[i][j] = dzdg[j];
    }

    for (int j = 1; j <= n; j++) {
      bij[i][j] = -dzdm[j];
    }

    //------ dres/dpsio
    aij[i][n + 1] = -1.0;

    gamu[i][1] = -res1;
    gamu[i][2] = -res2;
  }

  //---- set Kutta condition
  //-    res = gam(1) + gam[n]
  res = 0.0;

  for (int j = 1; j <= n + 1; j++) aij[n + 1][j] = 0.0;

  aij[n + 1][1] = 1.0;
  aij[n + 1][n] = 1.0;

  gamu[n + 1][1] = -res;
  gamu[n + 1][2] = -res;

  //---- set up Kutta condition (no direct source influence)
  for (int j = 1; j <= n; j++) bij[n + 1][j] = 0.0;

  if (sharp) {
    //----- set zero internal velocity in TE corner

    //----- set TE bisector angle
    const double ag1 = atan2(-dpoints_ds.row(1).y(), -dpoints_ds.row(1).x());
    const double ag2 = atanc(dpoints_ds.row(n).y(), dpoints_ds.row(n).x(), ag1);
    const double abis = 0.5 * (ag1 + ag2);

    const double cbis = cos(abis);
    const double sbis = sin(abis);

    //----- minimum panel length adjacent to TE
    const double ds1 = sqrt((points.row(1).x() - points.row(2).x()) * (points.row(1).x() - points.row(2).x()) + (points.row(1).y() - points.row(2).y()) * (points.row(1).y() - points.row(2).y()));
    const double ds2 = sqrt((points.row(n).x() - points.row(n - 1).x()) * (points.row(n).x() - points.row(n - 1).x()) +
               (points.row(n).y() - points.row(n - 1).y()) * (points.row(n).y() - points.row(n - 1).y()));
    const double dsmin = std::min(ds1, ds2);

    //---- distance of internal control point ahead of sharp TE
    //-    (fraction of smaller panel length adjacent to TE)
    const double bwt = 0.1;

    //----- control point on bisector just ahead of TE point
    const double xbis = point_te.x() - bwt * dsmin * cbis;
    const double ybis = point_te.y() - bwt * dsmin * sbis;

    //----- set velocity component along bisector line
    double qbis;
    psilin(0, xbis, ybis, -sbis, cbis, psi, qbis, false, true);

    //----- dres/dgamma
    for (int j = 1; j <= n; j++) aij[n][j] = dqdg[j];

    //----- -dres/dmass
    for (int j = 1; j <= n; j++) bij[n][j] = -dqdm[j];

    //----- dres/dpsio
    aij[n][n + 1] = 0.0;

    //----- -dres/duinf
    gamu[n][1] = -cbis;

    //----- -dres/dvinf
    gamu[n][2] = -sbis;
  }

  //---- lu-factor coefficient matrix aij
  ludcmp(n + 1, aij, aijpiv);
  lqaij = true;

  //---- solve system for the two vorticity distributions
  for (int iu = 0; iu < IQX; iu++)
    bbb[iu] = gamu[iu][1];  // techwinder : create a dummy array
  baksub(n + 1, aij, aijpiv, bbb);
  for (int iu = 0; iu < IQX; iu++) gamu[iu][1] = bbb[iu];

  for (int iu = 0; iu < IQX; iu++)
    bbb[iu] = gamu[iu][2];  // techwinder : create a dummy array
  baksub(n + 1, aij, aijpiv, bbb);
  for (int iu = 0; iu < IQX; iu++) gamu[iu][2] = bbb[iu];

  //---- set inviscid alpha=0,90 surface speeds for this geometry
  for (int i = 1; i <= n + 1; i++) {
    qinvu[i][1] = gamu[i][1];
    qinvu[i][2] = gamu[i][2];
  }

  lgamu = true;

  return true;
}

bool XFoil::baksub(int n, double a[IQX][IQX], const int indx[], double b[]) {
  double sum;
  int ii;
  ii = 0;
  for (int i = 1; i <= n; i++) {
    int ll = indx[i];
    sum = b[ll];
    b[ll] = b[i];
    if (ii != 0)
      for (int j = ii; j <= i - 1; j++) sum = sum - a[i][j] * b[j];
    else if (sum != 0.0)
      ii = i;

    b[i] = sum;
  }
  //
  for (int i = n; i >= 1; i--) {
    sum = b[i];
    if (i < n)
      for (int j = i + 1; j <= n; j++) sum = sum - a[i][j] * b[j];

    b[i] = sum / a[i][i];
  }
  //
  return true;
}

bool XFoil::hct(double hk, double msq, double &hc, double &hc_hk,
                double &hc_msq) {
  //---- density shape parameter    (from whitfield)
  hc = msq * (0.064 / (hk - 0.8) + 0.251);
  hc_hk = msq * (-.064 / (hk - 0.8) / (hk - 0.8));
  hc_msq = 0.064 / (hk - 0.8) + 0.251;

  return true;
}

bool XFoil::hkin(double h, double msq, double &hk, double &hk_h,
                 double &hk_msq) {
  //---- calculate kinematic shape parameter (assuming air)
  //     (from Whitfield )
  hk = (h - 0.29 * msq) / (1.0 + 0.113 * msq);
  hk_h = 1.0 / (1.0 + 0.113 * msq);
  hk_msq = (-.29 - 0.113 * (hk)) / (1.0 + 0.113 * msq);

  return true;
}

bool XFoil::hsl(double hk, double &hs, double &hs_hk, double &hs_rt,
                double &hs_msq) {
  //---- laminar hs correlation
  if (hk < 4.35) {
    double tmp = hk - 4.35;
    hs = 0.0111 * tmp * tmp / (hk + 1.0) -
         0.0278 * tmp * tmp * tmp / (hk + 1.0) + 1.528 -
         0.0002 * (tmp * hk) * (tmp * hk);
    hs_hk =
        0.0111 * (2.0 * tmp - tmp * tmp / (hk + 1.0)) / (hk + 1.0) -
        0.0278 * (3.0 * tmp * tmp - tmp * tmp * tmp / (hk + 1.0)) / (hk + 1.0) -
        0.0002 * 2.0 * tmp * hk * (tmp + hk);
  } else {
    hs = 0.015 * (hk - 4.35) * (hk - 4.35) / hk + 1.528;
    hs_hk = 0.015 * 2.0 * (hk - 4.35) / hk -
            0.015 * (hk - 4.35) * (hk - 4.35) / hk / hk;
  }

  hs_rt = 0.0;
  hs_msq = 0.0;

  return true;
}

bool XFoil::hst(double hk, double rt, double msq, double &hs, double &hs_hk,
                double &hs_rt, double &hs_msq) {
  double hsmin, dhsinf, rtz, rtz_rt, ho, ho_rt, fm;

  //---- turbulent hs correlation

  hsmin = 1.5;
  dhsinf = 0.015;

  //---- ###  12/4/94
  //---- limited rtheta dependence for rtheta < 200

  if (rt > 400.0) {
    ho = 3.0 + 400.0 / rt;
    ho_rt = -400.0 / rt / rt;
  } else {
    ho = 4.0;
    ho_rt = 0.0;
  }

  if (rt > 200.0) {
    rtz = rt;
    rtz_rt = 1.0;
  } else {
    rtz = 200.0;
    rtz_rt = 0.0;
  }

  if (hk < ho) {
    //----- attached branch
    //----- new correlation  29 nov 91
    //-     (from  arctan(y+) + schlichting  profiles)
    const double hr = (ho - hk) / (ho - 1.0);
    const double hr_hk = -1.0 / (ho - 1.0);
    const double hr_rt = (1.0 - hr) / (ho - 1.0) * ho_rt;
    hs = (2.0 - hsmin - 4.0 / rtz) * hr * hr * 1.5 / (hk + 0.5) + hsmin +
         4.0 / rtz;
    hs_hk =
        -(2.0 - hsmin - 4.0 / rtz) * hr * hr * 1.5 / (hk + 0.5) / (hk + 0.5) +
        (2.0 - hsmin - 4.0 / rtz) * hr * 2.0 * 1.5 / (hk + 0.5) * hr_hk;
    hs_rt = (2.0 - hsmin - 4.0 / rtz) * hr * 2.0 * 1.5 / (hk + 0.5) * hr_rt +
            (hr * hr * 1.5 / (hk + 0.5) - 1.0) * 4.0 / rtz / rtz * rtz_rt;
  } else {
    //----- separated branch
    const double grt = log(rtz);
    const double hdif = hk - ho;
    const double rtmp = hk - ho + 4.0 / grt;
    const double htmp = 0.007 * grt / rtmp / rtmp + dhsinf / hk;
    const double htmp_hk = -.014 * grt / rtmp / rtmp / rtmp - dhsinf / hk / hk;
    const double htmp_rt = -.014 * grt / rtmp / rtmp / rtmp *
                  (-ho_rt - 4.0 / grt / grt / rtz * rtz_rt) +
              0.007 / rtmp / rtmp / rtz * rtz_rt;
    hs = hdif * hdif * htmp + hsmin + 4.0 / rtz;
    hs_hk = hdif * 2.0 * htmp + hdif * hdif * htmp_hk;
    hs_rt = hdif * hdif * htmp_rt - 4.0 / rtz / rtz * rtz_rt +
            hdif * 2.0 * htmp * (-ho_rt);
  }

  //---- whitfield's minor additional compressibility correction
  fm = 1.0 + 0.014 * msq;
  hs = (hs + 0.028 * msq) / fm;
  hs_hk = (hs_hk) / fm;
  hs_rt = (hs_rt) / fm;
  hs_msq = 0.028 / fm - 0.014 * (hs) / fm;

  return true;
}

/** -----------------------------------------------------------
 *     sets  bl location -> panel location  pointer array ipan
 * -----------------------------------------------------------*/
bool XFoil::iblpan() {
  int iblmax, is, ibl;
  std::stringstream ss;

  //-- top surface first
  is = 1;

  ibl = 1;
  for (int i = ist; i >= 1; i--) {
    ibl = ibl + 1;
    ipan[ibl][is] = i;
    vti[ibl][is] = 1.0;
  }

  iblte[is] = ibl;
  nbl[is] = ibl;

  //-- bottom surface next
  is = 2;
  ibl = 1;
  for (int i = ist + 1; i <= n; i++) {
    ibl = ibl + 1;
    ipan[ibl][is] = i;
    vti[ibl][is] = -1.0;
  }

  //-- wake
  iblte[is] = ibl;

  for (int iw = 1; iw <= nw; iw++) {
    int i = n + iw;
    ibl = iblte[is] + iw;
    ipan[ibl][is] = i;
    vti[ibl][is] = -1.0;
  }

  nbl[is] = iblte[is] + nw;

  //-- upper wake pointers (for plotting only)
  for (int iw = 1; iw <= nw; iw++) {
    ipan[iblte[1] + iw][1] = ipan[iblte[2] + iw][2];
    vti[iblte[1] + iw][1] = 1.0;
  }
  iblmax = std::max(iblte[1], iblte[2]) + nw;
  if (iblmax > IVX) {
    ss << "iblpan :  ***  bl array overflow\n";
    ss << "Increase IVX to at least " << iblmax << "\n";
    writeString(ss.str(), true);
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
    for (int ibl = 2; ibl <= nbl[is]; ibl++) {
      iv++;
      isys[ibl][is] = iv;
    }
  }

  nsys = iv;
  if (nsys > 2 * IVX) {
    writeString("*** iblsys: bl system array overflow. ***", true);
    return false;
  }

  return true;
}

/** Loads the Foil's geometry in XFoil,
 *  calculates the normal vectors,
 *  and sets the results in current foil */
bool XFoil::initXFoilGeometry(int fn, const double *fx, const double *fy, double *fnx, double *fny) {

  MatrixX2d buffer_points = MatrixX2d::Zero(fn + 1, 2);
  for (int i = 0; i < fn; i++) {
    buffer_points.row(i + 1).x() = fx[i];
    buffer_points.row(i + 1).y() = fy[i];
  }

  if (!isValidFoilPointSize(buffer_points) || !isValidFoilAngles(buffer_points)) {
    writeString("Unrecognized foil format");
    return false;
  }

  abcopy(buffer_points);
  return true;
}

bool XFoil::initXFoilAnalysis(double Re, double alpha, double Mach,
                              double NCrit, double XtrTop, double XtrBot,
                              ReynoldsType reType, MachType maType, bool bViscous,
                              std::stringstream &outStream) {
  // Sets Analysis parameters in XFoil
  m_pOutStream = &outStream;

  lblini = false;
  lipan = false;

  reinf1 = Re;
  alfa = alpha * PI / 180.0;

  minf1 = Mach;
  reynolds_type = reType;
  mach_type = maType;
  lalfa = true;
  qinf = 1.0;

  lvisc = bViscous;

  acrit = NCrit;
  xstrip[1] = XtrTop;
  xstrip[2] = XtrBot;

  if (Mach > 0.000001) {
    if (!setMach()) {
      writeString(
          "... Invalid Analysis Settings\nCpCalc: local speed too large\n "
          "Compressibility corrections invalid ");
      return false;
    }
  }

  return true;
}

/** ------------------------------------------------------
 *     locates leading edge spline-parameter value sle
 *
 *     the defining condition is
 *
 *      (x-xte,y-yte) . (x',y') = 0     at  s = sle
 *
 *     i.e. the surface tangent is normal to the chord
 *     line connecting x(sle),y(sle) and the te point.
//------------------------------------------------------ */
bool XFoil::lefind(double &sle, MatrixX2d points, MatrixX2d dpoints_ds, VectorXd s, int n) {
  int i;
  double dseps;
  //---- convergence tolerance
  dseps = (s[n] - s[1]) * 0.00001;

  //---- set trailing edge point coordinates
  point_te.x() = 0.5 * (points.row(1).x() + points.row(n).x());
  point_te.y() = 0.5 * (points.row(1).y() + points.row(n).y());

  //---- get first guess for sle
  for (i = 3; i <= n - 2; i++) {
    const double dxte = points.row(i).x() - point_te.x();
    const double dyte = points.row(i).y() - point_te.y();
    const double dx = points.row(i + 1).x() - points.row(i).x();
    const double dy = points.row(i + 1).y() - points.row(i).y();
    const double dotp = dxte * dx + dyte * dy;
    if (dotp < 0.0) break;
  }

  sle = s[i];

  //---- check for sharp le case
  if (s[i] == s[i - 1]) return false;

  //---- newton iteration to get exact sle value
  for (int iter = 1; iter <= 50; iter++) {
    point_le.x() = spline::seval(sle, points.col(0), dpoints_ds.col(0), s, n);
    point_le.y() = spline::seval(sle, points.col(1), dpoints_ds.col(1), s, n);
    const double dxds = spline::deval(sle, points.col(0), dpoints_ds.col(0), s, n);
    const double dyds = spline::deval(sle, points.col(1), dpoints_ds.col(1), s, n);
    const double dxdd = spline::d2val(sle, points.col(0), dpoints_ds.col(0), s, n);
    const double dydd = spline::d2val(sle, points.col(1), dpoints_ds.col(1), s, n);

    Vector2d point_chord = point_le - point_te;

    const double xchord = point_le.x() - point_te.x();
    const double ychord = point_le.y() - point_te.y();

    //------ drive dot product between chord line and le tangent to zero
    const double res = xchord * dxds + ychord * dyds;
    const double ress = dxds * dxds + dyds * dyds + xchord * dxdd + ychord * dydd;

    //------ newton delta for sle
    double dsle = -res / ress;

    dsle = std::max(dsle, -0.02 * fabs(xchord + ychord));
    dsle = std::min(dsle, 0.02 * fabs(xchord + ychord));
    sle = sle + dsle;
    if (fabs(dsle) < dseps) return true;
  }

  sle = s[i];
  return true;
}

/**    *******************************************************
 *    *                                                     *
 *    *   factors a full nxn matrix into an lu form.        *
 *    *   subr. baksub can back-substitute it with some rhs.*
 *    *   assumes matrix is non-singular...                 *
 *    *    ...if it isn"t, a divide by zero will result.    *
 *    *                                                     *
 *    *   a is the matrix...                                *
 *    *     ...replaced with its lu factors.                *
 *    *                                                     *
 *    *                              mark drela  1988       *
 *    *******************************************************
 */

bool XFoil::ludcmp(int n, double a[IQX][IQX], int indx[IQX]) {
  int imax = 0;  // added techwinder
  int nvx = IQX;

  double vv[IQX];
  double dum, sum, aamax;
  if (n > nvx) {
    writeString("Stop ludcmp: array overflow. Increase nvx", true);
    return false;
  }

  for (int i = 1; i <= n; i++) {
    aamax = 0.0;
    for (int j = 1; j <= n; j++) aamax = std::max(fabs(a[i][j]), aamax);
    vv[i] = 1.0 / aamax;
  }

  for (int j = 1; j <= n; j++) {
    for (int i = 1; i <= j - 1; i++) {
      sum = a[i][j];
      for (int k = 1; k <= i - 1; k++) sum = sum - a[i][k] * a[k][j];
      a[i][j] = sum;
    }

    aamax = 0.0;

    for (int i = j; i <= n; i++) {
      sum = a[i][j];
      for (int k = 1; k <= j - 1; k++) sum = sum - a[i][k] * a[k][j];
      a[i][j] = sum;
      dum = (vv[i] * fabs(sum));
      if (dum >= aamax) {
        imax = i;
        aamax = dum;
      }
    }
    //		ASSERT(bimaxok);// to check imax has been initialized
    if (j != imax) {
      for (int k = 1; k <= n; k++) {
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      vv[imax] = vv[j];
    }

    indx[j] = imax;
    if (j != n) {
      dum = 1.0 / a[j][j];
      for (int i = j + 1; i <= n; i++) a[i][j] = a[i][j] * dum;
    }
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
bool XFoil::mrchdu() {
  std::stringstream ss;
  double vtmp[5][6], vztmp[5];
  memset(vtmp, 0, 30 * sizeof(double));
  memset(vztmp, 0, 5 * sizeof(double));
  double deps = 0.000005;
  int is = 0, ibl, ibm, itrold, iw = 0, itbl;  // icom

  double senswt = 0.0, thi, uei, dsi, cti, dswaki, ratlen = 0.0;
  double sens = 0.0, sennew = 0.0, dummy = 0.0, msq = 0.0, thm = 0.0, dsm = 0.0,
         uem = 0.0;
  double xsi, hklim = 0.0, dsw = 0.0;
  double ami = 0.0, dte = 0.0, cte = 0.0, tte = 0.0, ueref = 0.0, hkref = 0.0,
         dmax = 0.0;

  //---- constant controlling how far hk is allowed to deviate
  //    from the specified value.
  senswt = 1000.0;
  sens = 0.0;
  sennew = 0.0;

  for (is = 1; is <= 2; is++) {  // 2000

    //---- set forced transition arc length position
    xifset(is);

    //---- set leading edge pressure gradient parameter  x/u du/dx
    ibl = 2;
    bule = 1.0;

    //---- old transition station
    itrold = itran[is];

    tran = false;
    turb = false;
    itran[is] = iblte[is];
    //---- march downstream
    for (ibl = 2; ibl <= nbl[is]; ibl++) {
      ibm = ibl - 1;

      simi = ibl == 2;
      wake = ibl > iblte[is];

      //------ initialize current station to existing variables
      xsi = xssi[ibl][is];
      uei = uedg[ibl][is];
      thi = thet[ibl][is];
      dsi = dstr[ibl][is];

      //------ fixed bug   md 7 june 99
      if (ibl < itrold) {
        ami = ctau[ibl][is];  // ami must be initialized
        cti = 0.03;
      } else {
        cti = ctau[ibl][is];
        if (cti <= 0.0) cti = 0.03;
      }

      if (wake) {
        iw = ibl - iblte[is];
        dswaki = wgap[iw];
      } else
        dswaki = 0.0;

      if (ibl <= iblte[is])
        dsi = std::max(dsi - dswaki, 1.02000 * thi) + dswaki;
      if (ibl > iblte[is]) dsi = std::max(dsi - dswaki, 1.00005 * thi) + dswaki;

      //------ newton iteration loop for current station

      for (itbl = 1; itbl <= 25; itbl++) {  // 100

        //-------- assemble 10x3 linearized system for dctau, dth, dds, due, dxi
        //         at the previous "1" station and the current "2" station
        //         (the "1" station coefficients will be ignored)

        blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
        blkin();

        //-------- check for transition and set appropriate flags and things
        if ((!simi) && (!turb)) {
          trchek();
          ami = blData2.amplz;
          if (tran) itran[is] = ibl;
          if (!tran) itran[is] = ibl + 2;
        }
        if (ibl == iblte[is] + 1) {
          tte = thet[iblte[1]][1] + thet[iblte[2]][2];
          dte = dstr[iblte[1]][1] + dstr[iblte[2]][2] + ante;
          cte = (ctau[iblte[1]][1] * thet[iblte[1]][1] +
                 ctau[iblte[2]][2] * thet[iblte[2]][2]) /
                tte;
          tesys(cte, tte, dte);
        } else {
          blsys();
        }

        //-------- set stuff at first iteration...
        if (itbl == 1) {
          //--------- set "baseline" ue and hk for forming  ue(hk)  relation
          ueref = blData2.uz;
          hkref = blData2.hkz;

          //--------- if current point ibl was turbulent and is now laminar,
          // then...
          if (ibl < itran[is] && ibl >= itrold) {
            //---------- extrapolate baseline hk
            uem = uedg[ibl - 1][is];
            dsm = dstr[ibl - 1][is];
            thm = thet[ibl - 1][is];
            msq =
                uem * uem * hstinv / (gm1bl * (1.0 - 0.5 * uem * uem * hstinv));
            hkin(dsm / thm, msq, hkref, dummy, dummy);
          }

          //--------- if current point ibl was laminar, then...
          if (ibl < itrold) {
            //---------- reinitialize or extrapolate ctau if it's now turbulent
            if (tran) ctau[ibl][is] = 0.03;
            if (turb) ctau[ibl][is] = ctau[ibl - 1][is];
            if (tran || turb) {
              cti = ctau[ibl][is];
              blData2.sz = cti;
            }
          }
        }

        if (simi || ibl == iblte[is] + 1) {
          //--------- for similarity station or first wake point, prescribe ue
          vs2[4][1] = 0.0;
          vs2[4][2] = 0.0;
          vs2[4][3] = 0.0;
          vs2[4][4] = blData2.uz_uei;
          vsrez[4] = ueref - blData2.uz;
        } else {
          //******** calculate ue-hk characteristic slope
          for (int k = 1; k <= 4; k++) {
            vztmp[k] = vsrez[k];
            for (int l = 1; l <= 5; l++) vtmp[k][l] = vs2[k][l];
          }

          //--------- set unit dhk
          vtmp[4][1] = 0.0;
          vtmp[4][2] = blData2.hkz_tz;
          vtmp[4][3] = blData2.hkz_dz;
          vtmp[4][4] = blData2.hkz_uz * blData2.uz_uei;
          vztmp[4] = 1.0;

          //--------- calculate due response
          Gauss(4, vtmp, vztmp);

          //--------- set  senswt * (normalized due/dhk)
          sennew = senswt * vztmp[4] * hkref / ueref;
          if (itbl <= 5)
            sens = sennew;
          else if (itbl <= 15)
            sens = 0.5 * (sens + sennew);

          //--------- set prescribed ue-hk combination
          vs2[4][1] = 0.0;
          vs2[4][2] = blData2.hkz_tz * hkref;
          vs2[4][3] = blData2.hkz_dz * hkref;
          vs2[4][4] = (blData2.hkz_uz * hkref + sens / ueref) * blData2.uz_uei;
          vsrez[4] = -(hkref * hkref) * (blData2.hkz / hkref - 1.0) -
                     sens * (blData2.uz / ueref - 1.0);
        }

        //-------- solve newton system for current "2" station
        Gauss(4, vs2, vsrez);

        //-------- determine max changes and underrelax if necessary
        dmax = std::max(fabs(vsrez[2] / thi), fabs(vsrez[3] / dsi));
        if (ibl >= itran[is])
          dmax = std::max(dmax, fabs(vsrez[1] / (10.0 * cti)));

        rlx = 1.0;
        if (dmax > 0.3) rlx = 0.3 / dmax;

        //-------- update as usual
        if (ibl < itran[is]) ami = ami + rlx * vsrez[1];
        if (ibl >= itran[is]) cti = cti + rlx * vsrez[1];
        thi = thi + rlx * vsrez[2];
        dsi = dsi + rlx * vsrez[3];
        uei = uei + rlx * vsrez[4];

        //-------- eliminate absurd transients
        if (ibl >= itran[is]) {
          cti = std::min(cti, 0.30);
          cti = std::max(cti, 0.0000001);
        }

        if (ibl <= iblte[is])
          hklim = 1.02;
        else
          hklim = 1.00005;

        msq = uei * uei * hstinv / (gm1bl * (1.0 - 0.5 * uei * uei * hstinv));
        dsw = dsi - dswaki;
        dslim(dsw, thi, msq, hklim);
        dsi = dsw + dswaki;

        if (dmax <= deps) goto stop110;
      }

      ss << "     mrchdu: convergence failed at " << ibl << " ,  side " << is
         << ", res=" << std::setw(4) << std::fixed << std::setprecision(3)
         << dmax << "\n";
      writeString(ss.str(), true);
      ss.str("");

      if (dmax <= 0.1) goto stop109;
      //------ the current unconverged solution might still be reasonable...

      if (dmax > 0.1) {
        //------- the current solution is garbage --> extrapolate values instead
        if (ibl > 3) {
          if (ibl <= iblte[is]) {
            thi = thet[ibm][is] * sqrt(xssi[ibl][is] / xssi[ibm][is]);
            dsi = dstr[ibm][is] * sqrt(xssi[ibl][is] / xssi[ibm][is]);
            uei = uedg[ibm][is];
          } else {
            if (ibl == iblte[is] + 1) {
              cti = cte;
              thi = tte;
              dsi = dte;
              uei = uedg[ibm][is];
            } else {
              thi = thet[ibm][is];
              ratlen = (xssi[ibl][is] - xssi[ibm][is]) / (10.0 * dstr[ibm][is]);
              dsi = (dstr[ibm][is] + thi * ratlen) / (1.0 + ratlen);
              uei = uedg[ibm][is];
            }
          }
          if (ibl == itran[is]) cti = 0.05;
          if (ibl > itran[is]) cti = ctau[ibm][is];
        }
      }

    stop109:
      blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
      blkin();

      //------- check for transition and set appropriate flags and things
      if ((!simi) && (!turb)) {
        trchek();
        ami = blData2.amplz;
        if (tran) itran[is] = ibl;
        if (!tran) itran[is] = ibl + 2;
      }

      //------- set all other extrapolated values for current station
      if (ibl < itran[is]) blvar(1);
      if (ibl >= itran[is]) blvar(2);
      if (wake) blvar(3);

      if (ibl < itran[is]) blmid(1);
      if (ibl >= itran[is]) blmid(2);
      if (wake) blmid(3);

      //------ pick up here after the newton iterations
    stop110:
      sens = sennew;

      //------ store primary variables
      if (ibl < itran[is])
        ctau[ibl][is] = ami;
      else
        ctau[ibl][is] = cti;
      thet[ibl][is] = thi;
      dstr[ibl][is] = dsi;
      uedg[ibl][is] = uei;
      mass[ibl][is] = dsi * uei;
      ctq[ibl][is] = blData2.cqz;

      //------ set "1" variables to "2" variables for next streamwise station
      blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
      blkin();

      stepbl();

      //------ turbulent intervals will follow transition interval or te
      if (tran || ibl == iblte[is]) {
        turb = true;

        //------- save transition location
        tforce[is] = trforc;
      }

      tran = false;
      //			qApp->processEvents();
      if (s_bCancel) return false;
    }  // 1000 continue
  }    // 2000 continue
  return true;
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
  int is, ibl, ibm, iw, itbl;
  double msq, ratlen, dsw, hklim;
  double hlmax, htmax, xsi, uei, ucon, tsq, thi, ami, cti, dsi;
  double dswaki;
  double htest, hktest, dummy;
  double cst;
  double cte, dte, tte, dmax, hmax, htarg;
  cte = dte = tte = dmax = hmax = htarg = 0.0;

  //---- shape parameters for separation criteria
  hlmax = 3.8;
  htmax = 2.5;

  for (is = 1; is <= 2; is++) {  // 2000
    ss << "    Side " << is << " ...\n";
    writeString(ss.str());

    //---- set forced transition arc length position
    xifset(is);

    //---- initialize similarity station with thwaites' formula
    //	ibl = 2;
    xsi = xssi[2][is];
    uei = uedg[2][is];

    //      bule = log(uedg(ibl+1,is)/uei) / log(xssi(ibl+1,is)/xsi)
    //      bule = std::max( -.08 , bule )
    bule = 1.0;
    ucon = uei / pow(xsi, bule);
    tsq = 0.45 / (ucon * (5.0 * bule + 1.0) * reybl) * pow(xsi, (1.0 - bule));
    thi = sqrt(tsq);
    dsi = 2.2 * thi;
    ami = 0.0;

    //---- initialize ctau for first turbulent station
    cti = 0.03;

    tran = false;
    turb = false;
    itran[is] = iblte[is];

    //---- march downstream
    for (ibl = 2; ibl <= nbl[is]; ibl++) {  // 1000
      ibm = ibl - 1;
      iw = ibl - iblte[is];
      simi = (ibl == 2);
      wake = ibl > iblte[is];

      //------ prescribed quantities
      xsi = xssi[ibl][is];
      uei = uedg[ibl][is];

      if (wake) {
        iw = ibl - iblte[is];
        dswaki = wgap[iw];
      } else
        dswaki = 0.0;

      direct = true;

      //------ newton iteration loop for current station
      for (itbl = 1; itbl <= 25; itbl++) {  // 100

        //-------- assemble 10x3 linearized system for dctau, dth, dds, due, dxi
        //         at the previous "1" station and the current "2" station
        //         (the "1" station coefficients will be ignored)

        blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
        blkin();

        //-------- check for transition and set appropriate flags and things
        if ((!simi) && (!turb)) {
          trchek();
          ami = blData2.amplz;

          //--------- fixed bug   md 7 jun 99
          if (tran) {
            itran[is] = ibl;
            if (cti <= 0.0) {
              cti = 0.03;
              blData2.sz = cti;
            }
          } else
            itran[is] = ibl + 2;
        }

        if (ibl == iblte[is] + 1) {
          tte = thet[iblte[1]][1] + thet[iblte[2]][2];
          dte = dstr[iblte[1]][1] + dstr[iblte[2]][2] + ante;
          cte = (ctau[iblte[1]][1] * thet[iblte[1]][1] +
                 ctau[iblte[2]][2] * thet[iblte[2]][2]) /
                tte;
          tesys(cte, tte, dte);
        } else
          blsys();

        if (direct) {
          //--------- try direct mode (set due = 0 in currently empty 4th line)
          vs2[4][1] = 0.0;
          vs2[4][2] = 0.0;
          vs2[4][3] = 0.0;
          vs2[4][4] = 1.0;
          vsrez[4] = 0.0;
          //--------- solve newton system for current "2" station
          Gauss(4, vs2, vsrez);
          //--------- determine max changes and underrelax if necessary
          dmax = std::max(fabs(vsrez[2] / thi), fabs(vsrez[3] / dsi));
          if (ibl < itran[is]) dmax = std::max(dmax, fabs(vsrez[1] / 10.0));
          if (ibl >= itran[is]) dmax = std::max(dmax, fabs(vsrez[1] / cti));

          rlx = 1.0;
          if (dmax > 0.3) rlx = 0.3 / dmax;
          //--------- see if direct mode is not applicable
          if (ibl != iblte[is] + 1) {
            //---------- calculate resulting kinematic shape parameter hk
            msq =
                uei * uei * hstinv / (gm1bl * (1.0 - 0.5 * uei * uei * hstinv));
            htest = (dsi + rlx * vsrez[3]) / (thi + rlx * vsrez[2]);
            hkin(htest, msq, hktest, dummy, dummy);

            //---------- decide whether to do direct or inverse problem based on
            // hk
            if (ibl < itran[is]) hmax = hlmax;
            if (ibl >= itran[is]) hmax = htmax;
            direct = (hktest < hmax);
          }
          if (direct) {
            //---------- update as usual
            if (ibl >= itran[is]) cti = cti + rlx * vsrez[1];
            thi = thi + rlx * vsrez[2];
            dsi = dsi + rlx * vsrez[3];
          } else {
            //---------- set prescribed hk for inverse calculation at the
            // current station
            if (ibl < itran[is])
              //----------- laminar case: relatively slow increase in hk
              // downstream
              htarg = blData1.hkz + 0.03 * (blData2.xz - blData1.xz) / blData1.tz;
            else if (ibl == itran[is]) {
              //----------- transition interval: weighted laminar and turbulent
              // case
              htarg = blData1.hkz + (0.03 * (xt - blData1.xz) - 0.15 * (blData2.xz - xt)) / blData1.tz;
            } else if (wake) {
              //----------- turbulent wake case:
              //--          asymptotic wake behavior with approximate backward
              // euler
              cst = 0.03 * (blData2.xz - blData1.xz) / blData1.tz;
              blData2.hkz = blData1.hkz;
              blData2.hkz = blData2.hkz - (blData2.hkz + cst * (blData2.hkz - 1.0) * (blData2.hkz - 1.0) * (blData2.hkz - 1.0) -
                           blData1.hkz) /
                              (1.0 + 3.0 * cst * (blData2.hkz - 1.0) * (blData2.hkz - 1.0));
              blData2.hkz = blData2.hkz - (blData2.hkz + cst * (blData2.hkz - 1.0) * (blData2.hkz - 1.0) * (blData2.hkz - 1.0) -
                           blData1.hkz) /
                              (1.0 + 3.0 * cst * (blData2.hkz - 1.0) * (blData2.hkz - 1.0));
              blData2.hkz = blData2.hkz - (blData2.hkz + cst * (blData2.hkz - 1.0) * (blData2.hkz - 1.0) * (blData2.hkz - 1.0) -
                           blData1.hkz) /
                              (1.0 + 3.0 * cst * (blData2.hkz - 1.0) * (blData2.hkz - 1.0));
              htarg = blData2.hkz;
            } else
              htarg =
                  blData1.hkz - 0.15 * (blData2.xz - blData1.xz) /
                            blData1.tz;  //----------- turbulent case: relatively
                                     // fast decrease in hk downstream

            //---------- limit specified hk to something reasonable
            if (wake)
              htarg = std::max(htarg, 1.01);
            else
              htarg = std::max(htarg, hmax);

            ss.str("");
            ss << "     mrchue: inverse mode at " << ibl
               << "    hk=" << std::fixed << setprecision(3) << htarg << "\n";
            writeString(ss.str());
            ss.str("");

            //---------- try again with prescribed hk

            goto stop100;
          }
        } else {
          //-------- inverse mode (force hk to prescribed value htarg)
          vs2[4][1] = 0.0;
          vs2[4][2] = blData2.hkz_tz;
          vs2[4][3] = blData2.hkz_dz;
          vs2[4][4] = blData2.hkz_uz;
          vsrez[4] = htarg - blData2.hkz;
          Gauss(4, vs2, vsrez);

          dmax = std::max(fabs(vsrez[2] / thi), fabs(vsrez[3] / dsi));
          if (ibl >= itran[is]) dmax = std::max(dmax, fabs(vsrez[1] / cti));
          rlx = 1.0;
          if (dmax > 0.3) rlx = 0.3 / dmax;
          //--------- update variables
          if (ibl >= itran[is]) cti = cti + rlx * vsrez[1];
          thi = thi + rlx * vsrez[2];
          dsi = dsi + rlx * vsrez[3];
          uei = uei + rlx * vsrez[4];
        }
        //-------- eliminate absurd transients

        if (ibl >= itran[is]) {
          cti = std::min(cti, 0.30);
          cti = std::max(cti, 0.0000001);
        }
        if (ibl <= iblte[is])
          hklim = 1.02;
        else
          hklim = 1.00005;
        msq = uei * uei * hstinv / (gm1bl * (1.0 - 0.5 * uei * uei * hstinv));
        dsw = dsi - dswaki;
        dslim(dsw, thi, msq, hklim);
        dsi = dsw + dswaki;
        if (dmax <= 0.00001) goto stop110;
      stop100:
        int nothing = 1;

      }  // end itbl loop

      ss.str("");
      ss << "     mrchue: convergence failed at " << ibl << ",  side " << is
         << ", res =" << std::fixed << std::setprecision(3) << dmax << "\n";
      writeString(ss.str(), true);

      //------ the current unconverged solution might still be reasonable...
      if (dmax > 0.1) {
        //------- the current solution is garbage --> extrapolate values instead
        if (ibl > 3) {
          if (ibl <= iblte[is]) {
            thi = thet[ibm][is] * sqrt(xssi[ibl][is] / xssi[ibm][is]);
            dsi = dstr[ibm][is] * sqrt(xssi[ibl][is] / xssi[ibm][is]);
          } else {
            if (ibl == iblte[is] + 1) {
              cti = cte;
              thi = tte;
              dsi = dte;
            } else {
              thi = thet[ibm][is];
              ratlen = (xssi[ibl][is] - xssi[ibm][is]) / (10.0 * dstr[ibm][is]);
              dsi = (dstr[ibm][is] + thi * ratlen) / (1.0 + ratlen);
            }
          }
          if (ibl == itran[is]) cti = 0.05;
          if (ibl > itran[is]) cti = ctau[ibm][is];

          uei = uedg[ibl][is];

          if (ibl < nbl[is])
            uei = 0.5 * (uedg[ibl - 1][is] + uedg[ibl + 1][is]);
        }
      }
      // 109
      blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
      blkin();
      //------- check for transition and set appropriate flags and things
      if ((!simi) && (!turb)) {
        trchek();
        ami = blData2.amplz;
        if (tran) itran[is] = ibl;
        if (!tran) itran[is] = ibl + 2;
      }
      //------- set all other extrapolated values for current station
      if (ibl < itran[is]) blvar(1);
      if (ibl >= itran[is]) blvar(2);
      if (wake) blvar(3);
      if (ibl < itran[is]) blmid(1);
      if (ibl >= itran[is]) blmid(2);
      if (wake) blmid(3);
      //------ pick up here after the newton iterations
    stop110:
      //------ store primary variables
      if (ibl < itran[is]) ctau[ibl][is] = ami;
      if (ibl >= itran[is]) ctau[ibl][is] = cti;
      thet[ibl][is] = thi;
      dstr[ibl][is] = dsi;
      uedg[ibl][is] = uei;
      mass[ibl][is] = dsi * uei;
      ctq[ibl][is] = blData2.cqz;

      //------ set "1" variables to "2" variables for next streamwise station
      blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
      blkin();

      stepbl();

      //------ turbulent intervals will follow transition interval or te
      if (tran || ibl == iblte[is]) {
        turb = true;

        //------- save transition location
        tforce[is] = trforc;
      }

      tran = false;

      if (ibl == iblte[is]) {
        thi = thet[iblte[1]][1] + thet[iblte[2]][2];
        dsi = dstr[iblte[1]][1] + dstr[iblte[2]][2] + ante;
      }
    }  // 1000 continue : end ibl loop
  }    // 2000 continue : end is loop
  return true;
}

/** -------------------------------------------
 *      sets actual mach, reynolds numbers
 *      from unit-cl values and specified cls
 *      depending on matyp,retyp flags.
 * -------------------------------------------- */
bool XFoil::mrcl(double cls, double &m_cls, double &r_cls) {
  std::stringstream ss;
  double rrat, cla;
  cla = std::max(cls, 0.000001);

  switch (mach_type) {
    case MachType::CONSTANT: {
      minf = minf1;
      m_cls = 0.0;
      break;
    }
    case MachType::FIXED_LIFT: {
      minf = minf1 / sqrt(cla);
      m_cls = -0.5 * minf / cla;
      break;
    }
    case MachType::FIXED_LIFT_AND_DYNAMIC_PRESSURE: {
      minf = minf1;
      m_cls = 0.0;
      break;
    }
  }

  switch (reynolds_type) {
    case ReynoldsType::CONSTANT: {
      reinf = reinf1;
      r_cls = 0.0;
      break;
    }
    case ReynoldsType::FIXED_LIFT: {
      reinf = reinf1 / sqrt(cla);
      r_cls = -0.5 * reinf / cla;
      break;
    }
    case ReynoldsType::FIXED_LIFT_AND_DYNAMIC_PRESSURE: {
      reinf = reinf1 / cla;
      r_cls = -reinf / cla;
      break;
    }
  }

  if (minf >= 0.99) {
    // TRACE("      artificially limiting mach to  0.99\n");
    writeString("mrcl: Cl too low for chosen Mach(Cl) dependence\n", true);
    writeString("      artificially limiting mach to  0.99", true);
    minf = 0.99;
    m_cls = 0.0;
  }

  rrat = 1.0;
  if (reinf1 > 0.0) rrat = reinf / reinf1;

  if (rrat > 100.0) {
    // TRACE("     artificially limiting re to %f\n",reinf1*100.0);
    writeString("mrcl: cl too low for chosen Re(Cl) dependence\n", true);
    ss << "      artificially limiting Re to " << std::fixed
       << std::setprecision(0) << reinf1 * 100.0 << "\n";
    writeString(ss.str(), true);
    ss.str("");
    reinf = reinf1 * 100.0;
    r_cls = 0.0;
  }
  return true;
}

bool XFoil::ncalc(MatrixX2d points, VectorXd spline_length, int n, double xn[], double yn[]) {
  double sx, sy, smod;

  if (n <= 1) return false;
  xn = spline::splind(points.col(0), spline_length, n).data();
  yn = spline::splind(points.col(1), spline_length, n).data();
  for (int i = 1; i <= n; i++) {
    sx = yn[i];
    sy = -xn[i];
    smod = sqrt(sx * sx + sy * sy);
    xn[i] = sx / smod;
    yn[i] = sy / smod;
  }

  //---- average normal vectors at corner points
  for (int i = 1; i <= n - 1; i++) {
    if (spline_length[i] == spline_length[i + 1]) {
      sx = 0.5 * (xn[i] + xn[i + 1]);
      sy = 0.5 * (yn[i] + yn[i + 1]);
      smod = sqrt(sx * sx + sy * sy);
      xn[i] = sx / smod;
      yn[i] = sy / smod;
      xn[i + 1] = sx / smod;
      yn[i + 1] = sy / smod;
    }
  }

  return true;
}

/** -----------------------------------------------------------------------
 *	   Calculates current streamfunction psi at panel node or wake node
 *	   i due to freestream and all bound vorticity gam on the airfoil.
 *	   Sensitivities of psi with respect to alpha (z_alfa) and inverse
 *	   qspec dofs (z_qdof0,z_qdof1) which influence gam in inverse cases.
 *	   Also calculates the sensitivity vector dpsi/dgam (dzdg).
 *
 *	   If siglin=true, then psi includes the effects of the viscous
 *	   source distribution sig and the sensitivity vector dpsi/dsig
 *	   (dzdm) is calculated.
 *
 *	   If geolin=true, then the geometric sensitivity vector dpsi/dn
 *	   is calculated, where n is the normal motion of the jth node.
 *
 *			airfoil:  1   < i < n
 *			wake:	  n+1 < i < n+nw
 * ----------------------------------------------------------------------- */
bool XFoil::psilin(int iNode, double xi, double yi, double nxi, double nyi,
                   double &psi, double &psi_ni, bool geolin, bool siglin) {
  int io, jo, jm, jq, jp;

  double dxinv, psum, qtanm, scs, sds, rx1, rx2, sx, sy, dsio, dso, dsm, dsim;
  double sgn, x0, logr0, theta0, rs0, rs1, rs2, nxo, nyo,
      nxp, nyp, ry1, ry2;
  double ssum, sdif, psni, pdni, psx0, psx1, psx2, pdx0, pdx1, pdx2, psyy, pdyy,
      psis, psig, psid;
  double psigx1, psigx2, psigyy, pgamx1, pgamx2, pgamyy, psigni, pgamni;
  double gsum, gdif, gsum1, gsum2, gdif1, gdif2, pdif, dsp, dsip;
  double sigte1, sigte2, gamte1, gamte2, pgam;
  double apan, yy, logr12, logr22, x1i, x2i, yyi, x1o, x1p, x2o, x2p, yyo, yyp;
  double seps;

  //---- distance tolerance for determining if two points are the same
  seps = (spline_length[n] - spline_length[1]) * 0.00001;

  apan = yy = logr12 = logr22 = x1i = x2i = yyi = x1o = x1p = x2o = x2p = yyo =
      yyp = 0.0;

  io = iNode;

  cosa = cos(alfa);
  sina = sin(alfa);

  jp = 0;

  for (jo = 1; jo <= n; jo++) {
    dzdg[jo] = 0.0;
    dzdn[jo] = 0.0;
    dqdg[jo] = 0.0;
  }

  for (jo = 1; jo <= n; jo++) {
    dzdm[jo] = 0.0;
    dqdm[jo] = 0.0;
  }

  psi = 0.0;
  psi_ni = 0.0;

  qtan1 = 0.0;
  qtan2 = 0.0;
  qtanm = 0.0;

  if (sharp) {
    scs = 1.0;
    sds = 0.0;
  } else {
    scs = ante / dste;
    sds = aste / dste;
  }
  for (jo = 1; jo <= n; jo++) {
    // stop10
    jp = jo + 1;
    jm = jo - 1;
    jq = jp + 1;

    if (jo == 1)
      jm = jo;
    else {
      if (jo == n - 1)
        jq = jp;
      else {
        if (jo == n) {
          jp = 1;
          if ((points.row(jo).x() - points.row(jp).x()) * (points.row(jo).x() - points.row(jp).x()) +
                  (points.row(jo).y() - points.row(jp).y()) * (points.row(jo).y() - points.row(jp).y()) <
              seps * seps)
            goto stop12;
        }
      }
    }

    dso = sqrt((points.row(jo).x() - points.row(jp).x()) * (points.row(jo).x() - points.row(jp).x()) +
               (points.row(jo).y() - points.row(jp).y()) * (points.row(jo).y() - points.row(jp).y()));

    //------ skip null panel
    if (fabs(dso) < 1.0e-7) goto stop10;

    dsio = 1.0 / dso;

    apan = apanel[jo];

    rx1 = xi - points.row(jo).x();
    ry1 = yi - points.row(jo).y();
    rx2 = xi - points.row(jp).x();
    ry2 = yi - points.row(jp).y();

    sx = (points.row(jp).x() - points.row(jo).x()) / dso;
    sy = (points.row(jp).y() - points.row(jo).y()) / dso;

    blData1.xz = sx * rx1 + sy * ry1;
    blData2.xz = sx * rx2 + sy * ry2;
    yy = sx * ry1 - sy * rx1;

    rs1 = rx1 * rx1 + ry1 * ry1;
    rs2 = rx2 * rx2 + ry2 * ry2;

    //------ set reflection flag sgn to avoid branch problems with arctan
    if (io >= 1 && io <= n) {
      //------- no problem on airfoil surface
      sgn = 1.0;
    } else {
      //------- make sure arctan falls between  -/+  pi/2
      sgn = sign(1.0, yy);
    }

    //------ set log(r^2) and arctan(x/y), correcting for reflection if any
    if (io != jo && rs1 > 0.0) {
      logr12 = log(rs1);
      blData1.tz = atan2(sgn * blData1.xz, sgn * yy) + (0.5 - 0.5 * sgn) * PI;
    } else {
      logr12 = 0.0;
      blData1.tz = 0.0;
    }

    if (io != jp && rs2 > 0.0) {
      logr22 = log(rs2);
      blData2.tz = atan2(sgn * blData2.xz, sgn * yy) + (0.5 - 0.5 * sgn) * PI;
    } else {
      logr22 = 0.0;
      blData2.tz = 0.0;
    }

    x1i = sx * nxi + sy * nyi;
    x2i = sx * nxi + sy * nyi;
    yyi = sx * nyi - sy * nxi;

    if (geolin) {
      nxo = nx[jo];
      nyo = ny[jo];
      nxp = nx[jp];
      nyp = ny[jp];

      x1o = -((rx1 - blData1.xz * sx) * nxo + (ry1 - blData1.xz * sy) * nyo) / dso -
            (sx * nxo + sy * nyo);
      x1p = ((rx1 - blData1.xz * sx) * nxp + (ry1 - blData1.xz * sy) * nyp) / dso;
      x2o = -((rx2 - blData2.xz * sx) * nxo + (ry2 - blData2.xz * sy) * nyo) / dso;
      x2p = ((rx2 - blData2.xz * sx) * nxp + (ry2 - blData2.xz * sy) * nyp) / dso -
            (sx * nxp + sy * nyp);
      yyo = ((rx1 + blData1.xz * sy) * nyo - (ry1 - blData1.xz * sx) * nxo) / dso -
            (sx * nyo - sy * nxo);
      yyp = -((rx1 - blData1.xz * sy) * nyp - (ry1 + blData1.xz * sx) * nxp) / dso;
    }

    if (jo == n) goto stop11;

    if (siglin) {
      //------- set up midpoint quantities
      x0 = 0.5 * (blData1.xz + blData2.xz);
      rs0 = x0 * x0 + yy * yy;
      logr0 = log(rs0);
      theta0 = atan2(sgn * x0, sgn * yy) + (0.5 - 0.5 * sgn) * PI;

      //------- calculate source contribution to psi	for  1-0  half-panel
      dxinv = 1.0 / (blData1.xz - x0);
      psum = x0 * (theta0 - apan) - blData1.xz * (blData1.tz - apan) +
             0.5 * yy * (logr12 - logr0);
      pdif = ((blData1.xz + x0) * psum + rs1 * (blData1.tz - apan) - rs0 * (theta0 - apan) +
              (x0 - blData1.xz) * yy) *
             dxinv;

      psx1 = -(blData1.tz - apan);
      psx0 = theta0 - apan;
      psyy = 0.5 * (logr12 - logr0);

      pdx1 =
          ((blData1.xz + x0) * psx1 + psum + 2.0 * blData1.xz * (blData1.tz - apan) - pdif) * dxinv;
      pdx0 =
          ((blData1.xz + x0) * psx0 + psum - 2.0 * x0 * (theta0 - apan) + pdif) * dxinv;
      pdyy =
          ((blData1.xz + x0) * psyy + 2.0 * (x0 - blData1.xz + yy * (blData1.tz - theta0))) * dxinv;

      dsm = sqrt((points.row(jp).x() - points.row(jm).x()) * (points.row(jp).x() - points.row(jm).x()) +
                 (points.row(jp).y() - points.row(jm).y()) * (points.row(jp).y() - points.row(jm).y()));
      dsim = 1.0 / dsm;

      ssum = (sig[jp] - sig[jo]) / dso + (sig[jp] - sig[jm]) * dsim;
      sdif = (sig[jp] - sig[jo]) / dso - (sig[jp] - sig[jm]) * dsim;

      psi += qopi * (psum * ssum + pdif * sdif);

      //------- dpsi/dm
      dzdm[jm] += qopi * (-psum * dsim + pdif * dsim);
      dzdm[jo] += qopi * (-psum / dso - pdif / dso);
      dzdm[jp] += qopi * (psum * (dsio + dsim) + pdif * (dsio - dsim));

      //------- dpsi/dni
      psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
      pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;
      psi_ni = psi_ni + qopi * (psni * ssum + pdni * sdif);

      qtanm = qtanm + qopi * (psni * ssum + pdni * sdif);

      dqdm[jm] += qopi * (-psni * dsim + pdni * dsim);
      dqdm[jo] += qopi * (-psni / dso - pdni / dso);
      dqdm[jp] += qopi * (psni * (dsio + dsim) + pdni * (dsio - dsim));

      //------- calculate source contribution to psi	for  0-2  half-panel
      dxinv = 1.0 / (x0 - blData2.xz);
      psum = blData2.xz * (blData2.tz - apan) - x0 * (theta0 - apan) +
             0.5 * yy * (logr0 - logr22);
      pdif = ((x0 + blData2.xz) * psum + rs0 * (theta0 - apan) - rs2 * (blData2.tz - apan) +
              (blData2.xz - x0) * yy) *
             dxinv;

      psx0 = -(theta0 - apan);
      psx2 = blData2.tz - apan;
      psyy = 0.5 * (logr0 - logr22);

      pdx0 =
          ((x0 + blData2.xz) * psx0 + psum + 2.0 * x0 * (theta0 - apan) - pdif) * dxinv;
      pdx2 =
          ((x0 + blData2.xz) * psx2 + psum - 2.0 * blData2.xz * (blData2.tz - apan) + pdif) * dxinv;
      pdyy =
          ((x0 + blData2.xz) * psyy + 2.0 * (blData2.xz - x0 + yy * (theta0 - blData2.tz))) * dxinv;

      dsp = sqrt((points.row(jq).x() - points.row(jo).x()) * (points.row(jq).x() - points.row(jo).x()) +
                 (points.row(jq).y() - points.row(jo).y()) * (points.row(jq).y() - points.row(jo).y()));
      dsip = 1.0 / dsp;

      ssum = (sig[jq] - sig[jo]) * dsip + (sig[jp] - sig[jo]) / dso;
      sdif = (sig[jq] - sig[jo]) * dsip - (sig[jp] - sig[jo]) / dso;

      psi = psi + qopi * (psum * ssum + pdif * sdif);

      //------- dpsi/dm
      dzdm[jo] += qopi * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
      dzdm[jp] += qopi * (psum / dso - pdif / dso);
      dzdm[jq] += qopi * (psum * dsip + pdif * dsip);

      //------- dpsi/dni
      psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
      pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;
      psi_ni = psi_ni + qopi * (psni * ssum + pdni * sdif);

      qtanm = qtanm + qopi * (psni * ssum + pdni * sdif);

      dqdm[jo] += qopi * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
      dqdm[jp] += qopi * (psni / dso - pdni / dso);
      dqdm[jq] += qopi * (psni * dsip + pdni * dsip);
    }

    //------ calculate vortex panel contribution to psi
    dxinv = 1.0 / (blData1.xz - blData2.xz);
    psis = 0.5 * blData1.xz * logr12 - 0.5 * blData2.xz * logr22 + blData2.xz - blData1.xz +
           yy * (blData1.tz - blData2.tz);
    psid = ((blData1.xz + blData2.xz) * psis +
            0.5 * (rs2 * logr22 - rs1 * logr12 + blData1.xz * blData1.xz - blData2.xz * blData2.xz)) *
           dxinv;

    psx1 = 0.5 * logr12;
    psx2 = -.5 * logr22;
    psyy = blData1.tz - blData2.tz;

    pdx1 = ((blData1.xz + blData2.xz) * psx1 + psis - blData1.xz * logr12 - psid) * dxinv;
    pdx2 = ((blData1.xz + blData2.xz) * psx2 + psis + blData2.xz * logr22 + psid) * dxinv;
    pdyy = ((blData1.xz + blData2.xz) * psyy - yy * (logr12 - logr22)) * dxinv;

    gsum1 = gamu[jp][1] + gamu[jo][1];
    gsum2 = gamu[jp][2] + gamu[jo][2];
    gdif1 = gamu[jp][1] - gamu[jo][1];
    gdif2 = gamu[jp][2] - gamu[jo][2];

    gsum = gam[jp] + gam[jo];
    gdif = gam[jp] - gam[jo];

    psi += qopi * (psis * gsum + psid * gdif);

    //------ dpsi/dgam
    dzdg[jo] += qopi * (psis - psid);
    dzdg[jp] += qopi * (psis + psid);

    //------ dpsi/dni
    psni = psx1 * x1i + psx2 * x2i + psyy * yyi;
    pdni = pdx1 * x1i + pdx2 * x2i + pdyy * yyi;
    psi_ni += qopi * (gsum * psni + gdif * pdni);

    qtan1 += qopi * (gsum1 * psni + gdif1 * pdni);
    qtan2 += qopi * (gsum2 * psni + gdif2 * pdni);

    dqdg[jo] += qopi * (psni - pdni);
    dqdg[jp] += qopi * (psni + pdni);

    if (geolin) {
      //------- dpsi/dn
      dzdn[jo] += qopi * gsum * (psx1 * x1o + psx2 * x2o + psyy * yyo) +
                  qopi * gdif * (pdx1 * x1o + pdx2 * x2o + pdyy * yyo);
      dzdn[jp] += qopi * gsum * (psx1 * x1p + psx2 * x2p + psyy * yyp) +
                  qopi * gdif * (pdx1 * x1p + pdx2 * x2p + pdyy * yyp);
    }
  stop10:
    int nothing;
    nothing = 1;
  }

stop11:
  psig = 0.5 * yy * (logr12 - logr22) + blData2.xz * (blData2.tz - apan) -
         blData1.xz * (blData1.tz - apan);
  pgam =
      0.5 * blData1.xz * logr12 - 0.5 * blData2.xz * logr22 + blData2.xz - blData1.xz + yy * (blData1.tz - blData2.tz);

  psigx1 = -(blData1.tz - apan);
  psigx2 = blData2.tz - apan;
  psigyy = 0.5 * (logr12 - logr22);
  pgamx1 = 0.5 * logr12;
  pgamx2 = -0.5 * logr22;
  pgamyy = blData1.tz - blData2.tz;

  psigni = psigx1 * x1i + psigx2 * x2i + psigyy * yyi;
  pgamni = pgamx1 * x1i + pgamx2 * x2i + pgamyy * yyi;

  //---- TE panel source and vortex strengths
  sigte1 = 0.5 * scs * (gamu[jp][1] - gamu[jo][1]);
  sigte2 = 0.5 * scs * (gamu[jp][2] - gamu[jo][2]);
  gamte1 = -0.5 * sds * (gamu[jp][1] - gamu[jo][1]);
  gamte2 = -0.5 * sds * (gamu[jp][2] - gamu[jo][2]);

  sigte = 0.5 * scs * (gam[jp] - gam[jo]);
  gamte = -0.5 * sds * (gam[jp] - gam[jo]);

  //---- TE panel contribution to psi
  psi += hopi * (psig * sigte + pgam * gamte);

  //---- dpsi/dgam
  dzdg[jo] += -hopi * psig * scs * 0.5;
  dzdg[jp] += +hopi * psig * scs * 0.5;

  dzdg[jo] += +hopi * pgam * sds * 0.5;
  dzdg[jp] += -hopi * pgam * sds * 0.5;

  //---- dpsi/dni
  psi_ni += hopi * (psigni * sigte + pgamni * gamte);

  qtan1 += hopi * (psigni * sigte1 + pgamni * gamte1);
  qtan2 += hopi * (psigni * sigte2 + pgamni * gamte2);

  dqdg[jo] += -hopi * (psigni * 0.5 * scs - pgamni * 0.5 * sds);
  dqdg[jp] += +hopi * (psigni * 0.5 * scs - pgamni * 0.5 * sds);

  if (geolin) {
    //----- dpsi/dn
    dzdn[jo] += hopi * (psigx1 * x1o + psigx2 * x2o + psigyy * yyo) * sigte +
                hopi * (pgamx1 * x1o + pgamx2 * x2o + pgamyy * yyo) * gamte;
    dzdn[jp] += hopi * (psigx1 * x1p + psigx2 * x2p + psigyy * yyp) * sigte +
                hopi * (pgamx1 * x1p + pgamx2 * x2p + pgamyy * yyp) * gamte;

  }
stop12:

  //**** freestream terms
  psi += qinf * (cosa * yi - sina * xi);

  //---- dpsi/dn
  psi_ni = psi_ni + qinf * (cosa * nyi - sina * nxi);

  qtan1 += qinf * nyi;
  qtan2 += -qinf * nxi;

  // techwinder: removed image calculattion
  return false;
}

/** --------------------------------------------------------------------
 *	   Calculates current streamfunction psi and tangential velocity
 *	   qtan at panel node or wake node i due to freestream and wake
 *	   sources sig.  also calculates sensitivity vectors dpsi/dsig
 *	   (dzdm) and dqtan/dsig (dqdm).
 *
 *			airfoil:  1   < i < n
 *			wake:	  n+1 < i < n+nw
 *-------------------------------------------------------------------- */
bool XFoil::pswlin(int i, double xi, double yi, double nxi, double nyi,
                   double &psi, double &psi_ni) {
  double g1, g2, t1, t2;
  int io, jo;

  io = i;

  cosa = cos(alfa);
  sina = sin(alfa);

  for (jo = n + 1; jo <= n + nw; jo++) {
    dzdm[jo] = 0.0;
    dqdm[jo] = 0.0;
  }

  psi = 0.0;
  psi_ni = 0.0;

  for (jo = n + 1; jo <= n + nw - 1; jo++) {
    int jp = jo + 1;
    int jm = jo - 1;
    int jq = jp + 1;
    if (jo == n + 1) {
      jm = jo;
    } else {
      if (jo == n + nw - 1) jq = jp;
    }
    const double dso = sqrt((points.row(jo).x() - points.row(jp).x()) * (points.row(jo).x() - points.row(jp).x()) +
               (points.row(jo).y() - points.row(jp).y()) * (points.row(jo).y() - points.row(jp).y()));
    const double dsio = 1.0 / dso;

    const double apan = apanel[jo];

    const double rx1 = xi - points.row(jo).x();
    const double ry1 = yi - points.row(jo).y();
    const double rx2 = xi - points.row(jp).x();
    const double ry2 = yi - points.row(jp).y();

    const double sx = (points.row(jp).x() - points.row(jo).x()) * dsio;
    const double sy = (points.row(jp).y() - points.row(jo).y()) * dsio;

    blData1.xz = sx * rx1 + sy * ry1;
    blData2.xz = sx * rx2 + sy * ry2;
    const double yy = sx * ry1 - sy * rx1;
    const double rs1 = rx1 * rx1 + ry1 * ry1;
    const double rs2 = rx2 * rx2 + ry2 * ry2;

    double sgn = 1.0;

    if (io >= n + 1 && io <= n + nw) {
      sgn = 1.0;
    } else {
      sgn = sign(1.0, yy);
    }

    if (io != jo && rs1 > 0.0) {
      g1 = log(rs1);
      t1 = atan2(sgn * blData1.xz, sgn * yy) - (0.5 - 0.5 * sgn) * PI;
    } else {
      g1 = 0.0;
      t1 = 0.0;
    }

    if (io != jp && rs2 > 0.0) {
      g2 = log(rs2);
      t2 = atan2(sgn * blData2.xz, sgn * yy) - (0.5 - 0.5 * sgn) * PI;
    } else {
      g2 = 0.0;
      t2 = 0.0;
    }
    const double x1i = sx * nxi + sy * nyi;
    const double x2i = sx * nxi + sy * nyi;
    const double yyi = sx * nyi - sy * nxi;
    //------- set up midpoint quantities
    const double x0 = 0.5 * (blData1.xz + blData2.xz);
    const double rs0 = x0 * x0 + yy * yy;
    const double g0 = log(rs0);
    const double t0 = atan2(sgn * x0, sgn * yy) - (0.5 - 0.5 * sgn) * PI;

    //------- calculate source contribution to psi	for  1-0  half-panel
    double dxinv = 1.0 / (blData1.xz - x0);
    double psum = x0 * (t0 - apan) - blData1.xz * (t1 - apan) + 0.5 * yy * (g1 - g0);
    double pdif = ((blData1.xz + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan) +
            (x0 - blData1.xz) * yy) *
           dxinv;

    double psx1 = -(t1 - apan);
    double psx0 = t0 - apan;
    double psyy = 0.5 * (g1 - g0);

    double pdx1 = ((blData1.xz + x0) * psx1 + psum + 2.0 * blData1.xz * (t1 - apan) - pdif) * dxinv;
    double pdx0 = ((blData1.xz + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
    double pdyy = ((blData1.xz + x0) * psyy + 2.0 * (x0 - blData1.xz + yy * (t1 - t0))) * dxinv;

    const double dsm = sqrt((points.row(jp).x() - points.row(jm).x()) * (points.row(jp).x() - points.row(jm).x()) +
               (points.row(jp).y() - points.row(jm).y()) * (points.row(jp).y() - points.row(jm).y()));
    const double dsim = 1.0 / dsm;

    double ssum = (sig[jp] - sig[jo]) / dso + (sig[jp] - sig[jm]) * dsim;
    double sdif = (sig[jp] - sig[jo]) / dso - (sig[jp] - sig[jm]) * dsim;

    psi = psi + qopi * (psum * ssum + pdif * sdif);

    //------- dpsi/dm
    dzdm[jm] = dzdm[jm] + qopi * (-psum * dsim + pdif * dsim);
    dzdm[jo] = dzdm[jo] + qopi * (-psum / dso - pdif / dso);
    dzdm[jp] = dzdm[jp] + qopi * (psum * (dsio + dsim) + pdif * (dsio - dsim));

    //------- dpsi/dni
    double psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
    double pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;
    psi_ni = psi_ni + qopi * (psni * ssum + pdni * sdif);

    dqdm[jm] = dqdm[jm] + qopi * (-psni * dsim + pdni * dsim);
    dqdm[jo] = dqdm[jo] + qopi * (-psni / dso - pdni / dso);
    dqdm[jp] = dqdm[jp] + qopi * (psni * (dsio + dsim) + pdni * (dsio - dsim));

    //------- calculate source contribution to psi	for  0-2  half-panel
    dxinv = 1.0 / (x0 - blData2.xz);
    psum = blData2.xz * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
    pdif = ((x0 + blData2.xz) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan) +
            (blData2.xz - x0) * yy) *
           dxinv;

    psx0 = -(t0 - apan);
    const double psx2 = t2 - apan;
    psyy = 0.5 * (g0 - g2);

    pdx0 = ((x0 + blData2.xz) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv;
    const double pdx2 = ((x0 + blData2.xz) * psx2 + psum - 2.0 * blData2.xz * (t2 - apan) + pdif) * dxinv;
    pdyy = ((x0 + blData2.xz) * psyy + 2.0 * (blData2.xz - x0 + yy * (t0 - t2))) * dxinv;

    const double dsp = sqrt((points.row(jq).x() - points.row(jo).x()) * (points.row(jq).x() - points.row(jo).x()) +
               (points.row(jq).y() - points.row(jo).y()) * (points.row(jq).y() - points.row(jo).y()));
    const double dsip = 1.0 / dsp;

    ssum = (sig[jq] - sig[jo]) * dsip + (sig[jp] - sig[jo]) / dso;
    sdif = (sig[jq] - sig[jo]) * dsip - (sig[jp] - sig[jo]) / dso;

    psi = psi + qopi * (psum * ssum + pdif * sdif);

    //------- dpsi/dm
    dzdm[jo] = dzdm[jo] + qopi * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
    dzdm[jp] = dzdm[jp] + qopi * (psum / dso - pdif / dso);
    dzdm[jq] = dzdm[jq] + qopi * (psum * dsip + pdif * dsip);

    //------- dpsi/dni
    psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
    pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;
    psi_ni = psi_ni + qopi * (psni * ssum + pdni * sdif);

    dqdm[jo] = dqdm[jo] + qopi * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
    dqdm[jp] = dqdm[jp] + qopi * (psni / dso - pdni / dso);
    dqdm[jq] = dqdm[jq] + qopi * (psni * dsip + pdni * dsip);
  }

  return true;
}

/** -----------------------------------------------------
 * 	   calculates source panel influence coefficient
 * 	   matrix for current airfoil and wake geometry.
 * ------------------------------------------------------ */
bool XFoil::qdcalc() {
  int i, j, k, iu;
  double psi, psi_n, sum;
  double bbb[IQX];

  // TRACE("calculating source influence matrix ...\n");
  writeString("   Calculating source influence matrix ...\n");

  if (!ladij) {
    //----- calculate source influence matrix for airfoil surface if it doesn't
    // exist
    for (j = 1; j <= n; j++) {
      //------- multiply each dpsi/sig vector by inverse of factored dpsi/dgam
      // matrix
      for (iu = 0; iu < IQX; iu++)
        bbb[iu] = bij[iu][j];  // techwinder : create a dummy array
      baksub(n + 1, aij, aijpiv, bbb);
      for (iu = 0; iu < IQX; iu++) bij[iu][j] = bbb[iu];

      //------- store resulting dgam/dsig = dqtan/dsig vector
      for (i = 1; i <= n; i++) {
        dij[i][j] = bij[i][j];
      }
    }
    ladij = true;
  }

  //---- set up coefficient matrix of dpsi/dm on airfoil surface
  for (i = 1; i <= n; i++) {
    pswlin(i, points.row(i).x(), points.row(i).y(), nx[i], ny[i], psi, psi_n);
    for (j = n + 1; j <= n + nw; j++) {
      bij[i][j] = -dzdm[j];
    }
  }

  //---- set up kutta condition (no direct source influence)
  for (j = n + 1; j <= n + nw; j++) bij[n + 1][j] = 0.0;

  //---- sharp te gamma extrapolation also has no source influence
  if (sharp) {
    for (j = n + 1; j <= n + nw; j++) bij[n][j] = 0.0;
  }

  //---- multiply by inverse of factored dpsi/dgam matrix
  for (j = n + 1; j <= n + nw; j++) {
    //		baksub(iqx,n+1,aijpiv,j);
    for (iu = 0; iu < IQX; iu++)
      bbb[iu] = bij[iu][j];  // techwinder : create a dummy array

    baksub(n + 1, aij, aijpiv, bbb);
    for (iu = 0; iu < IQX; iu++) bij[iu][j] = bbb[iu];
  }
  //---- set the source influence matrix for the wake sources
  for (i = 1; i <= n; i++) {
    for (j = n + 1; j <= n + nw; j++) {
      dij[i][j] = bij[i][j];
    }
  }

  //**** now we need to calculate the influence of sources on the wake
  // velocities

  //---- calculate dqtan/dgam and dqtan/dsig at the wake points

  for (i = n + 1; i <= n + nw; i++) {
    int iw = i - n;
    //------ airfoil contribution at wake panel node
    psilin(i, points.row(i).x(), points.row(i).y(), nx[i], ny[i], psi, psi_n, false, true);
    for (j = 1; j <= n; j++) {
      cij[iw][j] = dqdg[j];
    }
    for (j = 1; j <= n; j++) {
      dij[i][j] = dqdm[j];
    }
    //------ wake contribution
    pswlin(i, points.row(i).x(), points.row(i).y(), nx[i], ny[i], psi, psi_n);
    for (j = n + 1; j <= n + nw; j++) {
      dij[i][j] = dqdm[j];
    }
  }

  //---- add on effect of all sources on airfoil vorticity which effects wake
  // qtan
  for (i = n + 1; i <= n + nw; i++) {
    int iw = i - n;

    //------ airfoil surface source contribution first
    for (j = 1; j <= n; j++) {
      sum = 0.0;
      for (k = 1; k <= n; k++) sum = sum + cij[iw][k] * dij[k][j];
      dij[i][j] = dij[i][j] + sum;
    }

    //------ wake source contribution next
    for (j = n + 1; j <= n + nw; j++) {
      sum = 0.0;
      for (k = 1; k <= n; k++) sum = sum + cij[iw][k] * bij[k][j];
      dij[i][j] = dij[i][j] + sum;
    }
  }

  //---- make sure first wake point has same velocity as trailing edge
  for (j = 1; j <= n + nw; j++) {
    dij[n + 1][j] = dij[n][j];
  }

  lwdij = true;
  return true;
}

/** -------------------------------------------------------
 *     sets inviscid panel tangential velocity for
 *      current alpha.
 * -------------------------------------------------------- */
bool XFoil::qiset() {
  cosa = cos(alfa);
  sina = sin(alfa);

  for (int i = 1; i <= n + nw; i++) {
    qinv[i] = cosa * qinvu[i][1] + sina * qinvu[i][2];
    qinv_a[i] = -sina * qinvu[i][1] + cosa * qinvu[i][2];
  }

  return true;
}

/** -------------------------------------------------------------
 *     sets panel viscous tangential velocity from viscous ue
 * -------------------------------------------------------------- */
bool XFoil::qvfue() {
  int is, ibl;
  for (is = 1; is <= 2; is++) {
    for (ibl = 2; ibl <= nbl[is]; ibl++) {
      int i = ipan[ibl][is];
      qvis[i] = vti[ibl][is] * uedg[ibl][is];
    }
  }

  return true;
}

/** ---------------------------------------------------------------
 *      sets inviscid tangential velocity for alpha = 0, 90
 *      on wake due to freestream and airfoil surface vorticity.
 * --------------------------------------------------------------- */
bool XFoil::qwcalc() {
  double psi, psi_ni;
  int i;

  //---- first wake point (same as te)
  qinvu[n + 1][1] = qinvu[n][1];
  qinvu[n + 1][2] = qinvu[n][2];

  //---- rest of wake
  for (i = n + 2; i <= n + nw; i++) {
    psilin(i, points.row(i).x(), points.row(i).y(), nx[i], ny[i], psi, psi_ni, false, false);
    qinvu[i][1] = qtan1;
    qinvu[i][2] = qtan2;
  }

  return true;
}

bool XFoil::restoreblData(int icom) {
  if (icom == 1) {
    blData1 = blsav[icom];
  } else if (icom == 2) {
    blData2 = blsav[icom];
  }
  return true;
}

bool XFoil::saveblData(int icom) {
  if (icom == 1) {
    blsav[icom] = blData1;
  } else {
    blsav[icom] = blData2;
  }
  return true;
}

bool XFoil::setbl() {
  //-------------------------------------------------
  //	   sets up the bl newton system coefficients for the current bl
  // variables
  //     and the edge velocities received from setup. the local bl system
  //     coefficients are then incorporated into the global newton system.
  //-------------------------------------------------

  std::stringstream ss;
  int i, ibl, iv, iw, j, js = 0, jv, jbl, is = 0;
  int ile1 = 0, ile2 = 0, ite1 = 0, ite2 = 0, jvte1 = 0, jvte2 = 0;
  double usav[IVX + 1][ISX];
  double u1_m[2 * IVX + 1], u2_m[2 * IVX + 1];
  double d1_m[2 * IVX + 1], d2_m[2 * IVX + 1];
  double ule1_m[2 * IVX + 1], ule2_m[2 * IVX + 1];
  double ute1_m[2 * IVX + 1], ute2_m[2 * IVX + 1];

  for (int i = 0; i < IVX + 1; i++) memset(usav[i], 0, ISX * sizeof(double));
  memset(u1_m, 0, (2 * IVX + 1) * sizeof(double));
  memset(u2_m, 0, (2 * IVX + 1) * sizeof(double));
  memset(d1_m, 0, (2 * IVX + 1) * sizeof(double));
  memset(d2_m, 0, (2 * IVX + 1) * sizeof(double));
  memset(ule1_m, 0, (2 * IVX + 1) * sizeof(double));
  memset(ule2_m, 0, (2 * IVX + 1) * sizeof(double));
  memset(ute1_m, 0, (2 * IVX + 1) * sizeof(double));
  memset(ute2_m, 0, (2 * IVX + 1) * sizeof(double));

  double msq_clmr = 0.0, mdi;
  double herat = 0.0, herat_ms = 0.0;

  double clmr = 0.0, ma_clmr = 0.0, re_clmr = 0.0;
  double ule1_a = 0.0, ule2_a = 0.0, u2_a, due2, dds2;
  double xsi, cti = 0.0, uei, thi, dsi, dswaki;
  double d2_a, d2_m2, d2_u2, dte_mte1, dte_ute1, dte_mte2, dte_ute2;
  double tte, cte, dte, dule1 = 0.0, dule2 = 0.0;
  double str, chx, chy, xtr, ytr, chsq;
  double xi_ule1, xi_ule2;
  double ami = 0.0, tte_tte1 = 0.0, tte_tte2 = 0.0, cte_tte1 = 0.0,
         cte_tte2 = 0.0, cte_cte1 = 0.0, cte_cte2 = 0.0;

  //---- set the cl used to define mach, reynolds numbers
  if (lalfa)
    clmr = cl;
  else
    clmr = clspec;

  cti = 0.0;  // techwinder added, otherwise variable is not initialized

  //---- set current minf(cl)
  mrcl(clmr, ma_clmr, re_clmr);
  msq_clmr = 2.0 * minf * ma_clmr;

  //---- set compressibility parameter tklam and derivative tk_msq
  comset();

  //---- set gas constant (= cp/cv)
  gambl = gamma;
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

  //---- sutherland's const./to	(assumes stagnation conditions are at stp)
  hvrat = 0.35;

  //---- set reynolds number based on freestream density, velocity, viscosity
  herat = 1.0 - 0.5 * qinfbl * qinfbl * hstinv;
  herat_ms = -0.5 * qinfbl * qinfbl * hstinv_ms;

  reybl = reinf * sqrt(herat * herat * herat) * (1.0 + hvrat) / (herat + hvrat);
  reybl_re = sqrt(herat * herat * herat) * (1.0 + hvrat) / (herat + hvrat);
  reybl_ms = reybl * (1.5 / herat - 1.0 / (herat + hvrat)) * herat_ms;

  amcrit = acrit;

  //---- save te thickness
  dwte = wgap[1];

  if (!lblini) {
    //----- initialize bl by marching with ue (fudge at separation)
    // TRACE(" initializing bl ...\n");
    writeString("   Initializing bl ...\n");

    mrchue();
    lblini = true;
  }

  //---- march bl with current ue and ds to establish transition
  mrchdu();

  for (is = 1; is <= 2; is++) {
    for (ibl = 2; ibl <= nbl[is]; ibl++) usav[ibl][is] = uedg[ibl][is];
  }

  ueset();

  for (is = 1; is <= 2; is++) {
    for (ibl = 2; ibl <= nbl[is]; ibl++) {
      double temp = usav[ibl][is];
      usav[ibl][is] = uedg[ibl][is];
      uedg[ibl][is] = temp;
    }
  }
  ile1 = ipan[2][1];
  ile2 = ipan[2][2];
  ite1 = ipan[iblte[1]][1];
  ite2 = ipan[iblte[2]][2];

  jvte1 = isys[iblte[1]][1];
  jvte2 = isys[iblte[2]][2];

  dule1 = uedg[2][1] - usav[2][1];
  dule2 = uedg[2][2] - usav[2][2];

  //---- set le and te ue sensitivities wrt all m values
  for (js = 1; js <= 2; js++) {
    for (jbl = 2; jbl <= nbl[js]; jbl++) {
      j = ipan[jbl][js];
      jv = isys[jbl][js];
      ule1_m[jv] = -vti[2][1] * vti[jbl][js] * dij[ile1][j];
      ule2_m[jv] = -vti[2][2] * vti[jbl][js] * dij[ile2][j];
      ute1_m[jv] = -vti[iblte[1]][1] * vti[jbl][js] * dij[ite1][j];
      ute2_m[jv] = -vti[iblte[2]][2] * vti[jbl][js] * dij[ite2][j];
    }
  }

  ule1_a = uinv_a[2][1];
  ule2_a = uinv_a[2][2];

  writeString(" \n");

  //*** go over each boundary layer/wake
  for (is = 1; is <= 2; is++) {
    //---- there is no station "1" at similarity, so zero everything out
    for (js = 1; js <= 2; js++) {
      for (jbl = 2; jbl <= nbl[js]; jbl++) {
        jv = isys[jbl][js];
        u1_m[jv] = 0.0;
        d1_m[jv] = 0.0;
      }
    }
    double u1_a = 0.0;
    double d1_a = 0.0;

    double due1 = 0.0;
    double dds1 = 0.0;

    //---- similarity station pressure gradient parameter  x/u du/dx
    ibl = 2;
    bule = 1.0;

    //---- set forced transition arc length position
    xifset(is);

    tran = false;
    turb = false;

    //**** sweep downstream setting up bl equation linearizations
    for (ibl = 2; ibl <= nbl[is]; ibl++) {
      iv = isys[ibl][is];

      simi = (ibl == 2);
      wake = (ibl > iblte[is]);
      tran = (ibl == itran[is]);
      turb = (ibl > itran[is]);

      i = ipan[ibl][is];

      //---- set primary variables for current station
      xsi = xssi[ibl][is];
      if (ibl < itran[is])
        ami = ctau[ibl][is];
      else
        cti = ctau[ibl][is];
      uei = uedg[ibl][is];
      thi = thet[ibl][is];
      mdi = mass[ibl][is];

      dsi = mdi / uei;

      if (wake) {
        iw = ibl - iblte[is];
        dswaki = wgap[iw];
      } else
        dswaki = 0.0;

      //---- set derivatives of dsi (= d2)
      d2_m2 = 1.0 / uei;
      d2_u2 = -dsi / uei;

      for (js = 1; js <= 2; js++) {
        for (jbl = 2; jbl <= nbl[js]; jbl++) {
          j = ipan[jbl][js];
          jv = isys[jbl][js];
          u2_m[jv] = -vti[ibl][is] * vti[jbl][js] * dij[i][j];
          d2_m[jv] = d2_u2 * u2_m[jv];
        }
      }
      d2_m[iv] = d2_m[iv] + d2_m2;

      u2_a = uinv_a[ibl][is];
      d2_a = d2_u2 * u2_a;

      //---- "forced" changes due to mismatch between uedg and
      // usav=uinv+dij*mass
      due2 = uedg[ibl][is] - usav[ibl][is];
      dds2 = d2_u2 * due2;

      blprv(xsi, ami, cti, thi, dsi, dswaki, uei);  // cti
      blkin();

      //---- check for transition and set tran, xt, etc. if found
      if (tran) {
        trchek();
        ami = blData2.amplz;
      }

      if (ibl == itran[is] && !tran) {
        // TRACE("setbl: xtr???  n1=%d n2=%d: \n", ampl1, ampl2);

        ss << "setbl: xtr???  n1=" << blData1.amplz << " n2=" << blData2.amplz << ":\n";
        writeString(ss.str());
        ss.str("");
      }

      //---- assemble 10x4 linearized system for dctau, dth, dds, due, dxi
      //	   at the previous "1" station and the current "2" station

      if (ibl == iblte[is] + 1) {
        //----- define quantities at start of wake, adding te base thickness to
        // dstar
        tte = thet[iblte[1]][1] + thet[iblte[2]][2];
        dte = dstr[iblte[1]][1] + dstr[iblte[2]][2] + ante;
        cte = (ctau[iblte[1]][1] * thet[iblte[1]][1] +
               ctau[iblte[2]][2] * thet[iblte[2]][2]) /
              tte;
        tesys(cte, tte, dte);

        tte_tte1 = 1.0;
        tte_tte2 = 1.0;
        dte_mte1 = 1.0 / uedg[iblte[1]][1];
        dte_ute1 = -dstr[iblte[1]][1] / uedg[iblte[1]][1];
        dte_mte2 = 1.0 / uedg[iblte[2]][2];
        dte_ute2 = -dstr[iblte[2]][2] / uedg[iblte[2]][2];
        cte_cte1 = thet[iblte[1]][1] / tte;
        cte_cte2 = thet[iblte[2]][2] / tte;
        cte_tte1 = (ctau[iblte[1]][1] - cte) / tte;
        cte_tte2 = (ctau[iblte[2]][2] - cte) / tte;

        //----- re-define d1 sensitivities wrt m since d1 depends on both te ds
        // values
        for (js = 1; js <= 2; js++) {
          for (jbl = 2; jbl <= nbl[js]; jbl++) {            
            jv = isys[jbl][js];
            d1_m[jv] = dte_ute1 * ute1_m[jv] + dte_ute2 * ute2_m[jv];
          }
        }
        d1_m[jvte1] = d1_m[jvte1] + dte_mte1;
        d1_m[jvte2] = d1_m[jvte2] + dte_mte2;

        //----- "forced" changes from  uedg --- usav=uinv+dij*mass	mismatch
        due1 = 0.0;
        dds1 = dte_ute1 * (uedg[iblte[1]][1] - usav[iblte[1]][1]) +
               dte_ute2 * (uedg[iblte[2]][2] - usav[iblte[2]][2]);
      } else {
        blsys();
      }

      //---- save wall shear and equil. max shear coefficient for plotting
      // output
      ctq[ibl][is] = blData2.cqz;

      //---- set xi sensitivities wrt le ue changes
      if (is == 1) {
        xi_ule1 = sst_go;
        xi_ule2 = -sst_gp;
      } else {
        xi_ule1 = -sst_go;
        xi_ule2 = sst_gp;
      }

      //---- stuff bl system coefficients into main jacobian matrix

      for (jv = 1; jv <= nsys; jv++) {
        vm[1][jv][iv] = vs1[1][3] * d1_m[jv] + vs1[1][4] * u1_m[jv] +
                        vs2[1][3] * d2_m[jv] + vs2[1][4] * u2_m[jv] +
                        (vs1[1][5] + vs2[1][5] + vsx[1]) *
                            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
      }

      vb[1][1][iv] = vs1[1][1];
      vb[1][2][iv] = vs1[1][2];

      va[1][1][iv] = vs2[1][1];
      va[1][2][iv] = vs2[1][2];

      if (lalfa)
        vdel[1][2][iv] = vsr[1] * re_clmr + vsm[1] * msq_clmr;
      else
        vdel[1][2][iv] = (vs1[1][4] * u1_a + vs1[1][3] * d1_a) +
                         (vs2[1][4] * u2_a + vs2[1][3] * d2_a) +
                         (vs1[1][5] + vs2[1][5] + vsx[1]) *
                             (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

      vdel[1][1][iv] = vsrez[1] + (vs1[1][4] * due1 + vs1[1][3] * dds1) +
                       (vs2[1][4] * due2 + vs2[1][3] * dds2) +
                       (vs1[1][5] + vs2[1][5] + vsx[1]) *
                           (xi_ule1 * dule1 + xi_ule2 * dule2);

      for (jv = 1; jv <= nsys; jv++) {
        vm[2][jv][iv] = vs1[2][3] * d1_m[jv] + vs1[2][4] * u1_m[jv] +
                        vs2[2][3] * d2_m[jv] + vs2[2][4] * u2_m[jv] +
                        (vs1[2][5] + vs2[2][5] + vsx[2]) *
                            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
      }
      vb[2][1][iv] = vs1[2][1];
      vb[2][2][iv] = vs1[2][2];

      va[2][1][iv] = vs2[2][1];
      va[2][2][iv] = vs2[2][2];

      if (lalfa)
        vdel[2][2][iv] = vsr[2] * re_clmr + vsm[2] * msq_clmr;
      else
        vdel[2][2][iv] = (vs1[2][4] * u1_a + vs1[2][3] * d1_a) +
                         (vs2[2][4] * u2_a + vs2[2][3] * d2_a) +
                         (vs1[2][5] + vs2[2][5] + vsx[2]) *
                             (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

      vdel[2][1][iv] = vsrez[2] + (vs1[2][4] * due1 + vs1[2][3] * dds1) +
                       (vs2[2][4] * due2 + vs2[2][3] * dds2) +
                       (vs1[2][5] + vs2[2][5] + vsx[2]) *
                           (xi_ule1 * dule1 + xi_ule2 * dule2);

      // memory overlap problem
      for (jv = 1; jv <= nsys; jv++) {
        vm[3][jv][iv] = vs1[3][3] * d1_m[jv] + vs1[3][4] * u1_m[jv] +
                        vs2[3][3] * d2_m[jv] + vs2[3][4] * u2_m[jv] +
                        (vs1[3][5] + vs2[3][5] + vsx[3]) *
                            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
      }

      vb[3][1][iv] = vs1[3][1];
      vb[3][2][iv] = vs1[3][2];

      va[3][1][iv] = vs2[3][1];
      va[3][2][iv] = vs2[3][2];

      if (lalfa)
        vdel[3][2][iv] = vsr[3] * re_clmr + vsm[3] * msq_clmr;
      else
        vdel[3][2][iv] = (vs1[3][4] * u1_a + vs1[3][3] * d1_a) +
                         (vs2[3][4] * u2_a + vs2[3][3] * d2_a) +
                         (vs1[3][5] + vs2[3][5] + vsx[3]) *
                             (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

      vdel[3][1][iv] = vsrez[3] + (vs1[3][4] * due1 + vs1[3][3] * dds1) +
                       (vs2[3][4] * due2 + vs2[3][3] * dds2) +
                       (vs1[3][5] + vs2[3][5] + vsx[3]) *
                           (xi_ule1 * dule1 + xi_ule2 * dule2);

      if (ibl == iblte[is] + 1) {
        //----- redefine coefficients for tte, dte, etc
        vz[1][1] = vs1[1][1] * cte_cte1;
        vz[1][2] = vs1[1][1] * cte_tte1 + vs1[1][2] * tte_tte1;
        vb[1][1][iv] = vs1[1][1] * cte_cte2;
        vb[1][2][iv] = vs1[1][1] * cte_tte2 + vs1[1][2] * tte_tte2;

        vz[2][1] = vs1[2][1] * cte_cte1;
        vz[2][2] = vs1[2][1] * cte_tte1 + vs1[2][2] * tte_tte1;
        vb[2][1][iv] = vs1[2][1] * cte_cte2;
        vb[2][2][iv] = vs1[2][1] * cte_tte2 + vs1[2][2] * tte_tte2;

        vz[3][1] = vs1[3][1] * cte_cte1;
        vz[3][2] = vs1[3][1] * cte_tte1 + vs1[3][2] * tte_tte1;
        vb[3][1][iv] = vs1[3][1] * cte_cte2;
        vb[3][2][iv] = vs1[3][1] * cte_tte2 + vs1[3][2] * tte_tte2;
      }

      //---- turbulent intervals will follow if currently at transition interval
      if (tran) {
        turb = true;

        //------ save transition location
        itran[is] = ibl;
        tforce[is] = trforc;

        //------ interpolate airfoil geometry to find transition x/c
        //		(for user output)
        if (is == 1)
          str = sst - xt;
        else
          str = sst + xt;

        chx = point_te.x() - point_le.x();
        chy = point_te.y() - point_le.y();
        chsq = chx * chx + chy * chy;
        xtr = spline::seval(str, points.col(0), dpoints_ds.col(0), spline_length, n);
        ytr = spline::seval(str, points.col(1), dpoints_ds.col(1), spline_length, n);
        xoctr[is] = ((xtr - point_le.x()) * chx + (ytr - point_le.y()) * chy) / chsq;
        yoctr[is] = ((ytr - point_le.y()) * chx - (xtr - point_le.x()) * chy) / chsq;
      }

      tran = false;

      if (ibl == iblte[is]) {
        //----- set "2" variables at te to wake correlations for next station

        turb = true;
        wake = true;
        blvar(3);
        blmid(3);
      }

      for (js = 1; js <= 2; js++) {
        for (jbl = 2; jbl <= nbl[js]; jbl++) {
          jv = isys[jbl][js];
          u1_m[jv] = u2_m[jv];
          d1_m[jv] = d2_m[jv];
        }
      }

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

    if (tforce[is]) {
      ss << "     Side " << is << ", forced transition at x/c = " << std::fixed
         << std::setprecision(4) << xoctr[is] << " " << itran[is] << "\n";
      writeString(ss.str());
      ss.str("");
    } else {
      ss << "     Side " << is << ",  free  transition at x/c = " << std::fixed
         << std::setprecision(4) << xoctr[is] << " " << itran[is] << "\n";
      writeString(ss.str());
      ss.str("");
    }

    //---- next airfoil side
  }

  return true;
}

bool XFoil::setexp(double spline_length[], double ds1, double smax, int nn) {
  //........................................................
  //     sets geometriy stretched array s:
  //
  //       s(i+1) - s(i)  =  r * [s(i) - s(i-1)]
  //
  //       s     (output)  array to be set
  //       ds1   (input)   first s increment:  spline_length[2] - spline_length[1]
  //       smax  (input)   final s value:      s(nn)
  //       nn    (input)   number of points
  //........................................................
  int nex, iter, n;
  double sigma, rnex, rni, aaa, bbb, ccc;
  double disc, ratio, sigman, res;
  double dresdr, dratio, ds;

  sigma = smax / ds1;
  nex = nn - 1;
  rnex = (double)nex;
  rni = 1.0 / rnex;

  //-- solve quadratic for initial geometric ratio guess
  aaa = rnex * (rnex - 1.0) * (rnex - 2.0) / 6.0;
  bbb = rnex * (rnex - 1.0) / 2.0;
  ccc = rnex - sigma;

  disc = bbb * bbb - 4.0 * aaa * ccc;
  disc = std::max(0.0, disc);

  if (nex <= 1) {
    writeString("setexp: cannot fill array.  n too small\n", true);
    return false;
  } else {
    if (nex == 2)
      ratio = -ccc / bbb + 1.0;
    else
      ratio = (-bbb + sqrt(disc)) / (2.0 * aaa) + 1.0;
  }
  if (ratio == 1.0) goto stop11;

  //-- newton iteration for actual geometric ratio
  for (iter = 1; iter <= 100; iter++) {
    sigman = (pow(ratio, (double)nex) - 1.0) / (ratio - 1.0);
    res = pow(sigman, rni) - pow(sigma, rni);
    dresdr = rni * pow(sigman, rni) *
             (rnex * pow(ratio, (double)(nex - 1)) - sigman) /
             (pow(ratio, (double)nex) - 1.0);

    dratio = -res / dresdr;
    ratio = ratio + dratio;

    if (fabs(dratio) < 1.0e-5) goto stop11;
  }

  writeString("Setexp: Convergence failed.  Continuing anyway ...\n", true);

  //-- set up stretched array using converged geometric ratio
stop11:
  spline_length[1] = 0.0;
  ds = ds1;
  for (n = 2; n <= nn; n++) {
    spline_length[n] = spline_length[n - 1] + ds;
    ds = ds * ratio;
  }
  return true;
}

bool XFoil::setMach() {
  mrcl(1.0, minf_cl, reinf_cl);
  comset();
  cpcalc(n, qinv, qinf, minf, cpi);
  if (lvisc) {
    cpcalc(n + nw, qvis, qinf, minf, cpv);
  }
  clcalc(xcmref, ycmref);
  cdcalc();
  lvconv = false;
  return true;
}

/** returns the absolute value of "a" x sign(b) */
double XFoil::sign(double a, double b) {
  if (b >= 0.0)
    return fabs(a);
  else
    return -fabs(a);
}

/**
 * 	   Calculates the "inverse" spline function s(x).
 * 	   Since s(x) can be multi-valued or not defined,
 * 	   this is not a "black-box" routine.  The calling
 * 	   program must pass via si a sufficiently good
 * 	   initial guess for s(xi).
 *
 * 	   xi	   specified x value	   (input)
 * 	   si	   calculated s(xi) value  (input,output)
 * 	   x,xs,s  usual spline arrays	   (input)
 */
bool XFoil::sinvrt(double &si, double xi, VectorXd x, VectorXd xs, VectorXd spline_length, int n) {
  int iter;
  double sisav;
  sisav = si;

  for (iter = 1; iter <= 10; iter++) {
    const double res = spline::seval(si, x, xs, spline_length, n) - xi;
    const double resp = spline::deval(si, x, xs, spline_length, n);
    const double ds = -res / resp;
    si = si + ds;
    if (fabs(ds / (spline_length[n] - spline_length[1])) < 1.0e-5) return true;
  }

  writeString("Sinvrt: spline inversion failed, input value returned\n", true);
  si = sisav;

  return false;
}

/**
 *      Converges to specified alpha.
 */
bool XFoil::specal() {
  double minf_clm, reinf_clm;
  double clm;
  int i, irlx, itcl;

  //---- calculate surface vorticity distributions for alpha = 0, 90 degrees
  if (!lgamu || !lqaij) ggcalc();

  cosa = cos(alfa);
  sina = sin(alfa);

  //---- superimpose suitably weighted  alpha = 0, 90  distributions
  for (i = 1; i <= n; i++) {
    gam[i] = cosa * gamu[i][1] + sina * gamu[i][2];
    gam_a[i] = -sina * gamu[i][1] + cosa * gamu[i][2];
  }
  psio = cosa * gamu[n + 1][1] + sina * gamu[n + 1][2];

  tecalc();
  qiset();

  //---- set initial guess for the newton variable clm
  clm = 1.0;

  //---- set corresponding  m(clm), re(clm)
  mrcl(clm, minf_clm, reinf_clm);
  comset();

  //---- set corresponding cl(m)
  clcalc(xcmref, ycmref);
  //---- iterate on clm
  bool bConv = false;
  for (itcl = 1; itcl <= 20; itcl++) {
    const double msq_clm = 2.0 * minf * minf_clm;
    const double dclm = (cl - clm) / (1.0 - cl_msq * msq_clm);

    const double clm1 = clm;
    rlx = 1.0;

    //------ under-relaxation loop to avoid driving m(cl) above 1
    for (irlx = 1; irlx <= 12; irlx++) {
      clm = clm1 + rlx * dclm;

      //-------- set new freestream mach m(clm)
      mrcl(clm, minf_clm, reinf_clm);

      //-------- if mach is ok, go do next newton iteration
      if (mach_type == MachType::CONSTANT || minf == 0.0 || minf_clm != 0.0) break;  // goto 91

      rlx = 0.5 * rlx;
    }

    //------ set new cl(m)
    comset();
    clcalc(xcmref, ycmref);

    if (fabs(dclm) <= 1.0e-6) {
      bConv = true;
      break;
    }
  }
  if (!bConv) {
    writeString("Specal:  MInf convergence failed\n", true);
    return false;
  }

  //---- set final mach, cl, cp distributions, and hinge moment
  mrcl(cl, minf_cl, reinf_cl);
  comset();
  clcalc(xcmref, ycmref);
  /*
          if (!cpcalc(n,qinv,qinf,minf,cpi)) return false;// no need to carry on
          if(lvisc) {
                  if(!cpcalc(n+nw,qvis,qinf,minf,cpv)) return false;// no need
     to carry on if(!cpcalc(n+nw,qinv,qinf,minf,cpi)) return false;// no need to
     carry on
          }
          else   if (!cpcalc(n,qinv,qinf,minf,cpi)) return false;// no need to
     carry on
  */
  cpcalc(n, qinv, qinf, minf, cpi);
  if (lvisc) {
    cpcalc(n + nw, qvis, qinf, minf, cpv);
    cpcalc(n + nw, qinv, qinf, minf, cpi);
  } else
    cpcalc(n, qinv, qinf, minf, cpi);

  // Added techwinder to get inviscid q after viscous calculation
  for (i = 1; i <= n; i++) {
    qgamm[i] = gam[i];
  }
  // end techwinder addition

  return true;
}

bool XFoil::speccl() {
  //-----------------------------------------
  //     converges to specified inviscid cl.
  //-----------------------------------------
  int i, ital;

  //---- calculate surface vorticity distributions for alpha = 0, 90 degrees
  if (!lgamu || !lqaij) ggcalc();

  //---- set freestream mach from specified cl -- mach will be held fixed
  mrcl(clspec, minf_cl, reinf_cl);
  comset();

  //---- current alpha is the initial guess for newton variable alfa
  cosa = cos(alfa);
  sina = sin(alfa);

  for (i = 1; i <= n; i++) {
    gam[i] = cosa * gamu[i][1] + sina * gamu[i][2];
    gam_a[i] = -sina * gamu[i][1] + cosa * gamu[i][2];
  }
  psio = cosa * gamu[n + 1][1] + sina * gamu[n + 1][2];

  //---- get corresponding cl, cl_alpha, cl_mach
  clcalc(xcmref, ycmref);

  //---- newton loop for alpha to get specified inviscid cl
  bool bConv = false;
  for (ital = 1; ital <= 20; ital++) {
    const double dalfa = (clspec - cl) / cl_alf;
    rlx = 1.0;

    alfa = alfa + rlx * dalfa;

    //------ set new surface speed distribution
    cosa = cos(alfa);
    sina = sin(alfa);
    for (i = 1; i <= n; i++) {
      gam[i] = cosa * gamu[i][1] + sina * gamu[i][2];
      gam_a[i] = -sina * gamu[i][1] + cosa * gamu[i][2];
    }
    psio = cosa * gamu[n + 1][1] + sina * gamu[n + 1][2];

    //------ set new cl(alpha)
    clcalc(xcmref, ycmref);

    if (fabs(dalfa) <= 1.0e-6) {
      bConv = true;
      break;
    }
  }
  if (!bConv) {
    writeString("Speccl:  cl convergence failed", true);
    return false;
  }

  //---- set final surface speed and cp distributions
  tecalc();
  qiset();

  if (lvisc) {
    cpcalc(n + nw, qvis, qinf, minf, cpv);
    cpcalc(n + nw, qinv, qinf, minf, cpi);

  } else {
    cpcalc(n, qinv, qinf, minf, cpi);
  }

  return true;
}

bool XFoil::stepbl() {
    blData1 = blData2;
    return true;
}

bool XFoil::stfind() {
  //-----------------------------------------
  //     locates stagnation point arc length
  //     location sst and panel index ist.
  //-----------------------------------------
  double dgam, ds;
  int i;
  bool bFound;

  bFound = false;
  for (i = 1; i <= n - 1; i++) {
    if (gam[i] >= 0.0 && gam[i + 1] < 0.0) {
      bFound = true;
      break;
    }
  }

  if (!bFound) {
    writeString("stfind: Stagnation point not found. Continuing ...\n", true);
    i = n / 2;
  }

  ist = i;
  dgam = gam[i + 1] - gam[i];
  ds = spline_length[i + 1] - spline_length[i];

  //---- evaluate so as to minimize roundoff for very small gam[i] or gam[i+1]
  if (gam[i] < -gam[i + 1])
    sst = spline_length[i] - ds * (gam[i] / dgam);
  else
    sst = spline_length[i + 1] - ds * (gam[i + 1] / dgam);

  //---- tweak stagnation point if it falls right on a node (very unlikely)
  if (sst <= spline_length[i]) sst = spline_length[i] + 0.0000001;
  if (sst >= spline_length[i + 1]) sst = spline_length[i + 1] - 0.0000001;

  sst_go = (sst - spline_length[i + 1]) / dgam;
  sst_gp = (spline_length[i] - sst) / dgam;

  return true;
}

bool XFoil::stmove() {
  //--------------------------------------------------
  //    moves stagnation point location to new panel.
  //---------------------------------------------------
  int ibl, istold, is;
  //-- locate new stagnation point arc length sst from gam distribution
  istold = ist;
  stfind();

  if (istold == ist) {
    //--- recalculate new arc length array
    xicalc();
  }
  else {
    //--- set new bl position -> panel position  pointers
    iblpan();

    //--- set new inviscid bl edge velocity uinv from qinv
    uicalc();

    //--- recalculate new arc length array
    xicalc();

    //--- set  bl position -> system line  pointers
    iblsys();

    if (ist > istold) {
      //---- increase in number of points on top side (is=1)
      int idif = ist - istold;

      itran[1] = itran[1] + idif;
      itran[2] = itran[2] - idif;

      //---- move top side bl variables downstream
      for (ibl = nbl[1]; ibl >= idif + 2; ibl--) {
        ctau[ibl][1] = ctau[ibl - idif][1];
        thet[ibl][1] = thet[ibl - idif][1];
        dstr[ibl][1] = dstr[ibl - idif][1];
        uedg[ibl][1] = uedg[ibl - idif][1];
      }

      //---- set bl variables between old and new stagnation point
      const double dudx = uedg[idif + 2][1] / xssi[idif + 2][1];
      for (ibl = idif + 1; ibl >= 2; ibl--) {
        ctau[ibl][1] = ctau[idif + 2][1];
        thet[ibl][1] = thet[idif + 2][1];
        dstr[ibl][1] = dstr[idif + 2][1];
        uedg[ibl][1] = dudx * xssi[ibl][1];
      }

      //---- move bottom side bl variables upstream
      for (ibl = 2; ibl <= nbl[2]; ibl++) {
        ctau[ibl][2] = ctau[ibl + idif][2];
        thet[ibl][2] = thet[ibl + idif][2];
        dstr[ibl][2] = dstr[ibl + idif][2];
        uedg[ibl][2] = uedg[ibl + idif][2];
      }
    } else {
      //---- increase in number of points on bottom side (is=2)
      int idif = istold - ist;

      itran[1] = itran[1] - idif;
      itran[2] = itran[2] + idif;

      //---- move bottom side bl variables downstream
      for (ibl = nbl[2]; ibl >= idif + 2; ibl--) {
        ctau[ibl][2] = ctau[ibl - idif][2];
        thet[ibl][2] = thet[ibl - idif][2];
        dstr[ibl][2] = dstr[ibl - idif][2];
        uedg[ibl][2] = uedg[ibl - idif][2];
      }

      //---- set bl variables between old and new stagnation point
      const double dudx = uedg[idif + 2][2] / xssi[idif + 2][2];
      for (ibl = idif + 1; ibl >= 2; ibl--) {
        ctau[ibl][2] = ctau[idif + 2][2];
        thet[ibl][2] = thet[idif + 2][2];
        dstr[ibl][2] = dstr[idif + 2][2];
        uedg[ibl][2] = dudx * xssi[ibl][2];
      }

      //---- move top side bl variables upstream
      for (ibl = 2; ibl <= nbl[1]; ibl++) {
        ctau[ibl][1] = ctau[ibl + idif][1];
        thet[ibl][1] = thet[ibl + idif][1];
        dstr[ibl][1] = dstr[ibl + idif][1];
        uedg[ibl][1] = uedg[ibl + idif][1];
      }
    }
  }

  //-- set new mass array since ue has been tweaked
  for (is = 1; is <= 2; is++) {
    for (ibl = 2; ibl <= nbl[is]; ibl++)
      mass[ibl][is] = dstr[ibl][is] * uedg[ibl][is];
  }

  return true;
}

bool XFoil::tecalc() {
  //-------------------------------------------
  //     calculates total and projected TE
  //     areas and TE panel strengths.
  //-------------------------------------------

  double scs, sds;
  //---- set te base vector and te bisector components
  double dxte = points.row(1).x() - points.row(n).x();
  double dyte = points.row(1).y() - points.row(n).y();
  double dxs = 0.5 * (-dpoints_ds.row(1).x() + dpoints_ds.row(n).x());
  double dys = 0.5 * (-dpoints_ds.row(1).y() + dpoints_ds.row(n).y());

  //---- normal and streamwise projected TE gap areas
  ante = dxs * dyte - dys * dxte;
  aste = dxs * dxte + dys * dyte;

  //---- total TE gap area
  dste = sqrt(dxte * dxte + dyte * dyte);

  sharp = dste < 0.0001 * chord;

  if (sharp) {
    scs = 1.0;
    sds = 0.0;
  } else {
    scs = ante / dste;
    sds = aste / dste;
  }

  //---- TE panel source and vorticity strengths
  sigte = 0.5 * (gam[1] - gam[n]) * scs;
  gamte = -.5 * (gam[1] - gam[n]) * sds;

  return true;
}

bool XFoil::tesys(double cte, double tte, double dte) {
  //--------------------------------------------------------
  //	   sets up "dummy" bl system between airfoil te point
  //	   and first wake point infinitesimally behind te.
  //--------------------------------------------------------

  for (int k = 1; k <= 4; k++) {
    vsrez[k] = 0.0;
    vsm[k] = 0.0;
    vsr[k] = 0.0;
    vsx[k] = 0.0;
    for (int l = 1; l <= 5; l++) {
      vs1[k][l] = 0.0;
      vs2[k][l] = 0.0;
    }
  }

  blvar(3);

  vs1[1][1] = -1.0;
  vs2[1][1] = 1.0;
  vsrez[1] = cte - blData2.sz;

  vs1[2][2] = -1.0;
  vs2[2][2] = 1.0;
  vsrez[2] = tte - blData2.tz;

  vs1[3][3] = -1.0;
  vs2[3][3] = 1.0;
  vsrez[3] = dte - blData2.dz - blData2.dwz;

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
  int itam = 0;
  double ax_hk1 = 0.0, ax_t1 = 0.0, ax_a1 = 0.0, ax_hk2 = 0.0, ax_t2 = 0.0,
         ax_rt2 = 0.0, ax_a2 = 0.0;
  double amplt, sfa, sfa_a1, sfa_a2, sfx;
  double sfx_x1, sfx_x2, sfx_xf;
  double tt, dt, ut, amsave;
  double ax = 0.0, ax_rt1 = 0.0, res = 0.0, res_a2 = 0.0;
  double da2 = 0.0, dxt = 0.0, tt_t1 = 0.0, dt_d1 = 0.0, ut_u1 = 0.0;
  double tt_t2 = 0.0, dt_d2 = 0.0, ut_u2 = 0.0, tt_a1 = 0.0, dt_a1 = 0.0;
  double ut_a1 = 0.0, tt_x1 = 0.0, dt_x1 = 0.0, ut_x1 = 0.0, tt_x2 = 0.0,
         dt_x2 = 0.0, ut_x2 = 0.0;
  double ax_d1 = 0.0, ax_u1 = 0.0, ax_x1 = 0.0, ax_d2 = 0.0, ax_u2 = 0.0,
         ax_x2 = 0.0, ax_ms = 0.0, ax_re = 0.0;
  double z_ax = 0.0, z_a1 = 0.0, z_t1 = 0.0, z_d1 = 0.0, z_u1 = 0.0, z_x1 = 0.0,
         z_a2 = 0.0, z_t2 = 0.0, z_d2 = 0.0, z_u2 = 0.0, z_x2 = 0.0, z_ms = 0.0,
         z_re = 0.0;
  double ax_at, ax_rtt, ax_tt, ax_hkt, amplt_a2, wf1, wf1_a1, wf1_a2, wf1_xf,
      wf1_x1, wf1_x2;
  double wf2, wf2_a1, wf2_a2, wf2_xf, wf2_x1, wf2_x2, xt_a2, dt_a2, tt_a2;
  double ut_a2, hkt, hkt_tt, hkt_dt, hkt_ut, hkt_ms, rtt_tt, rtt_ut, rtt_ms,
      rtt, rtt_re;
  double daeps = 0.00005;

  ax_at = ax_rtt = ax_tt = ax_hkt = amplt_a2 = wf1 = wf1_a1 = wf1_xf = wf1_x1 = wf1_x2 = 0.0;
  wf2 = wf2_a1 = wf2_xf = wf2_x1 = wf2_x2 = xt_a2 = dt_a2 = tt_a2 = 0.0;
  ut_a2 = hkt_tt = hkt_dt = hkt_ut = hkt_ms = rtt_tt = rtt_ut = rtt_ms = rtt_re = 0.0;

  //---- save variables and sensitivities at ibl ("2") for future restoration
  saveblData(2);

  //---- calculate average amplification rate ax over x1..x2 interval
  axset(blData1.hkz, blData1.tz, blData1.rtz, blData1.amplz, blData2.hkz, blData2.tz, blData2.rtz, blData2.amplz, amcrit, ax, ax_hk1,
        ax_t1, ax_rt1, ax_a1, ax_hk2, ax_t2, ax_rt2, ax_a2);

  //---- set initial guess for iterate n2 (ampl2) at x2
  blData2.amplz = blData1.amplz + ax * (blData2.xz - blData1.xz);

  //---- solve implicit system for amplification ampl2
  for (itam = 1; itam <= 30; itam++) {
    //---- define weighting factors wf1,wf2 for defining "t" quantities from 1,2
    if (blData2.amplz <= amcrit) {
      //------ there is no transition yet,  "t" is the same as "2"
      amplt = blData2.amplz;
      amplt_a2 = 1.0;
      sfa = 1.0;
      sfa_a1 = 0.0;
      sfa_a2 = 0.0;
    } else {
      //------ there is transition in x1..x2, "t" is set from n1, n2
      amplt = amcrit;
      amplt_a2 = 0.0;
      sfa = (amplt - blData1.amplz) / (blData2.amplz - blData1.amplz);
      sfa_a1 = (sfa - 1.0) / (blData2.amplz - blData1.amplz);
      sfa_a2 = (-sfa) / (blData2.amplz - blData1.amplz);
    }

    if (xiforc < blData2.xz) {
      sfx = (xiforc - blData1.xz) / (blData2.xz - blData1.xz);
      sfx_x1 = (sfx - 1.0) / (blData2.xz - blData1.xz);
      sfx_x2 = (-sfx) / (blData2.xz - blData1.xz);
      sfx_xf = 1.0 / (blData2.xz - blData1.xz);
    } else {
      sfx = 1.0;
      sfx_x1 = 0.0;
      sfx_x2 = 0.0;
      sfx_xf = 0.0;
    }

    //---- set weighting factor from free or forced transition
    if (sfa < sfx) {
      wf2 = sfa;
      wf2_a1 = sfa_a1;
      wf2_a2 = sfa_a2;
      wf2_x1 = 0.0;
      wf2_x2 = 0.0;
      wf2_xf = 0.0;
    } else {
      wf2 = sfx;
      wf2_a1 = 0.0;
      wf2_a2 = 0.0;
      wf2_x1 = sfx_x1;
      wf2_x2 = sfx_x2;
      wf2_xf = sfx_xf;
    }

    wf1 = 1.0 - wf2;
    wf1_a1 = -wf2_a1;
    wf1_a2 = -wf2_a2;
    wf1_x1 = -wf2_x1;
    wf1_x2 = -wf2_x2;
    wf1_xf = -wf2_xf;

    //---- interpolate bl variables to xt
    xt = blData1.xz * wf1 + blData2.xz * wf2;
    tt = blData1.tz * wf1 + blData2.tz * wf2;
    dt = blData1.dz * wf1 + blData2.dz * wf2;
    ut = blData1.uz * wf1 + blData2.uz * wf2;

    xt_a2 = blData1.xz * wf1_a2 + blData2.xz * wf2_a2;
    tt_a2 = blData1.tz * wf1_a2 + blData2.tz * wf2_a2;
    dt_a2 = blData1.dz * wf1_a2 + blData2.dz * wf2_a2;
    ut_a2 = blData1.uz * wf1_a2 + blData2.uz * wf2_a2;

    //---- temporarily set "2" variables from "t" for blkin
    blData2.xz = xt;
    blData2.tz = tt;
    blData2.dz = dt;
    blData2.uz = ut;

    //---- calculate laminar secondary "t" variables hkt, rtt
    blkin();

    hkt = blData2.hkz;
    hkt_tt = blData2.hkz_tz;
    hkt_dt = blData2.hkz_dz;
    hkt_ut = blData2.hkz_uz;
    hkt_ms = blData2.hkz_ms;

    rtt = blData2.rtz;
    rtt_tt = blData2.rtz_tz;
    rtt_ut = blData2.rtz_uz;
    rtt_ms = blData2.rtz_ms;
    rtt_re = blData2.rtz_re;

    //---- restore clobbered "2" variables, except for ampl2
    amsave = blData2.amplz;

    restoreblData(2);

    blData2.amplz = amsave;

    //---- calculate amplification rate ax over current x1-xt interval
    axset(blData1.hkz, blData1.tz, blData1.rtz, blData1.amplz, hkt, tt, rtt, amplt, amcrit, ax, ax_hk1,
          ax_t1, ax_rt1, ax_a1, ax_hkt, ax_tt, ax_rtt, ax_at);

    //---- punch out early if there is no amplification here
    if (ax <= 0.0) goto stop101;

    //---- set sensitivity of ax(a2)
    ax_a2 = (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_a2 +
            (ax_hkt * hkt_dt) * dt_a2 +
            (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_a2 + ax_at * amplt_a2;

    //---- residual for implicit ampl2 definition (amplification equation)
    res = blData2.amplz - blData1.amplz - ax * (blData2.xz - blData1.xz);
    res_a2 = 1.0 - ax_a2 * (blData2.xz - blData1.xz);

    da2 = -res / res_a2;

    rlx = 1.0;
    dxt = xt_a2 * da2;

    if (rlx * fabs(dxt / (blData2.xz - blData1.xz)) > 0.05) rlx = 0.05 * fabs((blData2.xz - blData1.xz) / dxt);

    if (rlx * fabs(da2) > 1.0) rlx = 1.0 * fabs(1.0 / da2);

    //---- check if converged
    if (fabs(da2) < daeps) goto stop101;

    if ((blData2.amplz > amcrit && blData2.amplz + rlx * da2 < amcrit) ||
        (blData2.amplz < amcrit && blData2.amplz + rlx * da2 > amcrit))
      //------ limited newton step so ampl2 doesn't step across amcrit either
      // way
      blData2.amplz = amcrit;
    else
      //------ regular newton step
      blData2.amplz = blData2.amplz + rlx * da2;
  }

  // TRACE("trchek2 - n2 convergence failed\n");
  writeString("trchek2 - n2 convergence failed\n", true);
  if (s_bCancel) return false;
stop101:

  //---- test for free or forced transition
  trfree = (blData2.amplz >= amcrit);
  trforc = (xiforc > blData1.xz) && (xiforc <= blData2.xz);

  //---- set transition interval flag
  tran = (trforc || trfree);

  if (!tran) return false;

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

  xt_x1 = wf1;
  tt_t1 = wf1;
  dt_d1 = wf1;
  ut_u1 = wf1;

  xt_x2 = wf2;
  tt_t2 = wf2;
  dt_d2 = wf2;
  ut_u2 = wf2;

  xt_a1 = blData1.xz * wf1_a1 + blData2.xz * wf2_a1;
  tt_a1 = blData1.tz * wf1_a1 + blData2.tz * wf2_a1;
  dt_a1 = blData1.dz * wf1_a1 + blData2.dz * wf2_a1;
  ut_a1 = blData1.uz * wf1_a1 + blData2.uz * wf2_a1;

  xt_x1 = blData1.xz * wf1_x1 + blData2.xz * wf2_x1 + xt_x1;
  tt_x1 = blData1.tz * wf1_x1 + blData2.tz * wf2_x1;
  dt_x1 = blData1.dz * wf1_x1 + blData2.dz * wf2_x1;
  ut_x1 = blData1.uz * wf1_x1 + blData2.uz * wf2_x1;

  xt_x2 = blData1.xz * wf1_x2 + blData2.xz * wf2_x2 + xt_x2;
  tt_x2 = blData1.tz * wf1_x2 + blData2.tz * wf2_x2;
  dt_x2 = blData1.dz * wf1_x2 + blData2.dz * wf2_x2;
  ut_x2 = blData1.uz * wf1_x2 + blData2.uz * wf2_x2;

  xt_xf = blData1.xz * wf1_xf + blData2.xz * wf2_xf;

  //---- at this point, ax = ax( hk1, t1, rt1, a1, hkt, tt, rtt, at )

  //---- set sensitivities of ax( t1 d1 u1 a1 t2 d2 u2 a2 ms re )
  ax_t1 = ax_hk1 * blData1.hkz_tz + ax_t1 + ax_rt1 * blData1.rtz_tz +
          (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_t1;
  ax_d1 = ax_hk1 * blData1.hkz_dz + (ax_hkt * hkt_dt) * dt_d1;
  ax_u1 = ax_hk1 * blData1.hkz_uz + ax_rt1 * blData1.rtz_uz +
          (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_u1;
  ax_a1 = ax_a1 +
          (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) *
              tt_a1
          + (ax_hkt * hkt_dt) * dt_a1 +
          (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_a1;
  ax_x1 = (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_x1 +
          (ax_hkt * hkt_dt) * dt_x1 +
          (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_x1;

  ax_t2 = (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_t2;
  ax_d2 = (ax_hkt * hkt_dt) * dt_d2;
  ax_u2 = (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_u2;
  ax_a2 = ax_at * amplt_a2 +
          (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) *
              tt_a2 
          + (ax_hkt * hkt_dt) * dt_a2 +
          (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_a2;
  ax_x2 = (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_x2 +
          (ax_hkt * hkt_dt) * dt_x2 +
          (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_x2;

  ax_ms = ax_hkt * hkt_ms + ax_rtt * rtt_ms + ax_hk1 * blData1.hkz_ms + ax_rt1 * blData1.rtz_ms;
  ax_re = ax_rtt * rtt_re + ax_rt1 * blData1.rtz_re;

  //---- set sensitivities of residual res
  z_ax = -(blData2.xz - blData1.xz);

  z_a1 = z_ax * ax_a1 - 1.0;
  z_t1 = z_ax * ax_t1;
  z_d1 = z_ax * ax_d1;
  z_u1 = z_ax * ax_u1;
  z_x1 = z_ax * ax_x1 + ax;

  z_a2 = z_ax * ax_a2 + 1.0;
  z_t2 = z_ax * ax_t2;
  z_d2 = z_ax * ax_d2;
  z_u2 = z_ax * ax_u2;
  z_x2 = z_ax * ax_x2 - ax;

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
  double bl1[5][6], bl2[5][6], blrez[5], blm[5], blr[5], blx[5], bt1[5][6],
      bt2[5][6], btrez[5], btm[5], btr[5], btx[5];
  double wf2, wf2_xt, wf2_a1, wf2_x1, wf2_x2, wf2_t1, wf2_t2;
  double wf2_d1, wf2_d2, wf2_u1, wf2_u2, wf2_ms, wf2_re, wf2_xf;
  double wf1, wf1_a1, wf1_x1, wf1_x2, wf1_t1, wf1_t2, wf1_d1, wf1_d2;
  double wf1_u1, wf1_u2, wf1_ms, wf1_re, wf1_xf;
  double tt, tt_a1, tt_x1, tt_x2, tt_t1, tt_t2, tt_d1, tt_d2, tt_u1, tt_u2;
  double tt_ms, tt_re, tt_xf, dt, dt_a1, dt_x1, dt_x2, dt_t1, dt_t2;
  double dt_d1, dt_d2, dt_u1, dt_u2, dt_ms, dt_re, dt_xf;
  double ut, ut_a1, ut_x1, ut_x2, ut_t1, ut_t2, ut_d1, ut_d2, ut_u1, ut_u2;
  double ut_ms, ut_re, ut_xf;
  double st, st_tt, st_dt, st_ut, st_ms, st_re, st_a1, st_x1, st_x2, st_t1,
      st_t2;
  double st_d1, st_d2, st_u1, st_u2, st_xf;
  double ctr, ctr_hk2;
  int k;
  //	double c1sav[74], c2sav[74];

  //---- save variables and sensitivities for future restoration
  //	for (int icom=1; icom<= ncom;icom++){
  //		c1sav[icom] = com1[icom];
  //		c2sav[icom] = com2[icom];
  //	}
  saveblData(1);
  saveblData(2);

  //---- weighting factors for linear interpolation to transition point
  wf2 = (xt - blData1.xz) / (blData2.xz - blData1.xz);
  wf2_xt = 1.0 / (blData2.xz - blData1.xz);

  wf2_a1 = wf2_xt * xt_a1;
  wf2_x1 = wf2_xt * xt_x1 + (wf2 - 1.0) / (blData2.xz - blData1.xz);
  wf2_x2 = wf2_xt * xt_x2 - wf2 / (blData2.xz - blData1.xz);
  wf2_t1 = wf2_xt * xt_t1;
  wf2_t2 = wf2_xt * xt_t2;
  wf2_d1 = wf2_xt * xt_d1;
  wf2_d2 = wf2_xt * xt_d2;
  wf2_u1 = wf2_xt * xt_u1;
  wf2_u2 = wf2_xt * xt_u2;
  wf2_ms = wf2_xt * xt_ms;
  wf2_re = wf2_xt * xt_re;
  wf2_xf = wf2_xt * xt_xf;

  wf1 = 1.0 - wf2;
  wf1_a1 = -wf2_a1;
  wf1_x1 = -wf2_x1;
  wf1_x2 = -wf2_x2;
  wf1_t1 = -wf2_t1;
  wf1_t2 = -wf2_t2;
  wf1_d1 = -wf2_d1;
  wf1_d2 = -wf2_d2;
  wf1_u1 = -wf2_u1;
  wf1_u2 = -wf2_u2;
  wf1_ms = -wf2_ms;
  wf1_re = -wf2_re;
  wf1_xf = -wf2_xf;

  //**** first,  do laminar part between x1 and xt

  //-----interpolate primary variables to transition point
  tt = blData1.tz * wf1 + blData2.tz * wf2;
  tt_a1 = blData1.tz * wf1_a1 + blData2.tz * wf2_a1;
  tt_x1 = blData1.tz * wf1_x1 + blData2.tz * wf2_x1;
  tt_x2 = blData1.tz * wf1_x2 + blData2.tz * wf2_x2;
  tt_t1 = blData1.tz * wf1_t1 + blData2.tz * wf2_t1 + wf1;
  tt_t2 = blData1.tz * wf1_t2 + blData2.tz * wf2_t2 + wf2;
  tt_d1 = blData1.tz * wf1_d1 + blData2.tz * wf2_d1;
  tt_d2 = blData1.tz * wf1_d2 + blData2.tz * wf2_d2;
  tt_u1 = blData1.tz * wf1_u1 + blData2.tz * wf2_u1;
  tt_u2 = blData1.tz * wf1_u2 + blData2.tz * wf2_u2;
  tt_ms = blData1.tz * wf1_ms + blData2.tz * wf2_ms;
  tt_re = blData1.tz * wf1_re + blData2.tz * wf2_re;
  tt_xf = blData1.tz * wf1_xf + blData2.tz * wf2_xf;

  dt = blData1.dz * wf1 + blData2.dz * wf2;
  dt_a1 = blData1.dz * wf1_a1 + blData2.dz * wf2_a1;
  dt_x1 = blData1.dz * wf1_x1 + blData2.dz * wf2_x1;
  dt_x2 = blData1.dz * wf1_x2 + blData2.dz * wf2_x2;
  dt_t1 = blData1.dz * wf1_t1 + blData2.dz * wf2_t1;
  dt_t2 = blData1.dz * wf1_t2 + blData2.dz * wf2_t2;
  dt_d1 = blData1.dz * wf1_d1 + blData2.dz * wf2_d1 + wf1;
  dt_d2 = blData1.dz * wf1_d2 + blData2.dz * wf2_d2 + wf2;
  dt_u1 = blData1.dz * wf1_u1 + blData2.dz * wf2_u1;
  dt_u2 = blData1.dz * wf1_u2 + blData2.dz * wf2_u2;
  dt_ms = blData1.dz * wf1_ms + blData2.dz * wf2_ms;
  dt_re = blData1.dz * wf1_re + blData2.dz * wf2_re;
  dt_xf = blData1.dz * wf1_xf + blData2.dz * wf2_xf;

  ut = blData1.uz * wf1 + blData2.uz * wf2;
  ut_a1 = blData1.uz * wf1_a1 + blData2.uz * wf2_a1;
  ut_x1 = blData1.uz * wf1_x1 + blData2.uz * wf2_x1;
  ut_x2 = blData1.uz * wf1_x2 + blData2.uz * wf2_x2;
  ut_t1 = blData1.uz * wf1_t1 + blData2.uz * wf2_t1;
  ut_t2 = blData1.uz * wf1_t2 + blData2.uz * wf2_t2;
  ut_d1 = blData1.uz * wf1_d1 + blData2.uz * wf2_d1;
  ut_d2 = blData1.uz * wf1_d2 + blData2.uz * wf2_d2;
  ut_u1 = blData1.uz * wf1_u1 + blData2.uz * wf2_u1 + wf1;
  ut_u2 = blData1.uz * wf1_u2 + blData2.uz * wf2_u2 + wf2;
  ut_ms = blData1.uz * wf1_ms + blData2.uz * wf2_ms;
  ut_re = blData1.uz * wf1_re + blData2.uz * wf2_re;
  ut_xf = blData1.uz * wf1_xf + blData2.uz * wf2_xf;

  //---- set primary "t" variables at xt  (really placed into "2" variables)
  blData2.xz = xt;
  blData2.tz = tt;
  blData2.dz = dt;
  blData2.uz = ut;

  blData2.amplz = amcrit;
  blData2.sz = 0.0;

  //---- calculate laminar secondary "t" variables
  blkin();
  blvar(1);

  //---- calculate x1-xt midpoint cfm value
  blmid(1);

  //=    at this point, all "2" variables are really "t" variables at xt

  //---- set up newton system for dam, dth, dds, due, dxi  at  x1 and xt
  bldif(1);

  //---- the current newton system is in terms of "1" and "t" variables,
  //-    so calculate its equivalent in terms of "1" and "2" variables.
  //-    in other words, convert residual sensitivities wrt "t" variables
  //-    into sensitivities wrt "1" and "2" variables.  the amplification
  //-    equation is unnecessary here, so the k=1 row is left empty.
  for (k = 2; k <= 3; k++) {
    blrez[k] = vsrez[k];
    blm[k] = vsm[k] + vs2[k][2] * tt_ms + vs2[k][3] * dt_ms +
             vs2[k][4] * ut_ms + vs2[k][5] * xt_ms;
    blr[k] = vsr[k] + vs2[k][2] * tt_re + vs2[k][3] * dt_re +
             vs2[k][4] * ut_re + vs2[k][5] * xt_re;
    blx[k] = vsx[k] + vs2[k][2] * tt_xf + vs2[k][3] * dt_xf +
             vs2[k][4] * ut_xf + vs2[k][5] * xt_xf;

    bl1[k][1] = vs1[k][1] + vs2[k][2] * tt_a1 + vs2[k][3] * dt_a1 +
                vs2[k][4] * ut_a1 + vs2[k][5] * xt_a1;
    bl1[k][2] = vs1[k][2] + vs2[k][2] * tt_t1 + vs2[k][3] * dt_t1 +
                vs2[k][4] * ut_t1 + vs2[k][5] * xt_t1;
    bl1[k][3] = vs1[k][3] + vs2[k][2] * tt_d1 + vs2[k][3] * dt_d1 +
                vs2[k][4] * ut_d1 + vs2[k][5] * xt_d1;
    bl1[k][4] = vs1[k][4] + vs2[k][2] * tt_u1 + vs2[k][3] * dt_u1 +
                vs2[k][4] * ut_u1 + vs2[k][5] * xt_u1;
    bl1[k][5] = vs1[k][5] + vs2[k][2] * tt_x1 + vs2[k][3] * dt_x1 +
                vs2[k][4] * ut_x1 + vs2[k][5] * xt_x1;

    bl2[k][1] = 0.0;
    bl2[k][2] = vs2[k][2] * tt_t2 + vs2[k][3] * dt_t2 + vs2[k][4] * ut_t2 +
                vs2[k][5] * xt_t2;
    bl2[k][3] = vs2[k][2] * tt_d2 + vs2[k][3] * dt_d2 + vs2[k][4] * ut_d2 +
                vs2[k][5] * xt_d2;
    bl2[k][4] = vs2[k][2] * tt_u2 + vs2[k][3] * dt_u2 + vs2[k][4] * ut_u2 +
                vs2[k][5] * xt_u2;
    bl2[k][5] = vs2[k][2] * tt_x2 + vs2[k][3] * dt_x2 + vs2[k][4] * ut_x2 +
                vs2[k][5] * xt_x2;
  }

  //**** second, set up turbulent part between xt and x2  ****

  //---- calculate equilibrium shear coefficient cqt at transition point
  blvar(2);

  //---- set initial shear coefficient value st at transition point
  //-    ( note that cq2, cq2_t2, etc. are really "cqt", "cqt_tt", etc.)

  ctr = 1.8 * exp(-3.3 / (blData2.hkz - 1.0));
  ctr_hk2 = ctr * 3.3 / (blData2.hkz - 1.0) / (blData2.hkz - 1.0);

  st = ctr * blData2.cqz;
  st_tt = ctr * blData2.cqz_tz + blData2.cqz * ctr_hk2 * blData2.hkz_tz;
  st_dt = ctr * blData2.cqz_dz + blData2.cqz * ctr_hk2 * blData2.hkz_dz;
  st_ut = ctr * blData2.cqz_uz + blData2.cqz * ctr_hk2 * blData2.hkz_uz;
  st_ms = ctr * blData2.cqz_ms + blData2.cqz * ctr_hk2 * blData2.hkz_ms;
  st_re = ctr * blData2.cqz_re;

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

  blData2.amplz = 0.0;
  blData2.sz = st;

  //---- recalculate turbulent secondary "t" variables using proper cti
  blvar(2);

  //---- set "1" variables to "t" variables and reset "2" variables
  //-    to their saved turbulent values
  //	for (icom=1; icom<= ncom; icom++){
  //		com1[icom] = com2[icom];
  //		com2[icom] = c2sav[icom];
  //	}
  stepbl();
  restoreblData(2);

  //---- calculate xt-x2 midpoint cfm value
  blmid(2);

  //---- set up newton system for dct, dth, dds, due, dxi  at  xt and x2
  bldif(2);

  //---- convert sensitivities wrt "t" variables into sensitivities
  //-    wrt "1" and "2" variables as done before for the laminar part
  for (k = 1; k <= 3; k++) {
    btrez[k] = vsrez[k];
    btm[k] = vsm[k] + vs1[k][1] * st_ms + vs1[k][2] * tt_ms +
             vs1[k][3] * dt_ms + vs1[k][4] * ut_ms + vs1[k][5] * xt_ms;
    btr[k] = vsr[k] + vs1[k][1] * st_re + vs1[k][2] * tt_re +
             vs1[k][3] * dt_re + vs1[k][4] * ut_re + vs1[k][5] * xt_re;
    btx[k] = vsx[k] + vs1[k][1] * st_xf + vs1[k][2] * tt_xf +
             vs1[k][3] * dt_xf + vs1[k][4] * ut_xf + vs1[k][5] * xt_xf;

    bt1[k][1] = vs1[k][1] * st_a1 + vs1[k][2] * tt_a1 + vs1[k][3] * dt_a1 +
                vs1[k][4] * ut_a1 + vs1[k][5] * xt_a1;
    bt1[k][2] = vs1[k][1] * st_t1 + vs1[k][2] * tt_t1 + vs1[k][3] * dt_t1 +
                vs1[k][4] * ut_t1 + vs1[k][5] * xt_t1;
    bt1[k][3] = vs1[k][1] * st_d1 + vs1[k][2] * tt_d1 + vs1[k][3] * dt_d1 +
                vs1[k][4] * ut_d1 + vs1[k][5] * xt_d1;
    bt1[k][4] = vs1[k][1] * st_u1 + vs1[k][2] * tt_u1 + vs1[k][3] * dt_u1 +
                vs1[k][4] * ut_u1 + vs1[k][5] * xt_u1;
    bt1[k][5] = vs1[k][1] * st_x1 + vs1[k][2] * tt_x1 + vs1[k][3] * dt_x1 +
                vs1[k][4] * ut_x1 + vs1[k][5] * xt_x1;

    bt2[k][1] = vs2[k][1];
    bt2[k][2] = vs2[k][2] + vs1[k][1] * st_t2 + vs1[k][2] * tt_t2 +
                vs1[k][3] * dt_t2 + vs1[k][4] * ut_t2 + vs1[k][5] * xt_t2;
    bt2[k][3] = vs2[k][3] + vs1[k][1] * st_d2 + vs1[k][2] * tt_d2 +
                vs1[k][3] * dt_d2 + vs1[k][4] * ut_d2 + vs1[k][5] * xt_d2;
    bt2[k][4] = vs2[k][4] + vs1[k][1] * st_u2 + vs1[k][2] * tt_u2 +
                vs1[k][3] * dt_u2 + vs1[k][4] * ut_u2 + vs1[k][5] * xt_u2;
    bt2[k][5] = vs2[k][5] + vs1[k][1] * st_x2 + vs1[k][2] * tt_x2 +
                vs1[k][3] * dt_x2 + vs1[k][4] * ut_x2 + vs1[k][5] * xt_x2;
  }

  //---- add up laminar and turbulent parts to get final system
  //-    in terms of honest-to-god "1" and "2" variables.
  vsrez[1] = btrez[1];
  vsrez[2] = blrez[2] + btrez[2];
  vsrez[3] = blrez[3] + btrez[3];
  vsm[1] = btm[1];
  vsm[2] = blm[2] + btm[2];
  vsm[3] = blm[3] + btm[3];
  vsr[1] = btr[1];
  vsr[2] = blr[2] + btr[2];
  vsr[3] = blr[3] + btr[3];
  vsx[1] = btx[1];
  vsx[2] = blx[2] + btx[2];
  vsx[3] = blx[3] + btx[3];
  for (int l = 1; l <= 5; l++) {
    vs1[1][l] = bt1[1][l];
    vs2[1][l] = bt2[1][l];
    vs1[2][l] = bl1[2][l] + bt1[2][l];
    vs2[2][l] = bl2[2][l] + bt2[2][l];
    vs1[3][l] = bl1[3][l] + bt1[3][l];
    vs2[3][l] = bl2[3][l] + bt2[3][l];
  }

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
  int i, is, ibl, j, js, jbl;
  double dui, ue_m;
  for (is = 1; is <= 2; is++) {
    for (ibl = 2; ibl <= nbl[is]; ibl++) {
      i = ipan[ibl][is];

      dui = 0.0;
      for (js = 1; js <= 2; js++) {
        for (jbl = 2; jbl <= nbl[js]; jbl++) {
          j = ipan[jbl][js];
          ue_m = -vti[ibl][is] * vti[jbl][js] * dij[i][j];
          dui = dui + ue_m * mass[jbl][js];
        }
      }

      uedg[ibl][is] = uinv[ibl][is] + dui;
    }
  }
  return true;
}

bool XFoil::uicalc() {
  //--------------------------------------------------------------
  //     sets inviscid ue from panel inviscid tangential velocity
  //--------------------------------------------------------------
  int i, ibl, is;

  for (is = 1; is <= 2; is++) {
    uinv[1][is] = 0.0;
    uinv_a[1][is] = 0.0;
    for (ibl = 2; ibl <= nbl[is]; ibl++) {
      i = ipan[ibl][is];
      uinv[ibl][is] = vti[ibl][is] * qinv[i];
      uinv_a[ibl][is] = vti[ibl][is] * qinv_a[i];
    }
  }

  return true;
}

bool XFoil::update() {
  //------------------------------------------------------------------
  //      adds on newton deltas to boundary layer variables.
  //      checks for excessive changes and underrelaxes if necessary.
  //      calculates max and rms changes.
  //      also calculates the change in the global variable "ac".
  //        if lalfa=true , "ac" is cl
  //        if lalfa=false, "ac" is alpha
  //------------------------------------------------------------------

  int i = 0, is = 0, iv, iw, j, js, jv, ibl, jbl, kbl = 0;
  double unew[IVX][3], u_ac[IVX][3];
  memset(unew, 0, IVX * 3 * sizeof(double));
  memset(u_ac, 0, IVX * 3 * sizeof(double));
  double qnew[IQX], q_ac[IQX];
  memset(qnew, 0, IQX * sizeof(double));
  memset(q_ac, 0, IQX * sizeof(double));
  double dalmax = 0.0, dalmin = 0.0, dclmax = 0.0, dclmin = 0.0;
  double dac = 0.0, dhi = 0.0, dlo = 0.0, dctau, dthet, dmass, duedg, ddstr;
  double dn1, dn2, dn3, dn4, rdn1, rdn2, rdn3, rdn4;
  double dswaki, hklim, msq, dsw;
  double dui, dui_ac, ue_m , uinv_ac, sa = 0.0, ca = 0.0,
         beta = 0.0, beta_msq = 0.0, bfac = 0.0, bfac_msq = 0.0;
  double clnew = 0.0, cl_a = 0.0, cl_ms = 0.0, cl_ac = 0.0, cginc = 0.0;
  double cpg1 = 0.0, cpg1_ms = 0.0, cpi_q = 0.0, cpc_cpi = 0.0, cpg1_ac = 0.0;
  std::string vmxbl;

  //---- max allowable alpha changes per iteration
  dalmax = 0.5 * dtor;
  dalmin = -0.5 * dtor;
  //---- max allowable cl change per iteration
  dclmax = 0.5;
  dclmin = -0.5;
  if (mach_type != MachType::CONSTANT) dclmin = std::max(-0.5, -0.9 * cl);
  hstinv =
      gamm1 * (minf / qinf) * (minf / qinf) / (1.0 + 0.5 * gamm1 * minf * minf);

  //--- calculate new ue distribution assuming no under-relaxation
  //--- also set the sensitivity of ue wrt to alpha or re
  for (is = 1; is <= 2; is++) {
    for (ibl = 2; ibl <= nbl[is]; ibl++) {
      i = ipan[ibl][is];
      dui = 0.0;
      dui_ac = 0.0;
      for (js = 1; js <= 2; js++) {
        for (jbl = 2; jbl <= nbl[js]; jbl++) {
          j = ipan[jbl][js];
          jv = isys[jbl][js];
          ue_m = -vti[ibl][is] * vti[jbl][js] * dij[i][j];
          dui = dui + ue_m * (mass[jbl][js] + vdel[3][1][jv]);
          dui_ac = dui_ac + ue_m * (-vdel[3][2][jv]);
        }
      }

      //------- uinv depends on "ac" only if "ac" is alpha
      if (lalfa)
        uinv_ac = 0.0;
      else
        uinv_ac = uinv_a[ibl][is];

      unew[ibl][is] = uinv[ibl][is] + dui;
      u_ac[ibl][is] = uinv_ac + dui_ac;
    }
  }

  //--- set new qtan from new ue with appropriate sign change

  for (is = 1; is <= 2; is++) {
    for (ibl = 2; ibl <= iblte[is]; ibl++) {
      i = ipan[ibl][is];
      qnew[i] = vti[ibl][is] * unew[ibl][is];
      q_ac[i] = vti[ibl][is] * u_ac[ibl][is];
    }
  }

  //--- calculate new cl from this new qtan
  sa = sin(alfa);
  ca = cos(alfa);

  beta = sqrt(1.0 - minf * minf);
  beta_msq = -0.5 / beta;

  bfac = 0.5 * minf * minf / (1.0 + beta);
  bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;

  clnew = 0.0;
  cl_a = 0.0;
  cl_ms = 0.0;
  cl_ac = 0.0;

  i = 1;
  cginc = 1.0 - (qnew[i] / qinf) * (qnew[i] / qinf);
  cpg1 = cginc / (beta + bfac * cginc);
  cpg1_ms = -cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

  cpi_q = -2.0 * qnew[i] / qinf / qinf;
  cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
  cpg1_ac = cpc_cpi * cpi_q * q_ac[i];

  for (i = 1; i <= n; i++) {
    int ip = i + 1;
    if (i == n) ip = 1;

    cginc = 1.0 - (qnew[ip] / qinf) * (qnew[ip] / qinf);
    const double cpg2 = cginc / (beta + bfac * cginc);
    const double cpg2_ms = -cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

    cpi_q = -2.0 * qnew[ip] / qinf / qinf;
    cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
    const double cpg2_ac = cpc_cpi * cpi_q * q_ac[ip];

    const double dx = (points.row(ip).x() - points.row(i).x()) * ca + (points.row(ip).y() - points.row(i).y()) * sa;
    const double dx_a = -(points.row(ip).x() - points.row(i).x()) * sa + (points.row(ip).y() - points.row(i).y()) * ca;

    const double ag = 0.5 * (cpg2 + cpg1);
    const double ag_ms = 0.5 * (cpg2_ms + cpg1_ms);
    const double ag_ac = 0.5 * (cpg2_ac + cpg1_ac);

    clnew = clnew + dx * ag;
    cl_a = cl_a + dx_a * ag;
    cl_ms = cl_ms + dx * ag_ms;
    cl_ac = cl_ac + dx * ag_ac;

    cpg1 = cpg2;
    cpg1_ms = cpg2_ms;
    cpg1_ac = cpg2_ac;
  }

  //--- initialize under-relaxation factor
  rlx = 1.0;

  if (lalfa) {
    //===== alpha is prescribed: ac is cl

    //---- set change in re to account for cl changing, since re = re(cl)
    dac = (clnew - cl) / (1.0 - cl_ac - cl_ms * 2.0 * minf * minf_cl);

    //---- set under-relaxation factor if re change is too large
    if (rlx * dac > dclmax) rlx = dclmax / dac;
    if (rlx * dac < dclmin) rlx = dclmin / dac;
  } else {
    //===== cl is prescribed: ac is alpha

    //---- set change in alpha to drive cl to prescribed value
    dac = (clnew - clspec) / (0.0 - cl_ac - cl_a);

    //---- set under-relaxation factor if alpha change is too large
    if (rlx * dac > dalmax) rlx = dalmax / dac;
    if (rlx * dac < dalmin) rlx = dalmin / dac;
  }
  rmsbl = 0.0;
  rmxbl = 0.0;
  dhi = 1.5;
  dlo = -.5;
  //--- calculate changes in bl variables and under-relaxation if needed

  for (is = 1; is <= 2; is++) {
    for (ibl = 2; ibl <= nbl[is]; ibl++) {
      iv = isys[ibl][is];
      //------- set changes without underrelaxation
      dctau = vdel[1][1][iv] - dac * vdel[1][2][iv];
      dthet = vdel[2][1][iv] - dac * vdel[2][2][iv];
      dmass = vdel[3][1][iv] - dac * vdel[3][2][iv];
      duedg = unew[ibl][is] + dac * u_ac[ibl][is] - uedg[ibl][is];
      ddstr = (dmass - dstr[ibl][is] * duedg) / uedg[ibl][is];
      //------- normalize changes
      if (ibl < itran[is])
        dn1 = dctau / 10.0;
      else
        dn1 = dctau / ctau[ibl][is];
      dn2 = dthet / thet[ibl][is];
      dn3 = ddstr / dstr[ibl][is];
      dn4 = fabs(duedg) / 0.25;
      //------- accumulate for rms change
      rmsbl = rmsbl + dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;
      //------- see if ctau needs underrelaxation
      rdn1 = rlx * dn1;
      if (fabs(dn1) > fabs(rmxbl)) {
        rmxbl = dn1;
        if (ibl < itran[is]) vmxbl = "n";
        if (ibl >= itran[is]) vmxbl = "c";
        imxbl = ibl;
        ismxbl = is;
      }
      if (rdn1 > dhi) rlx = dhi / dn1;
      if (rdn1 < dlo) rlx = dlo / dn1;
      //------- see if theta needs underrelaxation
      rdn2 = rlx * dn2;
      if (fabs(dn2) > fabs(rmxbl)) {
        rmxbl = dn2;
        vmxbl = "t";
        imxbl = ibl;
        ismxbl = is;
      }
      if (rdn2 > dhi) rlx = dhi / dn2;
      if (rdn2 < dlo) rlx = dlo / dn2;
      //------- see if dstar needs underrelaxation
      rdn3 = rlx * dn3;
      if (fabs(dn3) > fabs(rmxbl)) {
        rmxbl = dn3;
        vmxbl = "d";
        imxbl = ibl;
        ismxbl = is;
      }
      if (rdn3 > dhi) rlx = dhi / dn3;
      if (rdn3 < dlo) rlx = dlo / dn3;

      //------- see if ue needs underrelaxation
      rdn4 = rlx * dn4;
      if (fabs(dn4) > fabs(rmxbl)) {
        rmxbl = duedg;
        vmxbl = "u";
        imxbl = ibl;
        ismxbl = is;
      }
      if (rdn4 > dhi) rlx = dhi / dn4;
      if (rdn4 < dlo) rlx = dlo / dn4;
    }
  }

  //--- set true rms change
  rmsbl = sqrt(rmsbl / (4.0 * double(nbl[1] + nbl[2])));

  if (lalfa) {
    //---- set underrelaxed change in reynolds number from change in lift
    cl = cl + rlx * dac;
  } else {
    //---- set underrelaxed change in alpha
    alfa = alfa + rlx * dac;
    adeg = alfa / dtor;
  }

  //--- update bl variables with underrelaxed changes
  for (is = 1; is <= 2; is++) {
    for (ibl = 2; ibl <= nbl[is]; ibl++) {
      iv = isys[ibl][is];

      dctau = vdel[1][1][iv] - dac * vdel[1][2][iv];
      dthet = vdel[2][1][iv] - dac * vdel[2][2][iv];
      dmass = vdel[3][1][iv] - dac * vdel[3][2][iv];
      duedg = unew[ibl][is] + dac * u_ac[ibl][is] - uedg[ibl][is];
      ddstr = (dmass - dstr[ibl][is] * duedg) / uedg[ibl][is];

      ctau[ibl][is] = ctau[ibl][is] + rlx * dctau;
      thet[ibl][is] = thet[ibl][is] + rlx * dthet;
      dstr[ibl][is] = dstr[ibl][is] + rlx * ddstr;
      uedg[ibl][is] = uedg[ibl][is] + rlx * duedg;

      if (ibl > iblte[is]) {
        iw = ibl - iblte[is];
        dswaki = wgap[iw];
      } else
        dswaki = 0.0;
      //------- eliminate absurd transients
      if (ibl >= itran[is]) ctau[ibl][is] = std::min(ctau[ibl][is], 0.25);

      if (ibl <= iblte[is])
        hklim = 1.02;
      else
        hklim = 1.00005;

      msq = uedg[ibl][is] * uedg[ibl][is] * hstinv /
            (gamm1 * (1.0 - 0.5 * uedg[ibl][is] * uedg[ibl][is] * hstinv));
      dsw = dstr[ibl][is] - dswaki;
      dslim(dsw, thet[ibl][is], msq, hklim);
      dstr[ibl][is] = dsw + dswaki;

      //------- set new mass defect (nonlinear update)
      mass[ibl][is] = dstr[ibl][is] * uedg[ibl][is];
    }
  }

  //--- equate upper wake arrays to lower wake arrays
  for (kbl = 1; kbl <= nbl[2] - iblte[2]; kbl++) {
    ctau[iblte[1] + kbl][1] = ctau[iblte[2] + kbl][2];
    thet[iblte[1] + kbl][1] = thet[iblte[2] + kbl][2];
    dstr[iblte[1] + kbl][1] = dstr[iblte[2] + kbl][2];
    uedg[iblte[1] + kbl][1] = uedg[iblte[2] + kbl][2];
    ctq[iblte[1] + kbl][1] = ctq[iblte[2] + kbl][2];
  }

  return true;
}

bool XFoil::viscal() {
  ////--------------------------------------
  //     converges viscous operating point
  ////--------------------------------------

  //---- calculate wake trajectory from current inviscid solution if necessary
  if (!lwake) xyWake();

  //	---- set velocities on wake from airfoil vorticity for alpha=0, 90
  qwcalc();

  //	---- set velocities on airfoil and wake for initial alpha
  qiset();

  if (!lipan) {
    if (lblini) gamqv();

    //	----- locate stagnation point arc length position and panel index
    stfind();

    //	----- set  bl position -> panel position  pointers
    iblpan();

    //	----- calculate surface arc length array for current stagnation point
    // location
    xicalc();

    //	----- set  bl position -> system line  pointers
    iblsys();
  }

  //	---- set inviscid bl edge velocity uinv from qinv
  uicalc();

  if (!lblini) {
    //	----- set initial ue from inviscid ue
    for (int ibl = 1; ibl <= nbl[1]; ibl++) {
      uedg[ibl][1] = uinv[ibl][1];
    }
    for (int ibl = 1; ibl <= nbl[2]; ibl++) {
      uedg[ibl][2] = uinv[ibl][2];
    }
  }

  if (lvconv) {
    //	----- set correct cl if converged point exists
    qvfue();
    
    if (lvisc) {
      cpcalc(n + nw, qvis, qinf, minf, cpv);
      cpcalc(n + nw, qinv, qinf, minf, cpi);
    } else
      cpcalc(n, qinv, qinf, minf, cpi);

    gamqv();
    clcalc(xcmref, ycmref);
    cdcalc();
  }

  //	---- set up source influence matrix if it doesn't exist
  if (!lwdij || !ladij) qdcalc();

  return true;
}

bool XFoil::ViscalEnd() {

  cpcalc(n + nw, qinv, qinf, minf, cpi);
  cpcalc(n + nw, qvis, qinf, minf, cpv);

  return true;
}

bool XFoil::ViscousIter() {
  //	Performs one iteration
  std::stringstream ss;
  double eps1 = 0.0001;

  setbl();  //	------ fill newton system for bl variables

  blsolve();  //	------ solve newton system with custom solver

  update();  //	------ update bl variables

  if (lalfa) {  //	------- set new freestream mach, re from new cl
    mrcl(cl, minf_cl, reinf_cl);
    comset();
  } else {  //	------- set new inviscid speeds qinv and uinv for new alpha
    qiset();
    uicalc();
  }

  qvfue();   //	------ calculate edge velocities qvis(.) from uedg(..)
  gamqv();   //	------ set gam distribution from qvis
  stmove();  //	------ relocate stagnation point

  //	------ set updated cl,cd
  clcalc(xcmref, ycmref);
  cdcalc();

  //	------ display changes and test for convergence
  if (rlx < 1.0) {
    ss << "     rms:" << std::scientific << std::setprecision(2) << rmsbl
       << "   max:" << std::scientific << std::setprecision(2) << rmxbl
       << " at " << imxbl << " " << ismxbl << "   rlx:" << std::fixed
       << std::setprecision(3) << "\n";
  } else if (fabs(rlx - 1.0) < 0.001) {
    ss << "     rms:" << std::scientific << std::setprecision(2) << rmsbl
       << "   max:" << std::scientific << std::setprecision(2) << rmxbl
       << " at " << imxbl << " " << ismxbl << "\n";
  }

  if (rmsbl < eps1) {
    lvconv = true;
    avisc = alfa;
    mvisc = minf;
    writeString("----------CONVERGED----------\n\n", true);
  }

  return true;
}

bool XFoil::xicalc() {
  //-------------------------------------------------------------
  //     sets bl arc length array on each airfoil side and wake
  //-------------------------------------------------------------
  double telrat, crosp, dwdxte, aa, bb;
  int i, ibl, is, iw;
  is = 1;

  xssi[1][is] = 0.0;

  for (ibl = 2; ibl <= iblte[is]; ibl++) {
    i = ipan[ibl][is];
    xssi[ibl][is] = sst - spline_length[i];
  }

  is = 2;

  xssi[1][is] = 0.0;

  for (ibl = 2; ibl <= iblte[is]; ibl++) {
    i = ipan[ibl][is];
    xssi[ibl][is] = spline_length[i] - sst;
  }

  ibl = iblte[is] + 1;
  xssi[ibl][is] = xssi[ibl - 1][is];

  for (ibl = iblte[is] + 2; ibl <= nbl[is]; ibl++) {
    int i = ipan[ibl][is];
    xssi[ibl][is] =
        xssi[ibl - 1][is] + sqrt((points.row(i).x() - points.row(i - 1).x()) * (points.row(i).x() - points.row(i - 1).x()) +
                                 (points.row(i).y() - points.row(i - 1).y()) * (points.row(i).y() - points.row(i - 1).y()));
  }

  //---- trailing edge flap length to te gap ratio
  telrat = 2.50;

  //---- set up parameters for te flap cubics

  crosp =
      (dpoints_ds.row(1).x() * dpoints_ds.row(n).y() - dpoints_ds.row(1).y() * dpoints_ds.row(n).x()) /
      sqrt((dpoints_ds.row(1).x() * dpoints_ds.row(1).x() + dpoints_ds.row(1).y() * dpoints_ds.row(1).y()) * (dpoints_ds.row(n).x() * dpoints_ds.row(n).x() + dpoints_ds.row(n).y() * dpoints_ds.row(n).y()));
  dwdxte = crosp / sqrt(1.0 - crosp * crosp);

  //---- limit cubic to avoid absurd te gap widths
  dwdxte = std::max(dwdxte, -3.0 / telrat);
  dwdxte = std::min(dwdxte, 3.0 / telrat);

  aa = 3.0 + telrat * dwdxte;
  bb = -2.0 - telrat * dwdxte;

  if (sharp) {
    for (iw = 1; iw <= nw; iw++) wgap[iw] = 0.0;
  }

  else {
    //----- set te flap (wake gap) array
    is = 2;
    for (iw = 1; iw <= nw; iw++) {
      ibl = iblte[is] + iw;
      const double zn = 1.0 - (xssi[ibl][is] - xssi[iblte[is]][is]) / (telrat * ante);
      wgap[iw] = 0.0;
      if (zn >= 0.0) wgap[iw] = ante * (aa + bb * zn) * zn * zn;
    }
  }
  return true;
}

/** -----------------------------------------------------
 * 	   sets forced-transition bl coordinate locations.
 * ----------------------------------------------------- */
bool XFoil::xifset(int is) {
  std::stringstream ss;
  double chx, chy, chsq, str;

  if (xstrip[is] >= 1.0) {
    xiforc = xssi[iblte[is]][is];
    return false;
  }

  chx = point_te.x() - point_le.x();
  chy = point_te.y() - point_le.y();
  chsq = chx * chx + chy * chy;

  //---- calculate chord-based x/c, y/c
  for (int i = 1; i <= n; i++) {
    w1[i] = ((points.row(i).x() - point_le.x()) * chx + (points.row(i).y() - point_le.y()) * chy) / chsq;
    w2[i] = ((points.row(i).y() - point_le.y()) * chx - (points.row(i).x() - point_le.x()) * chy) / chsq;
  }

  w3 = spline::splind(w1, spline_length, n);
  w4 = spline::splind(w2, spline_length, n);

  if (is == 1) {
    //----- set approximate arc length of forced transition point for sinvrt
    str = sle + (spline_length[1] - sle) * xstrip[is];

    //----- calculate actual arc length
    sinvrt(str, xstrip[is], w1, w3, spline_length, n);

    //----- set bl coordinate value
    xiforc = std::min((sst - str), xssi[iblte[is]][is]);
  } else {
    //----- same for bottom side

    str = sle + (spline_length[n] - sle) * xstrip[is];
    sinvrt(str, xstrip[is], w1, w3, spline_length, n);
    xiforc = std::min((str - sst), xssi[iblte[is]][is]);
  }

  if (xiforc < 0.0) {
    ss << " ***  stagnation point is past trip on side " << is << "\n";
    writeString(ss.str());

    xiforc = xssi[iblte[is]][is];
  }

  return true;
}

bool XFoil::xyWake() {
  //-----------------------------------------------------
  //     sets wake coordinate array for current surface
  //     vorticity and/or mass source distributions.
  //-----------------------------------------------------
  double ds1, sx, sy, smod;
  double psi, psi_x, psi_y;
  //
  writeString("   Calculating wake trajectory ...\n", true);
  //
  //--- number of wake points
  nw = n / 8 + 2;
  if (nw > IWX) {
    writeString(
        " XYWake: array size (IWX) too small.\n  Last wake point index "
        "reduced.",
        true);
    nw = IWX;
  }

  ds1 = 0.5 * (spline_length[2] - spline_length[1] + spline_length[n] - spline_length[n - 1]);
  setexp(snew.data() + n, ds1, waklen * chord, nw);

  point_te.x() = 0.5 * (points.row(1).x() + points.row(n).x());
  point_te.y() = 0.5 * (points.row(1).y() + points.row(n).y());

  //-- set first wake point a tiny distance behind te
  int i = n + 1;
  sx = 0.5 * (dpoints_ds.row(n).y() - dpoints_ds.row(1).y());
  sy = 0.5 * (dpoints_ds.row(1).x() - dpoints_ds.row(n).x());
  smod = sqrt(sx * sx + sy * sy);
  nx[i] = sx / smod;
  ny[i] = sy / smod;
  points.row(i).x() = point_te.x() - 0.0001 * ny[i];
  points.row(i).y() = point_te.y() + 0.0001 * nx[i];
  spline_length[i] = spline_length[n];

  //---- calculate streamfunction gradient components at first point
  psilin(i, points.row(i).x(), points.row(i).y(), 1.0, 0.0, psi, psi_x, false, false);
  psilin(i, points.row(i).x(), points.row(i).y(), 0.0, 1.0, psi, psi_y, false, false);

  //---- set unit vector normal to wake at first point
  nx[i + 1] = -psi_x / sqrt(psi_x * psi_x + psi_y * psi_y);
  ny[i + 1] = -psi_y / sqrt(psi_x * psi_x + psi_y * psi_y);

  //---- set angle of wake panel normal
  apanel[i] = atan2(psi_y, psi_x);

  //---- set rest of wake points
  for (i = n + 2; i <= n + nw; i++) {
    const double ds = snew[i] - snew[i - 1];

    //------ set new point ds downstream of last point
    points.row(i).x() = points.row(i - 1).x() - ds * ny[i];
    points.row(i).y() = points.row(i - 1).y() + ds * nx[i];
    spline_length[i] = spline_length[i - 1] + ds;

    if (i != n + nw) {
      //------- calculate normal vector for next point
      psilin(i, points.row(i).x(), points.row(i).y(), 1.0, 0.0, psi, psi_x, false, false);
      psilin(i, points.row(i).x(), points.row(i).y(), 0.0, 1.0, psi, psi_y, false, false);

      nx[i + 1] = -psi_x / sqrt(psi_x * psi_x + psi_y * psi_y);
      ny[i + 1] = -psi_y / sqrt(psi_x * psi_x + psi_y * psi_y);

      //------- set angle of wake panel normal
      apanel[i] = atan2(psi_y, psi_x);
    }
  }

  //---- set wake presence flag and corresponding alpha
  lwake = true;
  awake = alfa;

  //---- old source influence matrix is invalid for the new wake geometry
  lwdij = false;

  return true;
}

bool XFoil::isValidFoilAngles(MatrixX2d points) {
  
  double max_angle = cang(points.middleRows(INDEX_START_WITH, points.rows() - INDEX_START_WITH));
  return max_angle <= angtol;
}

bool XFoil::isValidFoilPointSize(MatrixX2d points) {
  return points.rows() >= 3 + INDEX_START_WITH;
}
