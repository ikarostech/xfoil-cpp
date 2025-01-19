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

#include <cstring>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/StdVector"
using namespace Eigen;

#include "XFoil.h"
// determinant
double cross2(const Eigen::Vector2d& a, const Eigen::Vector2d& b)
{
  return a[0]*b[1] - a[1]*b[0];
}

#define PI 3.141592654

bool XFoil::s_bCancel = false;
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
  xstrip.top = 1.0;
  xstrip.bottom = 1.0;

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
  dtor = PI / 180.0;

  n = 0;  // so that current airfoil is not initialized

  memset(apanel, 0, sizeof(apanel));
  memset(blsav, 0, sizeof(blsav));
  
  bij = MatrixXd::Zero(IQX, IZX);
  memset(cij, 0, sizeof(cij));
  cpi = VectorXd::Zero(IVX);
  cpv = VectorXd::Zero(IVX);
  ctau.top = VectorXd::Zero(IVX);
  ctau.bottom = VectorXd::Zero(IVX);
  ctq.top = VectorXd::Zero(IVX);
  ctq.bottom = VectorXd::Zero(IVX);
  memset(dij, 0, sizeof(dij));
  dstr.top = VectorXd::Zero(IVX);
  dstr.bottom = VectorXd::Zero(IVX);
  
  ipan.top = VectorXi::Zero(IVX);
  ipan.bottom = VectorXi::Zero(IVX);
  isys.top = VectorXi::Zero(IVX);
  isys.bottom = VectorXi::Zero(IVX);
  itran.top = 0;
  itran.bottom = 0;
  mass.top = VectorXd::Zero(IVX);
  mass.bottom = VectorXd::Zero(IVX);
  normal_vectors = Matrix2Xd::Zero(2, IZX);
  gamu = Matrix2Xd::Zero(2, IQX);
  gam = Matrix2Xd::Zero(2, IQX);
  memset(qf0, 0, sizeof(qf0));
  memset(qf1, 0, sizeof(qf1));
  memset(qf2, 0, sizeof(qf2));
  memset(qf3, 0, sizeof(qf3));
  memset(qinv, 0, sizeof(qinv));
  qinvu = Matrix2Xd::Zero(2, IZX);
  memset(qinv_a, 0, sizeof(qinv_a));
  memset(qvis, 0, sizeof(qvis));
  spline_length.resize(IZX);
  snew = VectorXd::Zero(4 * IBX);
  thet.top = VectorXd::Zero(IVX);
  thet.bottom = VectorXd::Zero(IVX);
  uedg.top = VectorXd::Zero(IVX);
  uedg.bottom = VectorXd::Zero(IVX);
  uinv.top = VectorXd::Zero(IVX);
  uinv.bottom = VectorXd::Zero(IVX);
  uinv_a.top = VectorXd::Zero(IVX);
  uinv_a.bottom = VectorXd::Zero(IVX);
  vti.top = VectorXd::Zero(IVX);
  vti.bottom = VectorXd::Zero(IVX);
  points.resize(2, IZX);
  dpoints_ds.resize(2, IZX);
  
  xssi.top = VectorXd::Zero(IVX);
  xssi.bottom = VectorXd::Zero(IVX);
  
  memset(wgap, 0, sizeof(wgap));
  memset(va, 0, sizeof(va));
  memset(vb, 0, sizeof(vb));
  memset(vdel, 0, sizeof(vdel));
  memset(vm, 0, sizeof(vm));
  vs1 = Matrix<double, 4, 5>::Zero();
  vs2 = Matrix<double, 4, 5>::Zero();
  vsrez  = Vector<double, 4>::Zero();
  vsr = Vector<double, 4>::Zero();
  vsm = Vector<double, 4>::Zero();
  vsx = Vector<double, 4>::Zero();
  memset(vz, 0, sizeof(vz));
  
  // mdes
  memset(qgamm, 0, sizeof(qgamm));

  //---- default cp/cv (air)
  gamma = 1.4;
  gamm1 = gamma - 1.0;

  //---- set unity freestream speed
  qinf = 1.0;

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

  //---- circle plane array size (largest 2  + 1 that will fit array size)
  double ann = log(double((2 * IQX) - 1)) / log(2.0);
  int nn = int(ann + 0.00001);
  int tmp = 1;
  int nc1 = 0;
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
  cmref = Vector2d {0.25, 0.0};

  waklen = 1.0;

  // added techwinder : no wake yet
  nw = 0;

  // added techwinder : no flap yet
  hmom = 0.0;

  // added techwinder : fortran initializes to 0
  ist = 0;

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
  gm1bl = 0.0;
  bule = 0.0;
  xiforc = 0.0;
  amcrit = 0.0;

  alfa = 0.0;
  amax = 0.0;
  rmxbl = 0.0;
  rmsbl = 0.0;
  rlx = 0.0;
  ante = 0.0;
  clspec = 0.0;
  minf = 0.0;
  reinf = 0.0;
  minf_cl = 0.0;
  reinf_cl = 0.0;

  sle = 0.0;

  chord = 0.0;
  cl_alf = 0.0;
  cl_msq = 0.0;
  tklam = 0.0;
  tkl_msq = 0.0;
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
  minf_cl = getActualMach(1.0, mach_type);
  reinf_cl = getActualReynolds(1.0, reynolds_type);

  //---- set various compressibility parameters from minf
  comset();

  return true;
}

bool XFoil::abcopy(Matrix2Xd copyFrom) {

  if (n != copyFrom.cols() - 1) lblini = false;

  n = copyFrom.cols() - 1;
  for (int i = 1; i <= n; i++) {
    points.col(i) = copyFrom.col(i);
  }

  //---- strip out doubled points
  int r = 1;

  while (r < n) {
    r++;
    //FIXME double型の==比較
    if (points.col(r - 1) == points.col(r)) {
      for (int j = r; j <= n - 1; j++) {
        points.col(j) = points.col(j + 1);
      }
      n = n - 1;
    }
  }

  spline_length.segment(1, spline_length.size() - 1) = spline::scalc(points.middleCols(1, points.cols() - 1), n, spline_length.size() - 1);
  dpoints_ds.row(0) = spline::splind(points.row(0), spline_length, n);
  dpoints_ds.row(1) = spline::splind(points.row(1), spline_length, n);
  normal_vectors = ncalc(points, spline_length, n);
  lefind(sle, points.middleCols(INDEX_START_WITH, n), dpoints_ds.middleCols(INDEX_START_WITH, n), spline_length.segment(INDEX_START_WITH, n), n);
  point_le.x() = spline::seval(sle, points.row(0), dpoints_ds.row(0), spline_length, n);
  point_le.y() = spline::seval(sle, points.row(1), dpoints_ds.row(1), spline_length, n);
  point_te = 0.5 * (points.col(1) + points.col(n));
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

  //---- set angles of airfoil panels
  for (int i = 1; i <= n - 1; i++) {
    
    Vector2d s = points.col(i + 1) - points.col(i);

    //FIXME double型の==比較
    if (s.norm() == 0.0)
      apanel[i] = atan2(-normal_vectors.col(i).y(), -normal_vectors.col(i).x());
    else
      apanel[i] = atan2(s.x(), -s.y());
  }

  //---- TE panel
  if (sharp)
    apanel[n] = PI;
  else {
    Vector2d s = points.col(1) - points.col(n);
    apanel[n] = atan2(-s.x(), s.y()) + PI;
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
XFoil::AxResult XFoil::axset(double hk1, double t1, double rt1, double a1, double hk2,
                  double t2, double rt2, double a2, double acrit) {
  AxResult result;
  //
  //==========================
  //---- 2nd-order
  double axsq = 0.0;
  double axa = 0.0, axa_ax1 = 0.0, axa_ax2 = 0.0;
  double exn = 0.0, exn_a1 = 0.0, exn_a2 = 0.0, dax = 0.0, dax_a1 = 0.0,
         dax_a2 = 0.0, dax_t1 = 0.0, dax_t2 = 0.0;
  double f_arg = 0.0;  // ex arg

  EnvEnResult envEnResult1 = dampl(hk1, t1, rt1);
  EnvEnResult envEnResult2 = dampl(hk2, t2, rt2);

  //---- rms-average version (seems a little better on coarse grids)
  axsq = 0.5 * (envEnResult1.ax * envEnResult1.ax + envEnResult2.ax * envEnResult2.ax);
  if (axsq <= 0.0) {
    axa = 0.0;
    axa_ax1 = 0.0;
    axa_ax2 = 0.0;
  } else {
    axa = sqrt(axsq);
    axa_ax1 = 0.5 * envEnResult1.ax / axa;
    axa_ax2 = 0.5 * envEnResult2.ax / axa;
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

  result.ax = axa + dax;
  result.ax_hk1 = axa_ax1 * envEnResult1.ax_hk;
  result.ax_t1 = axa_ax1 * envEnResult1.ax_th + dax_t1;
  result.ax_rt1 = axa_ax1 * envEnResult1.ax_rt;
  result.ax_a1 = dax_a1;

  result.ax_hk2 = axa_ax2 * envEnResult2.ax_hk;
  result.ax_t2 = axa_ax2 * envEnResult2.ax_th + dax_t2;
  result.ax_rt2 = axa_ax2 * envEnResult2.ax_rt;
  result.ax_a2 = dax_a2;

  return result;
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
  double upw, upw_hl, upw_hd, upw_hk1, upw_hk2;
  double dxi, slog, scc, scc_usa;
  double rezc, z_ax, hr, hr_hka;
  double hl_hk1, hl_hk2, sa, cqa, cfa, hka, usa, rta, dea, da, ald;
  double gcc, hkc, hkc_hka, rezt, rezh;
  double btmp, hwa, ha, ma, xa, ta, xlog, ulog, tlog, hlog, ddlog;
  double z_cfx, z_ha, z_hwa, z_ma, z_xl, z_tl, z_cfm, z_t1, z_t2;
  double z_dix, z_hca, z_hl, z_hs1, z_hs2, z_di1, z_di2;
  double z_cfa, z_hka, z_da, z_sl, z_ul, z_dxi, z_usa, z_cqa, z_sa, z_dea;
  double z_upw, z_u1, z_u2, z_x1, z_x2;
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
    xlog = log(blData2.param.xz / blData1.param.xz);
    ulog = log(blData2.param.uz / blData1.param.uz);
    tlog = log(blData2.param.tz / blData1.param.tz);
    hlog = log(blData2.hsz.scalar / blData1.hsz.scalar);
    ddlog = 1.0;
  }
  
  vsrez = Vector<double, 4>::Zero();
  vsm = Vector<double, 4>::Zero();
  vsr = Vector<double, 4>::Zero();
  vsx = Vector<double, 4>::Zero();
  vs1 = Matrix<double, 4, 5>::Zero();
  vs2 = Matrix<double, 4, 5>::Zero();
  
  //---- set triggering constant for local upwinding
  hupwt = 1.0;

  hdcon = 5.0 * hupwt / blData2.hkz.scalar / blData2.hkz.scalar;
  hd_hk1 = 0.0;
  hd_hk2 = -hdcon * 2.0 / blData2.hkz.scalar;

  //---- use less upwinding in the wake
  if (ityp == 3) {
    hdcon = hupwt / blData2.hkz.scalar / blData2.hkz.scalar;
    hd_hk1 = 0.0;
    hd_hk2 = -hdcon * 2.0 / blData2.hkz.scalar;
  }
  //
  //---- local upwinding is based on local change in  log(hk-1)
  //-    (mainly kicks in at transition)
  f_arg = fabs((blData2.hkz.scalar - 1.0) / (blData1.hkz.scalar - 1.0));
  hl = log(f_arg);
  hl_hk1 = -1.0 / (blData1.hkz.scalar - 1.0);
  hl_hk2 = 1.0 / (blData2.hkz.scalar - 1.0);

  hlsq = std::min(hl * hl, 15.0);
  ehh = exp(-hlsq * hdcon);
  upw = 1.0 - 0.5 * ehh;
  upw_hl = ehh * hl * hdcon;
  upw_hd = 0.5 * ehh * hlsq;

  upw_hk1 = upw_hl * hl_hk1 + upw_hd * hd_hk1;
  upw_hk2 = upw_hl * hl_hk2 + upw_hd * hd_hk2;

  Vector3d upw1 = upw_hk1 * blData1.hkz.pos_vector();
  Vector3d upw2 = upw_hk2 * blData2.hkz.pos_vector();
  double upw_ms = upw_hk1 * blData1.hkz.ms() + upw_hk2 * blData2.hkz.ms();

  if (ityp == 0) {
    //***** le point -->  set zero amplification factor
    vs2(0, 0) = 1.0;
    vsr[0] = 0.0;
    vsrez[0] = -blData2.param.amplz;
  } else {
    if (ityp == 1) {
      //***** laminar part -->  set amplification equation
      //----- set average amplification ax over interval x1..x2

      AxResult ax_result = axset(blData1.hkz.scalar, blData1.param.tz, blData1.rtz.scalar, blData1.param.amplz, blData2.hkz.scalar, blData2.param.tz, blData2.rtz.scalar, blData2.param.amplz, amcrit);

      rezc = blData2.param.amplz - blData1.param.amplz - ax_result.ax * (blData2.param.xz - blData1.param.xz);
      z_ax = -(blData2.param.xz - blData1.param.xz);

      vs1(0, 0) = z_ax * ax_result.ax_a1 - 1.0;
      vs1(0, 1) = z_ax * (ax_result.ax_hk1 * blData1.hkz.t() + ax_result.ax_t1 + ax_result.ax_rt1 * blData1.rtz.t());
      vs1(0, 2) = z_ax * (ax_result.ax_hk1 * blData1.hkz.d());
      vs1(0, 3) = z_ax * (ax_result.ax_hk1 * blData1.hkz.u() + ax_result.ax_rt1 * blData1.rtz.u());
      vs1(0, 4) = ax_result.ax;
      vs2(0, 0) = z_ax * ax_result.ax_a2 + 1.0;
      vs2(0, 1) = z_ax * (ax_result.ax_hk2 * blData2.hkz.t() + ax_result.ax_t2 + ax_result.ax_rt2 * blData2.rtz.t());
      vs2(0, 2) = z_ax * (ax_result.ax_hk2 * blData2.hkz.d());
      vs2(0, 3) = z_ax * (ax_result.ax_hk2 * blData2.hkz.u() + ax_result.ax_rt2 * blData2.rtz.u());
      vs2(0, 4) = -ax_result.ax;
      vsm[0] = z_ax * (ax_result.ax_hk1 * blData1.hkz.ms() + ax_result.ax_rt1 * blData1.rtz.ms() + ax_result.ax_hk2 * blData2.hkz.ms() +
                       ax_result.ax_rt2 * blData2.rtz.ms());
      vsr[0] = z_ax * (ax_result.ax_rt1 * blData1.rtz.re() + ax_result.ax_rt2 * blData2.rtz.re());
      vsx[0] = 0.0;
      vsrez[0] = -rezc;
    } else {
      //***** turbulent part -->  set shear lag equation

      sa = (1.0 - upw) * blData1.param.sz + upw * blData2.param.sz;
      cqa = (1.0 - upw) * blData1.cqz.scalar + upw * blData2.cqz.scalar;
      cfa = (1.0 - upw) * blData1.cfz.scalar + upw * blData2.cfz.scalar;
      hka = (1.0 - upw) * blData1.hkz.scalar + upw * blData2.hkz.scalar;

      usa = 0.5 * (blData1.usz.scalar + blData2.usz.scalar);
      rta = 0.5 * (blData1.rtz.scalar + blData2.rtz.scalar);
      dea = 0.5 * (blData1.dez.scalar + blData2.dez.scalar);
      da = 0.5 * (blData1.param.dz + blData2.param.dz);

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

      slog = log(blData2.param.sz / blData1.param.sz);
      dxi = blData2.param.xz - blData1.param.xz;

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

      z_upw = z_cqa * (blData2.cqz.scalar - blData1.cqz.scalar) + z_sa * (blData2.param.sz - blData1.param.sz) + z_cfa * (blData2.cfz.scalar - blData1.cfz.scalar) +
              z_hka * (blData2.hkz.scalar - blData1.hkz.scalar);
      double z_de = 0.5 * z_dea;
      double z_us = 0.5 * z_usa;
      double z_d = 0.5 * z_da;
      z_u1 = -z_ul / blData1.param.uz;
      z_u2 = z_ul / blData2.param.uz;
      z_x1 = -z_dxi;
      z_x2 = z_dxi;
      z_s1 = (1.0 - upw) * z_sa - z_sl / blData1.param.sz;
      z_s2 = upw * z_sa + z_sl / blData2.param.sz;
      z_cq1 = (1.0 - upw) * z_cqa;
      z_cq2 = upw * z_cqa;
      z_cf1 = (1.0 - upw) * z_cfa;
      z_cf2 = upw * z_cfa;
      z_hk1 = (1.0 - upw) * z_hka;
      z_hk2 = upw * z_hka;

      vs1(0, 0) = z_s1;
      //vs1(0, 1) = z_us * blData1.usz.t();
      vs1(0, 2) = z_d;// + z_us * blData1.usz.d();
      vs1(0, 3) = z_u1;// + z_us * blData1.usz.u();
      vs1.col(0).segment(1,3) += (z_upw * upw1) + (z_de * blData1.dez.pos_vector()) + (z_us * blData1.usz.pos_vector());
      vs1(0, 4) = z_x1;
      vs2(0, 0) = z_s2;
      //vs2(0, 1) = z_de * blData2.dez.t() + z_us * blData2.usz.t();
      vs2(0, 2) = z_d;// + z_de * blData2.dez.d() + z_us * blData2.usz.d();
      vs2(0, 3) = z_u2;// + z_de * blData2.dez.u() + z_us * blData2.usz.u();
      vs2.col(0).segment(1,3) += (z_upw * upw2) + (z_de * blData2.dez.pos_vector()) + (z_us * blData2.usz.pos_vector());
      vs2(0, 4) = z_x2;
      vsm[0] = z_upw * upw_ms + z_de * blData1.dez.ms() + z_us * blData1.usz.ms() +
               z_de * blData2.dez.ms() + z_us * blData2.usz.ms();

      vs1(0, 1) += z_cq1 * blData1.cqz.t() + z_cf1 * blData1.cfz.t() + z_hk1 * blData1.hkz.t();
      vs1(0, 2) += z_cq1 * blData1.cqz.d() + z_cf1 * blData1.cfz.d() + z_hk1 * blData1.hkz.d();
      vs1(0, 3) += z_cq1 * blData1.cqz.u() + z_cf1 * blData1.cfz.u() + z_hk1 * blData1.hkz.u();
      // ↓ とんでもなく遅いので変数の方を変えるとかしないと無理そう
      //vs1.col(0).segment(1,3) += (z_cq1 * blData1.cqz.pos_vector()) + (z_cf1 * blData1.cfz.pos_vector()) + (z_hk1 * blData1.hkz.pos_vector());

      vs2(0, 1) += z_cq2 * blData2.cqz.t() + z_cf2 * blData2.cfz.t() + z_hk2 * blData2.hkz.t();
      vs2(0, 2) += z_cq2 * blData2.cqz.d() + z_cf2 * blData2.cfz.d() + z_hk2 * blData2.hkz.d();
      vs2(0, 3) += z_cq2 * blData2.cqz.u() + z_cf2 * blData2.cfz.u() + z_hk2 * blData2.hkz.u();

      vsm[0] += z_cq1 * blData1.cqz.ms() + z_cf1 * blData1.cfz.ms() + z_hk1 * blData1.hkz.ms() +
               z_cq2 * blData2.cqz.ms() + z_cf2 * blData2.cfz.ms() + z_hk2 * blData2.hkz.ms();
      vsr[0] =
          z_cq1 * blData1.cqz.re() + z_cf1 * blData1.cfz.re() + z_cq2 * blData2.cqz.re() + z_cf2 * blData2.cfz.re();
      vsx[0] = 0.0;
      vsrez[0] = -rezc;
    }
  }  // endif

  //**** set up momentum equation
  ha = 0.5 * (blData1.param.hz + blData2.param.hz);
  ma = 0.5 * (blData1.param.mz + blData2.param.mz);
  xa = 0.5 * (blData1.param.xz + blData2.param.xz);
  ta = 0.5 * (blData1.param.tz + blData2.param.tz);
  hwa = 0.5 * (blData1.param.dwz / blData1.param.tz + blData2.param.dwz / blData2.param.tz);

  //---- set cf term, using central value cfm for better accuracy in drag
  cfx = 0.50 * cfm * xa / ta + 0.25 * (blData1.cfz.scalar * blData1.param.xz / blData1.param.tz + blData2.cfz.scalar * blData2.param.xz / blData2.param.tz);
  cfx_xa = 0.50 * cfm / ta;
  cfx_ta = -.50 * cfm * xa / ta / ta;

  cfx_x1 = 0.25 * blData1.cfz.scalar / blData1.param.tz + cfx_xa * 0.5;
  cfx_x2 = 0.25 * blData2.cfz.scalar / blData2.param.tz + cfx_xa * 0.5;
  cfx_t1 = -.25 * blData1.cfz.scalar * blData1.param.xz / blData1.param.tz / blData1.param.tz + cfx_ta * 0.5;
  cfx_t2 = -.25 * blData2.cfz.scalar * blData2.param.xz / blData2.param.tz / blData2.param.tz + cfx_ta * 0.5;
  cfx_cf1 = 0.25 * blData1.param.xz / blData1.param.tz;
  cfx_cf2 = 0.25 * blData2.param.xz / blData2.param.tz;
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
      -z_tl / blData1.param.tz + z_cfx * cfx_t1 + z_hwa * 0.5 * (-blData1.param.dwz / blData1.param.tz / blData1.param.tz);
  z_t2 =
      z_tl / blData2.param.tz + z_cfx * cfx_t2 + z_hwa * 0.5 * (-blData2.param.dwz / blData2.param.tz / blData2.param.tz);
  z_x1 = -z_xl / blData1.param.xz + z_cfx * cfx_x1;
  z_x2 = z_xl / blData2.param.xz + z_cfx * cfx_x2;
  z_u1 = -z_ul / blData1.param.uz;
  z_u2 = z_ul / blData2.param.uz;

  vs1(1, 1) = 0.5 * z_ha * blData1.param.hz_tz + z_cfm * cfm_t1 + z_cf1 * blData1.cfz.t() + z_t1;
  vs1(1, 2) = 0.5 * z_ha * blData1.param.hz_dz + z_cfm * cfm_d1 + z_cf1 * blData1.cfz.d();
  vs1(1, 3) = 0.5 * z_ma * blData1.param.mz_uz + z_cfm * cfm_u1 + z_cf1 * blData1.cfz.u() + z_u1;
  vs1(1, 4) = z_x1;
  vs2(1, 1) = 0.5 * z_ha * blData2.param.hz_tz + z_cfm * cfm_t2 + z_cf2 * blData2.cfz.t() + z_t2;
  vs2(1, 2) = 0.5 * z_ha * blData2.param.hz_dz + z_cfm * cfm_d2 + z_cf2 * blData2.cfz.d();
  vs2(1, 3) = 0.5 * z_ma * blData2.param.mz_uz + z_cfm * cfm_u2 + z_cf2 * blData2.cfz.u() + z_u2;
  vs2(1, 4) = z_x2;

  vsm[1] = 0.5 * z_ma * blData1.param.mz_ms + z_cfm * cfm_ms + z_cf1 * blData1.cfz.ms() +
           0.5 * z_ma * blData2.param.mz_ms + z_cf2 * blData2.cfz.ms();
  vsr[1] = z_cfm * cfm_re + z_cf1 * blData1.cfz.re() + z_cf2 * blData2.cfz.re();
  vsx[1] = 0.0;
  vsrez[1] = -rezt;

  //**** set up shape parameter equation

  xot1 = blData1.param.xz / blData1.param.tz;
  xot2 = blData2.param.xz / blData2.param.tz;

  ha = 0.5 * (blData1.param.hz + blData2.param.hz);
  hsa = 0.5 * (blData1.hsz.scalar + blData2.hsz.scalar);
  hca = 0.5 * (blData1.hcz.scalar + blData2.hcz.scalar);
  hwa = 0.5 * (blData1.param.dwz / blData1.param.tz + blData2.param.dwz / blData2.param.tz);

  dix = (1.0 - upw) * blData1.diz.scalar * xot1 + upw * blData2.diz.scalar * xot2;
  cfx = (1.0 - upw) * blData1.cfz.scalar * xot1 + upw * blData2.cfz.scalar * xot2;
  dix_upw = blData2.diz.scalar * xot2 - blData1.diz.scalar * xot1;
  cfx_upw = blData2.cfz.scalar * xot2 - blData1.cfz.scalar * xot1;

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

  z_hs1 = -hca * ulog / hsa / hsa - z_hl / blData1.hsz.scalar;
  z_hs2 = -hca * ulog / hsa / hsa + z_hl / blData2.hsz.scalar;

  z_cf1 = (1.0 - upw) * z_cfx * xot1;
  z_cf2 = upw * z_cfx * xot2;
  z_di1 = (1.0 - upw) * z_dix * xot1;
  z_di2 = upw * z_dix * xot2;

  z_t1 = (1.0 - upw) * (z_cfx * blData1.cfz.scalar + z_dix * blData1.diz.scalar) * (-xot1 / blData1.param.tz);
  z_t2 = upw * (z_cfx * blData2.cfz.scalar + z_dix * blData2.diz.scalar) * (-xot2 / blData2.param.tz);
  z_x1 = (1.0 - upw) * (z_cfx * blData1.cfz.scalar + z_dix * blData1.diz.scalar) / blData1.param.tz - z_xl / blData1.param.xz;
  z_x2 = upw * (z_cfx * blData2.cfz.scalar + z_dix * blData2.diz.scalar) / blData2.param.tz + z_xl / blData2.param.xz;
  z_u1 = -z_ul / blData1.param.uz;
  z_u2 = z_ul / blData2.param.uz;

  z_t1 = z_t1 + z_hwa * 0.5 * (-blData1.param.dwz / blData1.param.tz / blData1.param.tz);
  z_t2 = z_t2 + z_hwa * 0.5 * (-blData2.param.dwz / blData2.param.tz / blData2.param.tz);

  vs1(2, 0) = z_di1 * blData1.diz.s();
  vs1(2, 1) = z_hs1 * blData1.hsz.t() + z_cf1 * blData1.cfz.t() + z_di1 * blData1.diz.t() + z_t1;
  vs1(2, 2) = z_hs1 * blData1.hsz.d() + z_cf1 * blData1.cfz.d() + z_di1 * blData1.diz.d();
  vs1(2, 3) = z_hs1 * blData1.hsz.u() + z_cf1 * blData1.cfz.u() + z_di1 * blData1.diz.u() + z_u1;
  vs1(2, 4) = z_x1;
  vs2(2, 0) = z_di2 * blData2.diz.s();
  vs2(2, 1) = z_hs2 * blData2.hsz.t() + z_cf2 * blData2.cfz.t() + z_di2 * blData2.diz.t() + z_t2;
  vs2(2, 2) = z_hs2 * blData2.hsz.d() + z_cf2 * blData2.cfz.d() + z_di2 * blData2.diz.d();
  vs2(2, 3) = z_hs2 * blData2.hsz.u() + z_cf2 * blData2.cfz.u() + z_di2 * blData2.diz.u() + z_u2;
  vs2(2, 4) = z_x2;
  vsm[2] = z_hs1 * blData1.hsz.ms() + z_cf1 * blData1.cfz.ms() + z_di1 * blData1.diz.ms() + z_hs2 * blData2.hsz.ms() +
           z_cf2 * blData2.cfz.ms() + z_di2 * blData2.diz.ms();
  vsr[2] = z_hs1 * blData1.hsz.re() + z_cf1 * blData1.cfz.re() + z_di1 * blData1.diz.re() + z_hs2 * blData2.hsz.re() +
           z_cf2 * blData2.cfz.re() + z_di2 * blData2.diz.re();

  vs1(2, 1) += 0.5 * (z_hca * blData1.hcz.t() + z_ha * blData1.param.hz_tz) + z_upw * upw1.x();
  vs1(2, 2) += 0.5 * (z_hca * blData1.hcz.d() + z_ha * blData1.param.hz_dz) + z_upw * upw1.y();
  vs1(2, 3) += 0.5 * (z_hca * blData1.hcz.u()) + z_upw * upw1.z();
  vs2(2, 1) += 0.5 * (z_hca * blData2.hcz.t() + z_ha * blData2.param.hz_tz) + z_upw * upw2.x();
  vs2(2, 2) += 0.5 * (z_hca * blData2.hcz.d() + z_ha * blData2.param.hz_dz) + z_upw * upw2.y();
  vs2(2, 3) += 0.5 * (z_hca * blData2.hcz.u()) + z_upw * upw2.z();

  vsm[2] = 0.5 * (z_hca * blData1.hcz.ms()) + z_upw * upw_ms + 0.5 * (z_hca * blData2.hcz.ms());

  vsx[2] = 0.0;
  vsrez[2] = -rezh;

  return true;
}

bool XFoil::blkin() {
  //----------------------------------------------------------
  //     calculates turbulence-independent secondary "2"
  //     variables from the primary "2" variables.
  //----------------------------------------------------------
  double tr2, herat, he_u2, he_ms, v2_he;
  //---- set edge mach number ** 2
  blData2.param.mz = blData2.param.uz * blData2.param.uz * hstinv / (gm1bl * (1.0 - 0.5 * blData2.param.uz * blData2.param.uz * hstinv));
  tr2 = 1.0 + 0.5 * gm1bl * blData2.param.mz;
  blData2.param.mz_uz = 2.0 * blData2.param.mz * tr2 / blData2.param.uz;
  blData2.param.mz_ms = blData2.param.uz * blData2.param.uz * tr2 / (gm1bl * (1.0 - 0.5 * blData2.param.uz * blData2.param.uz * hstinv)) * hstinv_ms;

  //---- set edge density (isentropic relation)
  blData2.param.rz = rstbl * pow(tr2, (-1.0 / gm1bl));
  blData2.param.rz_uz = -blData2.param.rz / tr2 * 0.5 * blData2.param.mz_uz;
  blData2.param.rz_ms = -blData2.param.rz / tr2 * 0.5 * blData2.param.mz_ms + rstbl_ms * pow(tr2, (-1.0 / gm1bl));

  //---- set shape parameter
  blData2.param.hz = blData2.param.dz / blData2.param.tz;
  blData2.param.hz_dz = 1.0 / blData2.param.tz;
  blData2.param.hz_tz = -blData2.param.hz / blData2.param.tz;

  //---- set edge static/stagnation enthalpy
  herat = 1.0 - 0.5 * blData2.param.uz * blData2.param.uz * hstinv;
  he_u2 = -blData2.param.uz * hstinv;
  he_ms = -0.5 * blData2.param.uz * blData2.param.uz * hstinv_ms;
  //---- set molecular viscosity
  v2_he = (1.5 / herat - 1.0 / (herat + hvrat));

  //---- set kinematic shape parameter
  boundary_layer::KineticShapeParameterResult hkin_result = boundary_layer::hkin(blData2.param.hz, blData2.param.mz);
  blData2.hkz.scalar = hkin_result.hk;

  blData2.hkz.u() = hkin_result.hk_msq * blData2.param.mz_uz;
  blData2.hkz.t() = hkin_result.hk_h * blData2.param.hz_tz;
  blData2.hkz.d() = hkin_result.hk_h * blData2.param.hz_dz;
  blData2.hkz.ms() = hkin_result.hk_msq * blData2.param.mz_ms;

  //---- set momentum thickness reynolds number
  blData2.rtz.scalar = blData2.param.rz * blData2.param.uz * blData2.param.tz / (sqrt(herat * herat * herat) * (1.0 + hvrat) / (herat + hvrat) / reybl);
  blData2.rtz.u() = blData2.rtz.scalar * (1.0 / blData2.param.uz + blData2.param.rz_uz / blData2.param.rz - v2_he * he_u2);
  blData2.rtz.t() = blData2.rtz.scalar / blData2.param.tz;
  blData2.rtz.ms() = blData2.rtz.scalar * (blData2.param.rz_ms / blData2.param.rz + (1 / reybl * reybl_ms - v2_he * he_ms));
  blData2.rtz.re() = blData2.rtz.scalar * (reybl_re / reybl);

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
    blData1.rtz = blData2.rtz;
    blData1.param.mz = blData2.param.mz;
    blData1.param.mz_uz = blData2.param.mz_uz;
    blData1.param.mz_ms = blData2.param.mz_ms;
  }

  //---- define stuff for midpoint cf
  hka = 0.5 * (blData1.hkz.scalar + blData2.hkz.scalar);
  rta = 0.5 * (blData1.rtz.scalar + blData2.rtz.scalar);
  ma = 0.5 * (blData1.param.mz + blData2.param.mz);

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
  cfm_u1 = 0.5 * (cfm_hka * blData1.hkz.u() + cfm_ma * blData1.param.mz_uz + cfm_rta * blData1.rtz.u());
  cfm_t1 = 0.5 * (cfm_hka * blData1.hkz.t() + cfm_rta * blData1.rtz.t());
  cfm_d1 = 0.5 * (cfm_hka * blData1.hkz.d());

  cfm_u2 = 0.5 * (cfm_hka * blData2.hkz.u() + cfm_ma * blData2.param.mz_uz + cfm_rta * blData2.rtz.u());
  cfm_t2 = 0.5 * (cfm_hka * blData2.hkz.t() + cfm_rta * blData2.rtz.t());
  cfm_d2 = 0.5 * (cfm_hka * blData2.hkz.d());

  cfm_ms = 0.5 * (cfm_hka * blData1.hkz.ms() + cfm_ma * blData1.param.mz_ms + cfm_rta * blData1.rtz.ms() +
                  cfm_hka * blData2.hkz.ms() + cfm_ma * blData2.param.mz_ms + cfm_rta * blData2.rtz.ms());
  cfm_re = 0.5 * (cfm_rta * blData1.rtz.re() + cfm_rta * blData2.rtz.re());

  return true;
}

/** ----------------------------------------------------------
 *     set bl primary "2" variables from parameter list
 *  ---------------------------------------------------------- */
bool XFoil::blprv(double xsi, double ami, double cti, double thi, double dsi,
                  double dswaki, double uei) {
  blData2.param.xz = xsi;
  blData2.param.amplz = ami;
  blData2.param.sz = cti;
  blData2.param.tz = thi;
  blData2.param.dz = dsi - dswaki;
  blData2.param.dwz = dswaki;

  blData2.param.uz = uei * (1.0 - tkbl) / (1.0 - tkbl * (uei / qinfbl) * (uei / qinfbl));
  blData2.param.uz_uei = (1.0 + tkbl * (2.0 * blData2.param.uz * uei / qinfbl / qinfbl - 1.0)) /
           (1.0 - tkbl * (uei / qinfbl) * (uei / qinfbl));
  blData2.param.uz_ms = (blData2.param.uz * (uei / qinfbl) * (uei / qinfbl) - uei) * tkbl_ms /
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

  int ivte1 = isys.top[iblte.top];
  //
  for (int iv = 1; iv <= nsys; iv++) {
    //
    int ivp = iv + 1;
    //
    //====== invert va[iv] block
    //
    //------ normalize first row
    double pivot = 1.0 / va[0][0][iv];
    va[0][1][iv] *= pivot;
    for (int l = iv; l <= nsys; l++) vm[0][l][iv] *= pivot;
    vdel[0][0][iv] *= pivot;
    vdel[0][1][iv] *= pivot;
    //
    //------ eliminate lower first column in va block
    for (int k = 1; k < 3; k++) {
      vtmp = va[k][0][iv];
      va[k][1][iv] -= vtmp * va[0][1][iv];
      for (int l = iv; l <= nsys; l++) vm[k][l][iv] -= vtmp * vm[0][l][iv];
      vdel[k][0][iv] -= vtmp * vdel[0][0][iv];
      vdel[k][1][iv] -= vtmp * vdel[0][1][iv];
    }
    //
    //------ normalize second row
    pivot = 1.0 / va[1][1][iv];
    for (int l = iv; l <= nsys; l++) vm[1][l][iv] *= pivot;
    vdel[1][0][iv] *= pivot;
    vdel[1][1][iv] *= pivot;
    //
    //------ eliminate lower second column in va block
    
    vtmp = va[2][1][iv];
    for (int l = iv; l <= nsys; l++) vm[2][l][iv] -= vtmp * vm[1][l][iv];
    vdel[2][0][iv] -= vtmp * vdel[1][0][iv];
    vdel[2][1][iv] -= vtmp * vdel[1][1][iv];

    //------ normalize third row
    pivot = 1.0 / vm[2][iv][iv];
    for (int l = ivp; l <= nsys; l++) vm[2][l][iv] *= pivot;
    vdel[2][0][iv] *= pivot;
    vdel[2][1][iv] *= pivot;
    //
    //
    //------ eliminate upper third column in va block
    double vtmp1 = vm[0][iv][iv];
    double vtmp2 = vm[1][iv][iv];
    for (int l = ivp; l <= nsys; l++) {
      vm[0][l][iv] -= vtmp1 * vm[2][l][iv];
      vm[1][l][iv] -= vtmp2 * vm[2][l][iv];
    }
    vdel[0][0][iv] -= vtmp1 * vdel[2][0][iv];
    vdel[1][0][iv] -= vtmp2 * vdel[2][0][iv];
    vdel[0][1][iv] -= vtmp1 * vdel[2][1][iv];
    vdel[1][1][iv] -= vtmp2 * vdel[2][1][iv];
    //
    //------ eliminate upper second column in va block
    vtmp = va[0][1][iv];
    for (int l = ivp; l <= nsys; l++) vm[0][l][iv] -= vtmp * vm[1][l][iv];

    vdel[0][0][iv] -= vtmp * vdel[1][0][iv];
    vdel[0][1][iv] -= vtmp * vdel[1][1][iv];
    //
    //
    if (iv != nsys) {
      //
      //====== eliminate vb(iv+1) block][ rows  1 -> 3
      for (int k = 0; k < 3; k++) {
        vtmp1 = vb[k][0][ivp];
        vtmp2 = vb[k][1][ivp];
        vtmp3 = vm[k][iv][ivp];
        for (int l = ivp; l <= nsys; l++)
          vm[k][l][ivp] -= (vtmp1 * vm[0][l][iv] + vtmp2 * vm[1][l][iv] +
                            vtmp3 * vm[2][l][iv]);
        vdel[k][0][ivp] -= (vtmp1 * vdel[0][0][iv] + vtmp2 * vdel[1][0][iv] +
                            vtmp3 * vdel[2][0][iv]);
        vdel[k][1][ivp] -= (vtmp1 * vdel[0][1][iv] + vtmp2 * vdel[1][1][iv] +
                            vtmp3 * vdel[2][1][iv]);
      }
      //
      if (iv == ivte1) {
        //------- eliminate vz block
        int ivz = isys.bottom[iblte.bottom + 1];
        //
        for (int k = 0; k < 3; k++) {
          vtmp1 = vz[k][0];
          vtmp2 = vz[k][1];
          for (int l = ivp; l <= nsys; l++) {
            vm[k][l][ivz] -= (vtmp1 * vm[0][l][iv] + vtmp2 * vm[1][l][iv]);
          }
          vdel[k][0][ivz] -= (vtmp1 * vdel[0][0][iv] + vtmp2 * vdel[1][0][iv]);
          vdel[k][1][ivz] -= (vtmp1 * vdel[0][1][iv] + vtmp2 * vdel[1][1][iv]);
        }
      }
      //
      if (ivp != nsys) {
        //
        //====== eliminate lower vm column
        for (int kv = iv + 2; kv <= nsys; kv++) {
          vtmp1 = vm[0][iv][kv];
          vtmp2 = vm[1][iv][kv];
          vtmp3 = vm[2][iv][kv];
          //
          if (fabs(vtmp1) > vaccel) {
            for (int l = ivp; l <= nsys; l++) vm[0][l][kv] -= vtmp1 * vm[2][l][iv];
            vdel[0][0][kv] -= vtmp1 * vdel[2][0][iv];
            vdel[0][1][kv] -= vtmp1 * vdel[2][1][iv];
          }
          //
          if (fabs(vtmp2) > vaccel) {
            for (int l = ivp; l <= nsys; l++) vm[1][l][kv] -= vtmp2 * vm[2][l][iv];
            vdel[1][0][kv] -= vtmp2 * vdel[2][0][iv];
            vdel[1][1][kv] -= vtmp2 * vdel[2][1][iv];
          }
          //
          if (fabs(vtmp3) > vaccel) {
            for (int l = ivp; l <= nsys; l++) vm[2][l][kv] -= vtmp3 * vm[2][l][iv];
            vdel[2][0][kv] -= vtmp3 * vdel[2][0][iv];
            vdel[2][1][kv] -= vtmp3 * vdel[2][1][iv];
          }
          //
        }
      }
    }
  }  // 1000

  //
  for (int iv = nsys; iv >= 2; iv--) {
    //------ eliminate upper vm columns
    vtmp = vdel[2][0][iv];
    for (int kv = iv - 1; kv >= 1; kv--) {
      vdel[0][0][kv] -= vm[0][iv][kv] * vtmp;
      vdel[1][0][kv] -= vm[1][iv][kv] * vtmp;
      vdel[2][0][kv] -= vm[2][iv][kv] * vtmp;
    }
    vtmp = vdel[2][1][iv];
    for (int kv = iv - 1; kv >= 1; kv--) {
      vdel[0][1][kv] -= vm[0][iv][kv] * vtmp;
      vdel[1][1][kv] -= vm[1][iv][kv] * vtmp;
      vdel[2][1][kv] -= vm[2][iv][kv] * vtmp;
    }
    //
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
bool XFoil::blvar(blData& ref, int ityp) {
  double hs2_hk2, hs2_rt2, hs2_m2;
  double us2_hs2, us2_hk2, us2_h2;
  double hkc, hkc_hk2, hkc_rt2, hkb, usb;
  double cq2_hs2, cq2_us2, cq2_hk2;
  double cq2_rt2, cq2_h2, cf2_hk2, cf2_rt2, cf2_m2;
  double cf2l, cf2l_hk2, cf2l_rt2, cf2l_m2;
  double cf2t, cf2t_hk2;
  double di2_hs2;
  double di2_us2, di2_cf2t, hmin, hm_rt2;
  double grt, fl, fl_hk2, fl_rt2, tfl;
  double dfac, df_fl, df_hk2, df_rt2;
  double di2l, de2_hk2, hdmax;

  //	double gbcon, gccon, ctcon, hkc2;//were are they initialized?
  if (ityp == 3) ref.hkz.scalar = std::max(ref.hkz.scalar, 1.00005);
  if (ityp != 3) ref.hkz.scalar = std::max(ref.hkz.scalar, 1.05000);

  //---- density thickness shape parameter     ( h** )
  boundary_layer::DensityShapeParameterResult hct_result = boundary_layer::hct(ref.hkz.scalar, ref.param.mz);
  ref.hcz.scalar = hct_result.hc;
  ref.hcz.pos_vector() = hct_result.hc_hk * ref.hkz.pos_vector();
  ref.hcz.u() += hct_result.hc_msq * ref.param.mz_uz;
  ref.hcz.ms()+= hct_result.hc_msq * ref.param.mz_ms;

  //---- set ke thickness shape parameter from  h - h*  correlations
  if (ityp == 1) {
    boundary_layer::ThicknessShapeParameterResult hsl_result = boundary_layer::hsl(ref.hkz.scalar);
    ref.hsz.scalar = hsl_result.hs;
    hs2_hk2 = hsl_result.hs_hk;
    hs2_rt2 = hsl_result.hs_rt;
    hs2_m2 = hsl_result.hs_msq;
  }
  else {
    boundary_layer::ThicknessShapeParameterResult hst_result = boundary_layer::hst(ref.hkz.scalar, ref.rtz.scalar, ref.param.mz);
    ref.hsz.scalar = hst_result.hs;
    hs2_hk2 = hst_result.hs_hk;
    hs2_rt2 = hst_result.hs_rt;
    hs2_m2 = hst_result.hs_msq;
  }

  ref.hsz.vector = hs2_hk2 * ref.hkz.vector + hs2_rt2 * ref.rtz.vector;
  ref.hsz.u() = ref.hsz.u() + hs2_m2 * ref.param.mz_uz;
  ref.hsz.ms() = ref.hsz.ms() + hs2_m2 * ref.param.mz_uz;

  //---- normalized slip velocity  us
  ref.usz.scalar = 0.5 * ref.hsz.scalar * (1.0 - (ref.hkz.scalar - 1.0) / (gbcon * ref.param.hz));
  us2_hs2 = 0.5 * (1.0 - (ref.hkz.scalar - 1.0) / (gbcon * ref.param.hz));
  us2_hk2 = 0.5 * ref.hsz.scalar * (-1.0 / (gbcon * ref.param.hz));
  us2_h2 = 0.5 * ref.hsz.scalar * (ref.hkz.scalar - 1.0) / (gbcon * ref.param.hz * ref.param.hz);

  ref.usz.vector = us2_hs2 * ref.hsz.vector + us2_hk2 * ref.hkz.vector;

  ref.usz.t() = ref.usz.t() + us2_h2 * ref.param.hz_tz;
  ref.usz.d() = ref.usz.d() + us2_h2 * ref.param.hz_dz;

  if (ityp <= 2 && ref.usz.scalar > 0.95) {
    //       write(*,*) 'blvar: us clamped:', us2
    ref.usz.scalar = 0.98;
    ref.usz.vector = Vector<double , 6>::Zero();
  }

  if (ityp == 3 && ref.usz.scalar > 0.99995) {
    //       write(*,*) 'blvar: wake us clamped:', us2
    ref.usz.scalar = 0.99995;
    ref.usz.vector = Vector<double , 6>::Zero();
  }

  //---- equilibrium wake layer shear coefficient (ctau)eq ** 1/2
  //   ...  new  12 oct 94
  hkc = ref.hkz.scalar - 1.0;
  hkc_hk2 = 1.0;
  hkc_rt2 = 0.0;
  if (ityp == 2) {
    const double gcc = gccon;
    hkc = ref.hkz.scalar - 1.0 - gcc / ref.rtz.scalar;
    hkc_hk2 = 1.0;
    hkc_rt2 = gcc / ref.rtz.scalar / ref.rtz.scalar;
    if (hkc < 0.01) {
      hkc = 0.01;
      hkc_hk2 = 0.0;
      hkc_rt2 = 0.0;
    }
  }

  hkb = ref.hkz.scalar - 1.0;
  usb = 1.0 - ref.usz.scalar;
  ref.cqz.scalar = sqrt(ctcon * ref.hsz.scalar * hkb * hkc * hkc / (usb * ref.param.hz * ref.hkz.scalar * ref.hkz.scalar));
  cq2_hs2 = ctcon * hkb * hkc * hkc / (usb * ref.param.hz * ref.hkz.scalar * ref.hkz.scalar) * 0.5 / ref.cqz.scalar;
  cq2_us2 =
      ctcon * ref.hsz.scalar * hkb * hkc * hkc / (usb * ref.param.hz * ref.hkz.scalar * ref.hkz.scalar) / usb * 0.5 / ref.cqz.scalar;
  cq2_hk2 = ctcon * ref.hsz.scalar * hkc * hkc / (usb * ref.param.hz * ref.hkz.scalar * ref.hkz.scalar) * 0.5 / ref.cqz.scalar -
            ctcon * ref.hsz.scalar * hkb * hkc * hkc / (usb * ref.param.hz * ref.hkz.scalar * ref.hkz.scalar * ref.hkz.scalar) * 2.0 *
                0.5 / ref.cqz.scalar +
            ctcon * ref.hsz.scalar * hkb * hkc / (usb * ref.param.hz * ref.hkz.scalar * ref.hkz.scalar) * 2.0 * 0.5 / ref.cqz.scalar *
                hkc_hk2;
  cq2_rt2 = ctcon * ref.hsz.scalar * hkb * hkc / (usb * ref.param.hz * ref.hkz.scalar * ref.hkz.scalar) * 2.0 * 0.5 / ref.cqz.scalar *
            hkc_rt2;
  cq2_h2 =
      -ctcon * ref.hsz.scalar * hkb * hkc * hkc / (usb * ref.param.hz * ref.hkz.scalar * ref.hkz.scalar) / ref.param.hz * 0.5 / ref.cqz.scalar;

  ref.cqz.vector = cq2_hs2 * ref.hsz.vector + cq2_us2 * ref.usz.vector + cq2_hk2 * ref.hkz.vector + cq2_rt2 * ref.rtz.vector;
  ref.cqz.t() = ref.cqz.t() + cq2_h2 * ref.param.hz_tz;
  ref.cqz.d() = ref.cqz.d() + cq2_h2 * ref.param.hz_dz;


  //---- set skin friction coefficient
  if (ityp == 3) {
    //----- wake
    ref.cfz.scalar = 0.0;
    cf2_hk2 = 0.0;
    cf2_rt2 = 0.0;
    cf2_m2 = 0.0;
  } else {
    if (ityp == 1) {
      //----- laminar
      C_f c_f = cfl(ref.hkz.scalar, ref.rtz.scalar);
      ref.cfz.scalar = c_f.cf;
      cf2_hk2 = c_f.hk;
      cf2_rt2 = c_f.rt;
      cf2_m2 = c_f.msq;
    }
    else {
      //----- turbulent
      C_f c_fl = cfl(ref.hkz.scalar, ref.rtz.scalar);
      cf2l = c_fl.cf;
      cf2l_hk2 = c_fl.hk;
      cf2l_rt2 = c_fl.rt;
      cf2l_m2 = c_fl.msq;
      C_f c_ft = cft(ref.hkz.scalar, ref.rtz.scalar, ref.param.mz);
      ref.cfz.scalar = c_ft.cf;
      cf2_hk2 = c_ft.hk;
      cf2_rt2 = c_ft.rt;
      cf2_m2 = c_ft.msq;
      if (cf2l > ref.cfz.scalar) {
        //------- laminar cf is greater than turbulent cf -- use laminar
        //-       (this will only occur for unreasonably small rtheta)
        ref.cfz.scalar = cf2l;
        cf2_hk2 = cf2l_hk2;
        cf2_rt2 = cf2l_rt2;
        cf2_m2 = cf2l_m2;
      }
    }
  }
  ref.cfz.vector = cf2_hk2 * ref.hkz.vector + cf2_rt2 * ref.rtz.vector;
  ref.cfz.u() += cf2_m2 * ref.param.mz_uz;
  ref.cfz.ms() += cf2_m2 * ref.param.mz_ms;

  //---- dissipation function    2 cd / h*
  if (ityp == 1) {
    //----- laminar
    DissipationResult dissipation_result = dil(ref.hkz.scalar, ref.rtz.scalar);
    ref.diz.scalar = dissipation_result.di;
    ref.diz.vector = dissipation_result.di_hk * ref.hkz.vector + dissipation_result.di_rt * ref.rtz.vector;
  } else {
    if (ityp == 2) {
      //----- turbulent wall contribution
      C_f c_ft = cft(ref.hkz.scalar, ref.rtz.scalar, ref.param.mz);
      cf2t = c_ft.cf;
      cf2t_hk2 = c_ft.hk;
      double cf2t_rt2 = c_ft.rt;
      double cf2t_m2 = c_ft.msq;
      Vector<double, 6> cf2t_vector = cf2t_hk2 * ref.hkz.vector + cf2t_rt2 * ref.rtz.vector + Vector<double, 6> {
        cf2t_m2 * ref.param.mz_uz,
        0,
        0,
        0,
        cf2t_m2 * ref.param.mz_ms,
        0
      };

      ref.diz.scalar = (0.5 * cf2t * ref.usz.scalar) * 2.0 / ref.hsz.scalar;
      di2_hs2 = -(0.5 * cf2t * ref.usz.scalar) * 2.0 / ref.hsz.scalar / ref.hsz.scalar;
      di2_us2 = (0.5 * cf2t) * 2.0 / ref.hsz.scalar;
      di2_cf2t = (0.5 * ref.usz.scalar) * 2.0 / ref.hsz.scalar;

      ref.diz.vector = di2_hs2 * ref.hsz.vector + di2_us2 * ref.usz.vector + di2_cf2t * cf2t_vector;

      //----- set minimum hk for wake layer to still exist
      grt = log(ref.rtz.scalar);
      hmin = 1.0 + 2.1 / grt;
      hm_rt2 = -(2.1 / grt / grt) / ref.rtz.scalar;

      //----- set factor dfac for correcting wall dissipation for very low hk
      fl = (ref.hkz.scalar - 1.0) / (hmin - 1.0);
      fl_hk2 = 1.0 / (hmin - 1.0);
      fl_rt2 = (-fl / (hmin - 1.0)) * hm_rt2;

      tfl = tanh(fl);
      dfac = 0.5 + 0.5 * tfl;
      df_fl = 0.5 * (1.0 - tfl * tfl);

      df_hk2 = df_fl * fl_hk2;
      df_rt2 = df_fl * fl_rt2;

      ref.diz.vector = dfac * ref.diz.vector + ref.diz.scalar * (df_hk2 * ref.hkz.vector + df_rt2 * ref.rtz.vector);
      ref.diz.re() -= df_rt2 * ref.rtz.re();

      ref.diz.scalar = ref.diz.scalar * dfac;
    } else {
      //----- zero wall contribution for wake
      ref.diz.scalar = 0.0;
      ref.diz.vector = Vector<double, 6>::Zero();
    }
  }
  //---- add on turbulent outer layer contribution
  if (ityp != 1) {
    double dd = ref.param.sz * ref.param.sz * (0.995 - ref.usz.scalar) * 2.0 / ref.hsz.scalar;
    double dd_hs2 = -ref.param.sz * ref.param.sz * (0.995 - ref.usz.scalar) * 2.0 / ref.hsz.scalar / ref.hsz.scalar;
    double dd_us2 = -ref.param.sz * ref.param.sz * 2.0 / ref.hsz.scalar;
    const double dd_s2 = ref.param.sz * 2.0 * (0.995 - ref.usz.scalar) * 2.0 / ref.hsz.scalar;

    ref.diz.scalar = ref.diz.scalar + dd;
    ref.diz.s() = dd_s2;
    ref.diz.vector = ref.diz.vector + dd_hs2 * ref.hsz.vector + dd_us2 * ref.usz.vector;

    //----- add laminar stress contribution to outer layer cd
    dd = 0.15 * (0.995 - ref.usz.scalar) * (0.995 - ref.usz.scalar) / ref.rtz.scalar * 2.0 / ref.hsz.scalar;
    dd_us2 = -0.15 * (0.995 - ref.usz.scalar) * 2.0 / ref.rtz.scalar * 2.0 / ref.hsz.scalar;
    dd_hs2 = -dd / ref.hsz.scalar;
    const double dd_rt2 = -dd / ref.rtz.scalar;

    ref.diz.scalar = ref.diz.scalar + dd;
    ref.diz.vector = ref.diz.vector + dd_hs2 * ref.hsz.vector + dd_us2 * ref.usz.vector + dd_rt2 * ref.rtz.vector;    
  }

  if (ityp == 2) {
    DissipationResult dissipation_result = dil(ref.hkz.scalar, ref.rtz.scalar);
    if (dissipation_result.di > ref.diz.scalar) {
      //------- laminar cd is greater than turbulent cd -- use laminar
      //-       (this will only occur for unreasonably small rtheta)
      ref.diz.scalar = dissipation_result.di;
      ref.diz.vector = dissipation_result.di_hk * ref.hkz.vector + dissipation_result.di_rt * ref.rtz.vector;
    }
  }

  if (ityp == 3) {
    //------ laminar wake cd
    DissipationResult dissipation_result = dilw(ref.hkz.scalar, ref.rtz.scalar);
    if (di2l > ref.diz.scalar) {
      //------- laminar wake cd is greater than turbulent cd -- use laminar
      //-       (this will only occur for unreasonably small rtheta)
      ref.diz.scalar = dissipation_result.di;
      ref.diz.vector = dissipation_result.di_hk * ref.hkz.vector + dissipation_result.di_rt * ref.rtz.vector;
    }

    //----- double dissipation for the wake (two wake halves)
    ref.diz.scalar = ref.diz.scalar * 2.0;
    ref.diz.vector *= 2;
  }

  //---- bl thickness (delta) from simplified green's correlation
  ref.dez.scalar = (3.15 + 1.72 / (ref.hkz.scalar - 1.0)) * ref.param.tz + ref.param.dz;
  de2_hk2 = (-1.72 / (ref.hkz.scalar - 1.0) / (ref.hkz.scalar - 1.0)) * ref.param.tz;
  ref.dez.vector = de2_hk2 * ref.hkz.vector;
  ref.dez.t() += (3.15 + 1.72 / (ref.hkz.scalar - 1.0));
  ref.dez.d() += 1.0;
  

  hdmax = 12.0;
  if (ref.dez.scalar > hdmax * ref.param.tz) {
    ref.dez.scalar = hdmax * ref.param.tz;
    ref.dez.vector = Vector<double, 6>::Zero();
    ref.dez.t() = hdmax;
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
    blvar(blData2, 3);
    blmid(3);
  } else {
    if (turb || tran) {
      blvar(blData2, 2);
      blmid(2);
    } else {
      blvar(blData2, 1);
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
    vs2 += vs1;
    vs1 = Matrix<double, 4, 5>::Zero();
  }

  //---- change system over into incompressible uei and mach
  for (int k = 0; k < 4; k++) {
    //------ residual derivatives wrt compressible uec
    double res_u1 = vs1(k, 3);
    double res_u2 = vs2(k, 3);
    double res_ms = vsm[k];

    //------ combine with derivatives of compressible  u1,u2 = uec(uei m)
    vs1(k, 3) *= blData1.param.uz_uei;
    vs2(k, 3) *= blData2.param.uz_uei;
    vsm[k] = res_u1 * blData1.param.uz_ms + res_u2 * blData2.param.uz_ms + res_ms;
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
double XFoil::cang(Matrix2Xd points) {
  
  double max_angle = 0;
  //---- go over each point, calculating corner angle
  for (int i = 1; i < points.cols() - 1; i++) {
    Vector2d delta_former = points.col(i) - points.col(i - 1);
    Vector2d delta_later =  points.col(i) - points.col(i + 1);

    double sin = cross2(delta_later, delta_former) / delta_former.norm() / delta_later.norm();
    double delta_angle = asin(sin) * 180.0 / PI;

    max_angle = max(fabs(delta_angle), max_angle);
  }
  return max_angle;
}

bool XFoil::cdcalc() {
  if (!(lvisc && lblini)) {
    cd = 0.0;
    return true;
  }

  //---- set variables at the end of the wake
  double thwake = thet.bottom[nbl.bottom];
  double urat = uedg.bottom[nbl.bottom] / qinf;
  double uewake = uedg.bottom[nbl.bottom] * (1.0 - tklam) / (1.0 - tklam * urat * urat);
  double shwake = dstr.bottom[nbl.bottom] / thet.bottom[nbl.bottom];

  //---- extrapolate wake to downstream infinity using squire-young relation
  //      (reduces errors of the wake not being long enough)
  cd = 2.0 * thwake * pow((uewake / qinf), (0.5 * (5.0 + shwake)));

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

bool XFoil::clcalc(Vector2d ref) {
  
  //-----------------------------------------------------------
  //	   integrates surface pressures to get cl and cm.
  //	   integrates skin friction to get cdf.
  //	   calculates dcl/dalpha for prescribed-cl routines.
  //-----------------------------------------------------------

  xcp = 0.0;

  const double beta = sqrt(1.0 - minf * minf);
  const double beta_msq = -0.5 / beta;

  const double bfac = 0.5 * minf * minf / (1.0 + beta);
  const double bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;

  cl = 0.0;
  cm = 0.0;

  cl_alf = 0.0;
  cl_msq = 0.0;
  
  double cginc = 1.0 - (gam(0, 1) / qinf) * (gam(0, 1) / qinf);
  double cpg1 = cginc / (beta + bfac * cginc);
  double cpg1_msq = -cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

  double cpi_gam = -2.0 * gam(0, 1) / qinf / qinf;
  double cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
  double cpg1_alf = cpc_cpi * cpi_gam * gam(1, 1);

  for (int i = 1; i <= n; i++) {
    int ip = i + 1;
    if (i == n) ip = 1;

    cginc = 1.0 - (gam(0, ip) / qinf) * (gam(0, ip) / qinf);
    double cpg2 = cginc / (beta + bfac * cginc);
    double cpg2_msq = -cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

    cpi_gam = -2.0 * gam(0, ip) / qinf / qinf;
    cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
    double cpg2_alf = cpc_cpi * cpi_gam * gam(1, ip);

    Matrix2d rotateMatrix = Matrix2d {
      {cos(alfa), sin(alfa)},
      {-sin(alfa), cos(alfa)}
    };
    const Vector2d dpoint = rotateMatrix * (points.col(ip) - points.col(i));
    const double dg = cpg2 - cpg1;

    const Vector2d apoint = rotateMatrix * ((points.col(ip) + points.col(i)) / 2 + ref);
    const double ag = 0.5 * (cpg2 + cpg1);

    const double dx_alf = cross2(points.col(ip) - points.col(i), rotateMatrix.row(0));
    const double ag_alf = 0.5 * (cpg2_alf + cpg1_alf);
    const double ag_msq = 0.5 * (cpg2_msq + cpg1_msq);

    cl = cl + dpoint.x() * ag;
    cm = cm - dpoint.dot(ag * apoint + dg * dpoint / 12.0);

    xcp += dpoint.x() * ag * (points.col(ip).x() + points.col(i).x()) / 2.0;

    cl_alf = cl_alf + dpoint.x() * ag_alf + ag * dx_alf;
    cl_msq = cl_msq + dpoint.x() * ag_msq;

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
  return true;
}

/** ---------------------------------------------
 *      sets compressible cp from speed.
 * ---------------------------------------------- */
VectorXd XFoil::cpcalc(int n, const double q[], double qinf, double minf) {
  VectorXd cp = VectorXd::Zero(n + INDEX_START_WITH);
  bool denneg = false;
  double beta, bfac;

  beta = sqrt(1.0 - MathUtil::pow(minf, 2));
  bfac = 0.5 * MathUtil::pow(minf, 2) / (1.0 + beta);

  for (int i = 1; i <= n; i++) {
    const double cpinc = 1.0 - (q[i] / qinf) * (q[i] / qinf);
    const double den = beta + bfac * cpinc;
    cp[i] = cpinc / den;
    if (den <= 0.0) denneg = true;
  }

  if (denneg) {
    writeString(
        "CpCalc: local speed too larger\n Compressibility corrections invalid\n");
  }

  return cp;
}

void XFoil::writeString(std::string str) {
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
XFoil::EnvEnResult XFoil::dampl(double hk, double th, double rt) {
  EnvEnResult result;
  double dgr = 0.08;

  const double hmi = 1.0 / (hk - 1.0);
  const double hmi_hk = -hmi * hmi;

  //---- log10(critical rth) - h   correlation for falkner-skan profiles
  const double aa = 2.492 * pow(hmi, 0.43);
  const double aa_hk = (aa / hmi) * 0.43 * hmi_hk;
  const double bb = tanh(14.0 * hmi - 9.24);
  const double bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk;
  const double grcrit = aa + 0.7 * (bb + 1.0);
  const double grc_hk = aa_hk + 0.7 * bb_hk;
  const double gr = log10(rt);
  const double gr_rt = 1.0 / (2.3025851 * rt);
  if (gr < grcrit - dgr) {
    //----- no amplification for rtheta < rcrit
    result.ax = 0.0;
    result.ax_hk = 0.0;
    result.ax_th = 0.0;
    result.ax_rt = 0.0;
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

    result.ax = (af * dadr / th) * rfac;
    result.ax_hk = (af_hk * dadr / th + af * dadr_hk / th) * rfac +
            (af * dadr / th) * rfac_hk;
    result.ax_th = -(result.ax) / th;
    result.ax_rt = (af * dadr / th) * rfac_rt;
  }

  return result;
}

/** Laminar dissipation function  ( 2 cd/h* )     (from Falkner-Skan)*/
XFoil::DissipationResult XFoil::dil(double hk, double rt) {
  DissipationResult result;
  if (hk < 4.0) {
    result.di = (0.00205 * pow((4.0 - hk), 5.5) + 0.207) / rt;
    result.di_hk = (-.00205 * 5.5 * pow((4.0 - hk), 4.5)) / rt;
  } else {
    double hkb = hk - 4.0;
    double den = 1.0 + 0.02 * hkb * hkb;
    result.di = (-.0016 * hkb * hkb / den + 0.207) / rt;
    result.di_hk =
        (-.0016 * 2.0 * hkb * (1.0 / den - 0.02 * hkb * hkb / den / den)) / rt;
  }
  result.di_rt = -(result.di) / rt;

  return result;
}

XFoil::DissipationResult XFoil::dilw(double hk, double rt) {
  DissipationResult result;
  boundary_layer::ThicknessShapeParameterResult hsl_result = boundary_layer::hsl(hk);
  //---- laminar wake dissipation function  ( 2 cd/h* )
  double rcd = 1.10 * (1.0 - 1.0 / hk) * (1.0 - 1.0 / hk) / hk;
  double rcd_hk = -1.10 * (1.0 - 1.0 / hk) * 2.0 / hk / hk / hk - rcd / hk;

  result.di = 2.0 * rcd / (hsl_result.hs * rt);
  result.di_hk = 2.0 * rcd_hk / (hsl_result.hs * rt) - ((result.di) / hsl_result.hs) * hsl_result.hs_hk;
  result.di_rt = -(result.di) / rt - ((result.di) / hsl_result.hs) * hsl_result.hs_rt;

  return result;
}

bool XFoil::dslim(double &dstr, double thet, double msq, double hklim) {
  const double h = (dstr) / thet;

  boundary_layer::KineticShapeParameterResult hkin_result = boundary_layer::hkin(h, msq);

  const double dh = std::max(0.0, hklim - hkin_result.hk) / hkin_result.hk_h;
  dstr = (dstr) + dh * thet;

  return true;
}

bool XFoil::gamqv() {
  for (int i = 1; i <= n; i++) {
    gam(0, i) = qvis[i];
    gam(1, i) = qinv_a[i];
  }

  return true;
}

bool XFoil::getxyf(Matrix2Xd points, Matrix2Xd dpoints_ds, VectorXd s,
                   int n, double &tops, double &bots, double xf, double &yf) {
  double topy, boty, yrel;

  tops = s[1] + (points.col(1).x() - xf);
  bots = s[n] - (points.col(n).x() - xf);
  spline::sinvrt(tops, xf, points.row(0), dpoints_ds.row(0), s, n);
  spline::sinvrt(bots, xf, points.row(0), dpoints_ds.row(0), s, n);
  topy = spline::seval(tops, points.row(1), dpoints_ds.row(1), s, n);
  boty = spline::seval(bots, points.row(1), dpoints_ds.row(1), s, n);

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

  double res;

  writeString("   Calculating unit vorticity distributions ...\n");

  MatrixXd dpsi_dgam = MatrixXd::Zero(n + 1, n + 1);

  Matrix2Xd psi = Matrix2Xd::Zero(2, n + 1);

  //---- set up matrix system for  psi = psio  on airfoil surface.
  //-    the unknowns are (dgamma)i and dpsio.
  for (int i = 1; i <= n; i++) {
    //------ calculate psi and dpsi/dgamma array for current node
    PsiResult psi_result = psilin(i, points.col(i), normal_vectors.col(i), true);

    const Vector2d res = qinf  * Vector2d {points.col(i).y() , -points.col(i).x()};

    //------ dres/dgamma
    for (int j = 0; j < n; j++) {
      dpsi_dgam(i - 1, j) = psi_result.dzdg[j + INDEX_START_WITH];
    }

    for (int j = 0; j < n; j++) {
      bij(i - INDEX_START_WITH, j) = -psi_result.dzdm[j + INDEX_START_WITH];
    }

    //------ dres/dpsio
    dpsi_dgam(i - 1, n) = -1.0;

    psi.col(i - 1) = -res;
  }

  //---- set Kutta condition
  //-    res = gam(1) + gam[n]
  res = 0.0;

  for (int j = 0; j < n + 1; j++) {
    dpsi_dgam(n, j) = 0;
  }

  dpsi_dgam(n, 0) = 1;
  dpsi_dgam(n, n - 1) = 1;

  psi.col(n).x() = -res;
  psi.col(n).y() = -res;

  //---- set up Kutta condition (no direct source influence)
  for (int j = 0; j < n; j++) bij(n, j) = 0.0;

  if (sharp) {
    //----- set zero internal velocity in TE corner

    //----- set TE bisector angle
    const double ag1 = atan2(-dpoints_ds.col(1).y(), -dpoints_ds.col(1).x());
    const double ag2 = atanc(dpoints_ds.col(n).y(), dpoints_ds.col(n).x(), ag1);
    const double abis = 0.5 * (ag1 + ag2);

    Vector2d bis_vector {cos(abis), sin(abis)};

    //----- minimum panel length adjacent to TE
    const double dsmin = std::min((points.col(1) - points.col(2)).norm(), (points.col(n) - points.col(n - 1)).norm());

    //---- distance of internal control point ahead of sharp TE
    //-    (fraction of smaller panel length adjacent to TE)
    const double bwt = 0.1;

    //----- control point on bisector just ahead of TE point
    const Vector2d bis = point_te - bwt * dsmin * bis_vector;
    const Vector2d normal_bis {-bis_vector.y(), bis_vector.x()};

    //----- set velocity component along bisector line
    PsiResult psi_result = psilin(0, bis, normal_bis, true);
    
    //----- dres/dgamma
    for (int j = 0; j < n; j++) {
      dpsi_dgam(n - 1, j) = psi_result.dqdg[j + INDEX_START_WITH];
    }

    //----- -dres/dmass
    for (int j = 0; j < n; j++) bij(n - 1, j) = -psi_result.dqdm[j];

    //----- dres/dpsio
    dpsi_dgam(n - 1, n);

    //----- -dres/duinf -dres/dvinf
    psi.col(n - 1) = -bis_vector;
  }

  //---- lu-factor coefficient matrix aij
  psi_gamma_lu = FullPivLU<MatrixXd>(dpsi_dgam);
  lqaij = true;
  VectorXd gamu_temp(n + 1);
  //---- solve system for the two vorticity distributions

  gamu_temp = psi_gamma_lu.solve(psi.row(0).transpose());
  
  for (int iu = 0; iu <= n; iu++) {
    gamu.col(iu).x() = gamu_temp[iu];
  }

  gamu_temp = psi_gamma_lu.solve(psi.row(1).transpose());
  for (int iu = 0; iu <= n; iu++) {
    gamu.col(iu).y() = gamu_temp[iu];
  }

  //---- set inviscid alpha=0,90 surface speeds for this geometry
  for (int i = 1; i <= n + 1; i++) {
    qinvu.col(i) = gamu.col(i - INDEX_START_WITH);
  }

  lgamu = true;

  return true;
}

/** -----------------------------------------------------------
 *     sets  bl location -> panel location  pointer array ipan
 * -----------------------------------------------------------*/
bool XFoil::iblpan() {
  int iblmax, ibl;
  std::stringstream ss;

  //-- top surface first
  ibl = 1;
  for (int i = ist; i >= 1; i--) {
    ibl = ibl + 1;
    ipan.top[ibl] = i;
    vti.top[ibl] = 1.0;
  }

  iblte.top = ibl;
  nbl.top = ibl;

  //-- bottom surface next
  ibl = 1;
  for (int i = ist + 1; i <= n; i++) {
    ibl = ibl + 1;
    ipan.bottom[ibl] = i;
    vti.bottom[ibl] = -1.0;
  }

  //-- wake
  iblte.bottom = ibl;

  for (int iw = 1; iw <= nw; iw++) {
    int i = n + iw;
    ibl = iblte.bottom + iw;
    ipan.bottom[ibl] = i;
    vti.bottom[ibl] = -1.0;
  }

  nbl.bottom = iblte.bottom + nw;

  //-- upper wake pointers (for plotting only)
  for (int iw = 1; iw <= nw; iw++) {
    ipan.top[iblte.top + iw] = ipan.bottom[iblte.bottom + iw];
    vti.top[iblte.top + iw] = 1.0;
  }
  iblmax = std::max(iblte.top, iblte.bottom) + nw;
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
    for (int ibl = 2; ibl <= nbl.get(is); ibl++) {
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

/** Loads the Foil's geometry in XFoil,
 *  calculates the normal vectors,
 *  and sets the results in current foil */
bool XFoil::initXFoilGeometry(int fn, const double *fx, const double *fy) {

  Matrix2Xd buffer_points = Matrix2Xd::Zero(2, fn + 1);
  for (int i = 0; i < fn; i++) {
    buffer_points.col(i + 1).x() = fx[i];
    buffer_points.col(i + 1).y() = fy[i];
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
  xstrip.top = XtrTop;
  xstrip.bottom = XtrBot;

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
bool XFoil::lefind(double &sle, Matrix2Xd points, Matrix2Xd dpoints_ds, VectorXd s, int n) {
  int i;
  double dseps;
  //---- convergence tolerance
  dseps = (s[n - 1] - s[0]) * 0.00001;

  //---- set trailing edge point coordinates
  point_te = 0.5 * (points.col(0) + points.col(n - 1));

  //---- get first guess for sle
  for (i = 2; i < n - 2; i++) {
    const Vector2d dpoint_te = points.col(i) - point_te;
    const Vector2d dpoint = points.col(i + 1) - points.col(i);
    const double dotp = dpoint_te.dot(dpoint);
    if (dotp < 0.0) break;
  }

  sle = s[i];

  //---- check for sharp le case
  if (s[i] == s[i - 1]) return false;

  //---- newton iteration to get exact sle value
  for (int iter = 1; iter <= 50; iter++) {
    point_le.x() = spline::seval(sle, points.row(0), dpoints_ds.row(0), s, n);
    point_le.y() = spline::seval(sle, points.row(1), dpoints_ds.row(1), s, n);
    const Vector2d dpoint_ds = {
      spline::deval(sle, points.row(0), dpoints_ds.row(0), s, n),
      spline::deval(sle, points.row(1), dpoints_ds.row(1), s, n)
    };
    const Vector2d dpoint_dd = {
      spline::d2val(sle, points.row(0), dpoints_ds.row(0), s, n),
      spline::d2val(sle, points.row(1), dpoints_ds.row(1), s, n)
    };

    Vector2d chord = point_le - point_te;

    //------ drive dot product between chord line and le tangent to zero
    const double res = chord.dot(dpoint_ds);
    const double ress = dpoint_ds.dot(dpoint_ds) + chord.dot(dpoint_dd);

    //------ newton delta for sle
    double dsle = -res / ress;

    dsle = std::max(dsle, -0.02 * fabs(chord.x() + chord.y()));
    dsle = std::min(dsle, 0.02 * fabs(chord.x() + chord.y()));
    sle = sle + dsle;
    if (fabs(dsle) < dseps) return true;
  }

  sle = s[i];
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
  
  double deps = 0.000005;
  int is = 0, ibl, ibm, itrold, iw = 0, itbl;  // icom

  double senswt = 0.0, thi, uei, dsi, cti, dswaki, ratlen = 0.0;
  double sens = 0.0, sennew = 0.0, msq = 0.0, thm = 0.0, dsm = 0.0,
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
    xiforc = xifset(is);

    //---- set leading edge pressure gradient parameter  x/u du/dx
    ibl = 2;
    bule = 1.0;

    //---- old transition station
    itrold = itran.get(is);

    tran = false;
    turb = false;
    itran.get(is) = iblte.get(is);
    //---- march downstream
    for (ibl = 2; ibl <= nbl.get(is); ibl++) {
      ibm = ibl - 1;

      simi = ibl == 2;
      wake = ibl > iblte.get(is);

      //------ initialize current station to existing variables
      xsi = xssi.get(is)[ibl];
      uei = uedg.get(is)[ibl];
      thi = thet.get(is)[ibl];
      dsi = dstr.get(is)[ibl];

      //------ fixed bug   md 7 june 99
      if (ibl < itrold) {
        ami = ctau.get(is)[ibl];  // ami must be initialized
        cti = 0.03;
      } else {
        cti = ctau.get(is)[ibl];
        if (cti <= 0.0) cti = 0.03;
      }

      if (wake) {
        iw = ibl - iblte.get(is);
        dswaki = wgap[iw];
      } else
        dswaki = 0.0;

      if (ibl <= iblte.get(is))
        dsi = std::max(dsi - dswaki, 1.02000 * thi) + dswaki;
      if (ibl > iblte.get(is)) dsi = std::max(dsi - dswaki, 1.00005 * thi) + dswaki;

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
          ami = blData2.param.amplz;
          if (tran) itran.get(is) = ibl;
          if (!tran) itran.get(is) = ibl + 2;
        }
        if (ibl == iblte.get(is) + 1) {
          tte = thet.top[iblte.top] + thet.bottom[iblte.bottom];
          dte = dstr.top[iblte.top] + dstr.bottom[iblte.bottom] + ante;
          cte = (ctau.top[iblte.top] * thet.top[iblte.top] +
                 ctau.bottom[iblte.bottom] * thet.bottom[iblte.bottom]) /
                tte;
          tesys(cte, tte, dte);
        } else {
          blsys();
        }

        //-------- set stuff at first iteration...
        if (itbl == 1) {
          //--------- set "baseline" ue and hk for forming  ue(hk)  relation
          ueref = blData2.param.uz;
          hkref = blData2.hkz.scalar;

          //--------- if current point ibl was turbulent and is now laminar,
          // then...
          if (ibl < itran.get(is) && ibl >= itrold) {
            //---------- extrapolate baseline hk
            uem = uedg.get(is)[ibl - 1];
            dsm = dstr.get(is)[ibl - 1];
            thm = thet.get(is)[ibl - 1];
            msq =
                uem * uem * hstinv / (gm1bl * (1.0 - 0.5 * uem * uem * hstinv));
            boundary_layer::KineticShapeParameterResult hkin_result = boundary_layer::hkin(dsm / thm, msq);
            hkref = hkin_result.hk;
          }

          //--------- if current point ibl was laminar, then...
          if (ibl < itrold) {
            //---------- reinitialize or extrapolate ctau if it's now turbulent
            if (tran) ctau.get(is)[ibl] = 0.03;
            if (turb) ctau.get(is)[ibl] = ctau.get(is)[ibl - 1];
            if (tran || turb) {
              cti = ctau.get(is)[ibl];
              blData2.param.sz = cti;
            }
          }
        }

        if (simi || ibl == iblte.get(is) + 1) {
          //--------- for similarity station or first wake point, prescribe ue
          vs2(3, 0) = 0.0;
          vs2(3, 1) = 0.0;
          vs2(3, 2) = 0.0;
          vs2(3, 3) = blData2.param.uz_uei;
          vsrez[3] = ueref - blData2.param.uz;
        } else {
          //******** calculate ue-hk characteristic slope
          
          //--------- set unit dhk
          vs2(3, 0) = 0.0;
          vs2(3, 1) = blData2.hkz.t();
          vs2(3, 2) = blData2.hkz.d();
          vs2(3, 3) = blData2.hkz.u() * blData2.param.uz_uei;
          vsrez[3] = 1.0;

          //--------- calculate due response
          double delta_sen = vs2.block(0, 0, 4, 4).fullPivLu().solve(vsrez)[3];

          //--------- set  senswt * (normalized due/dhk)
          sennew = senswt * delta_sen * hkref / ueref;
          if (itbl <= 5)
            sens = sennew;
          else if (itbl <= 15)
            sens = 0.5 * (sens + sennew);

          //--------- set prescribed ue-hk combination
          vs2(3, 0) = 0.0;
          vs2(3, 1) = blData2.hkz.t() * hkref;
          vs2(3, 2) = blData2.hkz.d() * hkref;
          vs2(3, 3) = (blData2.hkz.u() * hkref + sens / ueref) * blData2.param.uz_uei;
          vsrez[3] = -(hkref * hkref) * (blData2.hkz.scalar / hkref - 1.0) -
                     sens * (blData2.param.uz / ueref - 1.0);
        }

        //-------- solve newton system for current "2" station
        vsrez = vs2.block(0, 0, 4, 4).fullPivLu().solve(vsrez);

        //-------- determine max changes and underrelax if necessary
        dmax = std::max(fabs(vsrez[1] / thi), fabs(vsrez[2] / dsi));
        if (ibl >= itran.get(is))
          dmax = std::max(dmax, fabs(vsrez[0] / (10.0 * cti)));

        rlx = 1.0;
        if (dmax > 0.3) rlx = 0.3 / dmax;

        //-------- update as usual
        if (ibl < itran.get(is)) ami = ami + rlx * vsrez[0];
        if (ibl >= itran.get(is)) cti = cti + rlx * vsrez[0];
        thi = thi + rlx * vsrez[1];
        dsi = dsi + rlx * vsrez[2];
        uei = uei + rlx * vsrez[3];

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

        if (dmax <= deps) goto stop110;
      }

      ss << "     mrchdu: convergence failed at " << ibl << " ,  side " << is
         << ", res=" << std::setw(4) << std::fixed << std::setprecision(3)
         << dmax << "\n";
      writeString(ss.str());
      ss.str("");

      if (dmax <= 0.1) goto stop109;
      //------ the current unconverged solution might still be reasonable...

      if (dmax > 0.1) {
        //------- the current solution is garbage --> extrapolate values instead
        if (ibl > 3) {
          if (ibl <= iblte.get(is)) {
            thi = thet.get(is)[ibm] * sqrt(xssi.get(is)[ibl] / xssi.get(is)[ibm]);
            dsi = dstr.get(is)[ibm] * sqrt(xssi.get(is)[ibl] / xssi.get(is)[ibm]);
            uei = uedg.get(is)[ibm];
          } else {
            if (ibl == iblte.get(is) + 1) {
              cti = cte;
              thi = tte;
              dsi = dte;
              uei = uedg.get(is)[ibm];
            } else {
              thi = thet.get(is)[ibm];
              ratlen = (xssi.get(is)[ibl] - xssi.get(is)[ibm]) / (10.0 * dstr.get(is)[ibm]);
              dsi = (dstr.get(is)[ibm] + thi * ratlen) / (1.0 + ratlen);
              uei = uedg.get(is)[ibm];
            }
          }
          if (ibl == itran.get(is)) cti = 0.05;
          if (ibl > itran.get(is)) cti = ctau.get(is)[ibm];
        }
      }

    stop109:
      blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
      blkin();

      //------- check for transition and set appropriate flags and things
      if ((!simi) && (!turb)) {
        trchek();
        ami = blData2.param.amplz;
        if (tran) itran.get(is) = ibl;
        if (!tran) itran.get(is) = ibl + 2;
      }

      //------- set all other extrapolated values for current station
      if (ibl < itran.get(is)) blvar(blData2, 1);
      if (ibl >= itran.get(is)) blvar(blData2, 2);
      if (wake) blvar(blData2, 3);

      if (ibl < itran.get(is)) blmid(1);
      if (ibl >= itran.get(is)) blmid(2);
      if (wake) blmid(3);

      //------ pick up here after the newton iterations
    stop110:
      sens = sennew;

      //------ store primary variables
      if (ibl < itran.get(is))
        ctau.get(is)[ibl] = ami;
      else
        ctau.get(is)[ibl] = cti;
      thet.get(is)[ibl] = thi;
      dstr.get(is)[ibl] = dsi;
      uedg.get(is)[ibl] = uei;
      mass.get(is)[ibl] = dsi * uei;
      ctq.get(is)[ibl] = blData2.cqz.scalar;

      //------ set "1" variables to "2" variables for next streamwise station
      blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
      blkin();

      stepbl();

      //------ turbulent intervals will follow transition interval or te
      if (tran || ibl == iblte.get(is)) {
        turb = true;
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

  double msq, ratlen, dsw, hklim;
  double xsi, uei, ucon, tsq, thi, ami, cti, dsi;
  double dswaki;
  double htest, hktest;
  double cst;
  double cte, dte, tte, dmax, hmax, htarg = 0.0;
  //double cte = dte = tte = dmax = hmax = htarg = 0.0;

  //---- shape parameters for separation criteria
  const double hlmax = 3.8;
  const double htmax = 2.5;

  for (int is = 1; is <= 2; is++) {  // 2000
    ss << "    Side " << is << " ...\n";
    writeString(ss.str());

    //---- set forced transition arc length position
    xiforc = xifset(is);

    //---- initialize similarity station with thwaites' formula
    //	ibl = 2;
    xsi = xssi.get(is)[2];
    uei = uedg.get(is)[2];

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
    itran.get(is) = iblte.get(is);

    //---- march downstream
    for (int ibl = 2; ibl <= nbl.get(is); ibl++) {  // 1000
      int ibm = ibl - 1;
      int iw = ibl - iblte.get(is);
      simi = (ibl == 2);
      wake = ibl > iblte.get(is);

      //------ prescribed quantities
      xsi = xssi.get(is)[ibl];
      uei = uedg.get(is)[ibl];

      if (wake) {
        iw = ibl - iblte.get(is);
        dswaki = wgap[iw];
      } else
        dswaki = 0.0;

      direct = true;

      //------ newton iteration loop for current station
      for (int itbl = 1; itbl <= 25; itbl++) {  // 100

        //-------- assemble 10x3 linearized system for dctau, dth, dds, due, dxi
        //         at the previous "1" station and the current "2" station
        //         (the "1" station coefficients will be ignored)

        blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
        blkin();

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
          tte = thet.top[iblte.top] + thet.bottom[iblte.bottom];
          dte = dstr.top[iblte.top] + dstr.bottom[iblte.bottom] + ante;
          cte = (ctau.top[iblte.top] * thet.top[iblte.top] +
                 ctau.bottom[iblte.bottom] * thet.bottom[iblte.bottom]) /
                tte;
          tesys(cte, tte, dte);
        } else
          blsys();

        if (direct) {
          //--------- try direct mode (set due = 0 in currently empty 4th line)
          vs2(3, 0) = 0.0;
          vs2(3, 1) = 0.0;
          vs2(3, 2) = 0.0;
          vs2(3, 3) = 1.0;
          vsrez[3] = 0.0;
          //--------- solve newton system for current "2" station
          vsrez = vs2.block(0, 0, 4, 4).fullPivLu().solve(vsrez);
          //--------- determine max changes and underrelax if necessary
          dmax = std::max(fabs(vsrez[1] / thi), fabs(vsrez[2] / dsi));
          if (ibl < itran.get(is)) dmax = std::max(dmax, fabs(vsrez[0] / 10.0));
          if (ibl >= itran.get(is)) dmax = std::max(dmax, fabs(vsrez[0] / cti));

          rlx = 1.0;
          if (dmax > 0.3) rlx = 0.3 / dmax;
          //--------- see if direct mode is not applicable
          if (ibl != iblte.get(is) + 1) {
            //---------- calculate resulting kinematic shape parameter hk
            msq =
                uei * uei * hstinv / (gm1bl * (1.0 - 0.5 * uei * uei * hstinv));
            htest = (dsi + rlx * vsrez[2]) / (thi + rlx * vsrez[1]);
            boundary_layer::KineticShapeParameterResult hkin_result = boundary_layer::hkin(htest, msq);
            hktest = hkin_result.hk;

            //---------- decide whether to do direct or inverse problem based on
            // hk
            if (ibl < itran.get(is)) hmax = hlmax;
            if (ibl >= itran.get(is)) hmax = htmax;
            direct = (hktest < hmax);
          }
          if (direct) {
            //---------- update as usual
            if (ibl >= itran.get(is)) cti = cti + rlx * vsrez[0];
            thi = thi + rlx * vsrez[1];
            dsi = dsi + rlx * vsrez[2];
          } else {
            //---------- set prescribed hk for inverse calculation at the
            // current station
            if (ibl < itran.get(is))
              //----------- laminar case: relatively slow increase in hk
              // downstream
              htarg = blData1.hkz.scalar + 0.03 * (blData2.param.xz - blData1.param.xz) / blData1.param.tz;
            else if (ibl == itran.get(is)) {
              //----------- transition interval: weighted laminar and turbulent
              // case
              htarg = blData1.hkz.scalar + (0.03 * (xt - blData1.param.xz) - 0.15 * (blData2.param.xz - xt)) / blData1.param.tz;
            } else if (wake) {
              //----------- turbulent wake case:
              //--          asymptotic wake behavior with approximate backward
              // euler
              
              cst = 0.03 * (blData2.param.xz - blData1.param.xz) / blData1.param.tz;
              auto euler = [] (double hk2, double hk1, double cst) -> double {
                return hk2 - (hk2 + cst * pow(hk2 - 1 , 3) - hk1) / (1 + 3 * cst * pow(hk2 - 1, 2));
              };
              
              blData2.hkz.scalar = blData1.hkz.scalar;
              
              for (int i=0; i<3; i++) {
                blData2.hkz.scalar = euler(blData2.hkz.scalar, blData1.hkz.scalar, cst);
              }
              htarg = blData2.hkz.scalar;
            } else
              htarg =
                  blData1.hkz.scalar - 0.15 * (blData2.param.xz - blData1.param.xz) /
                            blData1.param.tz;  //----------- turbulent case: relatively
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

            continue;
          }
        } else {
          //-------- inverse mode (force hk to prescribed value htarg)
          vs2(3, 0) = 0.0;
          vs2(3, 1) = blData2.hkz.t();
          vs2(3, 2) = blData2.hkz.d();
          vs2(3, 3) = blData2.hkz.u();
          vsrez[3] = htarg - blData2.hkz.scalar;
          vsrez = vs2.block(0, 0, 4, 4).fullPivLu().solve(vsrez);

          dmax = std::max(fabs(vsrez[1] / thi), fabs(vsrez[2] / dsi));
          if (ibl >= itran.get(is)) dmax = std::max(dmax, fabs(vsrez[0] / cti));
          rlx = 1.0;
          if (dmax > 0.3) rlx = 0.3 / dmax;
          //--------- update variables
          if (ibl >= itran.get(is)) cti = cti + rlx * vsrez[0];
          thi = thi + rlx * vsrez[1];
          dsi = dsi + rlx * vsrez[2];
          uei = uei + rlx * vsrez[3];
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
        if (dmax <= 0.00001) goto stop110;

      }  // end itbl loop

      ss.str("");
      ss << "     mrchue: convergence failed at " << ibl << ",  side " << is
         << ", res =" << std::fixed << std::setprecision(3) << dmax << "\n";
      writeString(ss.str());

      //------ the current unconverged solution might still be reasonable...
      if (dmax > 0.1) {
        //------- the current solution is garbage --> extrapolate values instead
        if (ibl > 3) {
          if (ibl <= iblte.get(is)) {
            thi = thet.get(is)[ibm] * sqrt(xssi.get(is)[ibl] / xssi.get(is)[ibm]);
            dsi = dstr.get(is)[ibm] * sqrt(xssi.get(is)[ibl] / xssi.get(is)[ibm]);
          } else {
            if (ibl == iblte.get(is) + 1) {
              cti = cte;
              thi = tte;
              dsi = dte;
            } else {
              thi = thet.get(is)[ibm];
              ratlen = (xssi.get(is)[ibl] - xssi.get(is)[ibm]) / (10.0 * dstr.get(is)[ibm]);
              dsi = (dstr.get(is)[ibm] + thi * ratlen) / (1.0 + ratlen);
            }
          }
          if (ibl == itran.get(is)) cti = 0.05;
          if (ibl > itran.get(is)) cti = ctau.get(is)[ibm];

          uei = uedg.get(is)[ibl];

          if (ibl < nbl.get(is))
            uei = 0.5 * (uedg.get(is)[ibl - 1] + uedg.get(is)[ibl + 1]);
        }
      }
      // 109
      blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
      blkin();
      //------- check for transition and set appropriate flags and things
      if ((!simi) && (!turb)) {
        trchek();
        ami = blData2.param.amplz;
        if (tran) itran.get(is) = ibl;
        if (!tran) itran.get(is) = ibl + 2;
      }
      //------- set all other extrapolated values for current station
      if (ibl < itran.get(is)) blvar(blData2, 1);
      if (ibl >= itran.get(is)) blvar(blData2, 2);
      if (wake) blvar(blData2, 3);
      if (ibl < itran.get(is)) blmid(1);
      if (ibl >= itran.get(is)) blmid(2);
      if (wake) blmid(3);
      //------ pick up here after the newton iterations
    stop110:
      //------ store primary variables
      if (ibl < itran.get(is)) ctau.get(is)[ibl] = ami;
      if (ibl >= itran.get(is)) ctau.get(is)[ibl] = cti;
      thet.get(is)[ibl] = thi;
      dstr.get(is)[ibl] = dsi;
      uedg.get(is)[ibl] = uei;
      mass.get(is)[ibl] = dsi * uei;
      ctq.get(is)[ibl] = blData2.cqz.scalar;

      //------ set "1" variables to "2" variables for next streamwise station
      blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
      blkin();

      stepbl();

      //------ turbulent intervals will follow transition interval or te
      if (tran || ibl == iblte.get(is)) {
        turb = true;
      }

      tran = false;

      if (ibl == iblte.get(is)) {
        thi = thet.top[iblte.top] + thet.bottom[iblte.bottom];
        dsi = dstr.top[iblte.top] + dstr.bottom[iblte.bottom] + ante;
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
double XFoil::getActualMach(double cls, MachType mach_type) {
  const double cla = std::max(cls, 0.000001);
  switch (mach_type) {
    case MachType::CONSTANT: {
      minf = minf1;
      return 0.0;
    }
    case MachType::FIXED_LIFT: {
      minf = minf1 / sqrt(cla);
      return -0.5 * minf / cla;
    }
    case MachType::FIXED_LIFT_AND_DYNAMIC_PRESSURE: {
      minf = minf1;
      return 0.0;
    }
    default:
      return 0;
  }
}
double XFoil::getActualReynolds(double cls, ReynoldsType reynolds_type) {
  const double cla = std::max(cls, 0.000001);
  switch (reynolds_type) {
    case ReynoldsType::CONSTANT: {
      reinf = reinf1;
      return 0.0;
    }
    case ReynoldsType::FIXED_LIFT: {
      reinf = reinf1 / sqrt(cla);
      return -0.5 * reinf / cla;
    }
    case ReynoldsType::FIXED_LIFT_AND_DYNAMIC_PRESSURE: {
      reinf = reinf1 / cla;
      return -reinf / cla;
    }
    default:
      return 0;
  }
}

Matrix2Xd XFoil::ncalc(Matrix2Xd points, VectorXd spline_length, int n) {

  Matrix2Xd normal_vector = Matrix2Xd::Zero(2, IZX);
  normal_vector.row(0) = spline::splind(points.row(0), spline_length, n);
  normal_vector.row(1) = spline::splind(points.row(0), spline_length, n);
  for (int i = 1; i <= n; i++) {
    Vector2d temp = normal_vector.col(i);
    normal_vector.col(i).x() = temp.normalized().y();
    normal_vector.col(i).y() = temp.normalized().x();
  }

  return normal_vector;
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
 *			airfoil:  1   < i < n
 *			wake:	  n+1 < i < n+nw
 * ----------------------------------------------------------------------- */
XFoil::PsiResult XFoil::psilin(int iNode, Vector2d point, Vector2d normal_vector, bool siglin) {
  PsiResult psi_result;
  
  //---- distance tolerance for determining if two points are the same
  const double seps = (spline_length[n] - spline_length[1]) * 0.00001;

  psi_result.psi = 0.0;
  psi_result.psi_ni = 0.0;

  psi_result.qtan = Vector2d::Zero();

  for (int jo = 0; jo < n; jo++) {
    int jp = (jo + 1) % n;
    double dso = (points.col(jo + INDEX_START_WITH) - points.col(jp + INDEX_START_WITH)).norm();

    //------ skip null panel
    if (fabs(dso) < 1.0e-7) continue;

    Vector2d r1 = point - points.col(jo + INDEX_START_WITH);
    Vector2d r2 = point - points.col(jp + INDEX_START_WITH);
    Vector2d s = (points.col(jp + INDEX_START_WITH) - points.col(jo + INDEX_START_WITH)).normalized();

    blData1.param.xz = s.dot(r1);
    blData2.param.xz = s.dot(r2);
    double yy = cross2(s, r1);

    //FIXME normを使うと正しく計算されない。計算精度の問題？
    double rs1 = r1.dot(r1);
    double rs2 = r2.dot(r2);

    //------ set reflection flag sgn to avoid branch problems with arctan
    double sgn;
    if (iNode >= 1 && iNode <= n) {
      //------- no problem on airfoil surface
      sgn = 1.0;
    } else {
      //------- make sure arctan falls between  -/+  pi/2
      sgn = sign(1.0, yy);
    }

    //------ set log(r^2) and arctan(x/y), correcting for reflection if any
    double logr12;
    if (iNode != jo + INDEX_START_WITH && rs1 > 0.0) {
      logr12 = log(rs1);
      blData1.param.tz = atan2(sgn * blData1.param.xz, sgn * yy) + (0.5 - 0.5 * sgn) * PI;
    } else {
      logr12 = 0.0;
      blData1.param.tz = 0.0;
    }
    double logr22;
    if (iNode != jp + INDEX_START_WITH && rs2 > 0.0) {
      logr22 = log(rs2);
      blData2.param.tz = atan2(sgn * blData2.param.xz, sgn * yy) + (0.5 - 0.5 * sgn) * PI;
    } else {
      logr22 = 0.0;
      blData2.param.tz = 0.0;
    }

    double x1i = s.dot(normal_vector);
    double x2i = s.dot(normal_vector);
    double yyi = cross2(s, normal_vector);
    if (jo + INDEX_START_WITH == n) break;
    if (siglin) {
      PsiResult sig_result = psisig(iNode, jo, point, normal_vector);
      psi_result = PsiResult::sum(psi_result, sig_result);
    }

    //------ calculate vortex panel contribution to psi
    double dxinv = 1.0 / (blData1.param.xz - blData2.param.xz);
    double psis = 0.5 * blData1.param.xz * logr12 - 0.5 * blData2.param.xz * logr22 + blData2.param.xz - blData1.param.xz +
           yy * (blData1.param.tz - blData2.param.tz);
    double psid = ((blData1.param.xz + blData2.param.xz) * psis +
            0.5 * (rs2 * logr22 - rs1 * logr12 + blData1.param.xz * blData1.param.xz - blData2.param.xz * blData2.param.xz)) *
           dxinv;

    double psx1 = 0.5 * logr12;
    double psx2 = -.5 * logr22;
    double psyy = blData1.param.tz - blData2.param.tz;

    double pdx1 = ((blData1.param.xz + blData2.param.xz) * psx1 + psis - blData1.param.xz * logr12 - psid) * dxinv;
    double pdx2 = ((blData1.param.xz + blData2.param.xz) * psx2 + psis + blData2.param.xz * logr22 + psid) * dxinv;
    double pdyy = ((blData1.param.xz + blData2.param.xz) * psyy - yy * (logr12 - logr22)) * dxinv;

    Vector2d gsum_vector = gamu.col(jp) + gamu.col(jo);
    Vector2d gdif_vector = gamu.col(jp) - gamu.col(jo);

    double gsum = gam(0, jp + INDEX_START_WITH) + gam(0, jo + INDEX_START_WITH);
    double gdif = gam(0, jp + INDEX_START_WITH) - gam(0, jo + INDEX_START_WITH);

    psi_result.psi += (1 / (4 * PI)) * (psis * gsum + psid * gdif);

    //------ dpsi/dgam
    psi_result.dzdg[jo + INDEX_START_WITH] += (1 / (4 * PI)) * (psis - psid);
    psi_result.dzdg[jp + INDEX_START_WITH] += (1 / (4 * PI)) * (psis + psid);

    //------ dpsi/dni
    double psni = psx1 * x1i + psx2 * x2i + psyy * yyi;
    double pdni = pdx1 * x1i + pdx2 * x2i + pdyy * yyi;
    psi_result.psi_ni += (1 / (4 * PI)) * (gsum * psni + gdif * pdni);

    psi_result.qtan += (1 / (4 * PI)) * (psni * gsum_vector + pdni * gdif_vector);

    psi_result.dqdg[jo + INDEX_START_WITH] += (1 / (4 * PI)) * (psni - pdni);
    psi_result.dqdg[jp + INDEX_START_WITH] += (1 / (4 * PI)) * (psni + pdni);
  
  }
  if ((points.col(n) - points.col(1)).norm() > seps) {
    PsiResult te_result = psi_te(iNode, point, normal_vector);
    psi_result = PsiResult::sum(psi_result, te_result);
  }

  //**** freestream terms
  Vector2d rotateVector = {cos(alfa), sin(alfa)};
  psi_result.psi += qinf * cross2(rotateVector, point);

  //---- dpsi/dn
  psi_result.psi_ni += qinf * cross2(rotateVector, normal_vector);

  psi_result.qtan.x() += qinf * normal_vector.y();
  psi_result.qtan.y() += -qinf * normal_vector.x();

  return psi_result;
}

XFoil::PsiResult XFoil::psisig(int iNode, int jo, Vector2d point, Vector2d normal_vector) {
  PsiResult psi_result;
  psi_result.psi = 0;
  psi_result.psi_ni = 0;

  int io = iNode;

  int jp = (jo == n - 1) ? 0 : jo + 1;
  int jm = std::max(0, jo - 1);
  int jq = (jo >= n - 2) ? jp : jp + 1;

  double dso = (points.col(jo + INDEX_START_WITH) - points.col(jp + INDEX_START_WITH)).norm();


  double dsio = 1.0 / dso;

  double apan = apanel[jo + INDEX_START_WITH];

  Vector2d r1 = point - points.col(jo + INDEX_START_WITH);
  Vector2d r2 = point - points.col(jp + INDEX_START_WITH);
  Vector2d s = (points.col(jp + INDEX_START_WITH)- points.col(jo + INDEX_START_WITH)).normalized();

  double x1 = s.dot(r1);
  double x2 = s.dot(r2);
  double yy = cross2(s, r1);

  //FIXME normを使うと正しく計算されない。計算精度の問題？
  double rs1 = r1.dot(r1);
  double rs2 = r2.dot(r2);

  //------ set reflection flag sgn to avoid branch problems with arctan
  double sgn;
  if (io >= 1 && io <= n) {
    //------- no problem on airfoil surface
    sgn = 1.0;
  } else {
    //------- make sure arctan falls between  -/+  pi/2
    sgn = sign(1.0, yy);
  }
  double logr12, t1;
  //------ set log(r^2) and arctan(x/y), correcting for reflection if any
  if (io != jo && rs1 > 0.0) {
    logr12 = log(rs1);
    t1 = atan2(sgn * x1, sgn * yy) + (0.5 - 0.5 * sgn) * PI;
  } else {
    logr12 = 0.0;
    t1 = 0.0;
  }
  double logr22, t2;
  if (io != jp && rs2 > 0.0) {
    logr22 = log(rs2);
    t2 = atan2(sgn * x2, sgn * yy) + (0.5 - 0.5 * sgn) * PI;
  } else {
    logr22 = 0.0;
    t2 = 0.0;
  }

  double x1i = s.dot(normal_vector);
  double x2i = s.dot(normal_vector);
  double yyi = cross2(s, normal_vector);

  //------- set up midpoint quantities
  double x0 = 0.5 * (x1 + x2);
  double rs0 = x0 * x0 + yy * yy;
  double logr0 = log(rs0);
  double theta0 = atan2(sgn * x0, sgn * yy) + (0.5 - 0.5 * sgn) * PI;

  //------- calculate source contribution to psi	for  1-0  half-panel
  double dxinv = 1.0 / (x1 - x0);
  double psum = x0 * (theta0 - apan) - x1 * (t1 - apan) +
          0.5 * yy * (logr12 - logr0);
  double pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (theta0 - apan) +
          (x0 - x1) * yy) *
          dxinv;

  double psx1 = -(t1- apan);
  double psx0 = theta0 - apan;
  double psyy = 0.5 * (logr12 - logr0);

  double pdx1 =
      ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
  double pdx0 =
      ((x1 + x0) * psx0 + psum - 2.0 * x0 * (theta0 - apan) + pdif) * dxinv;
  double pdyy =
      ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - theta0))) * dxinv;

  const double dsm = (points.col(jp + INDEX_START_WITH) - points.col(jm + INDEX_START_WITH)).norm();
  double dsim = 1.0 / dsm;

  //------- dpsi/dm
  psi_result.dzdm[jm + INDEX_START_WITH] += (1 / (4 * PI)) * (-psum * dsim + pdif * dsim);
  psi_result.dzdm[jo + INDEX_START_WITH] += (1 / (4 * PI)) * (-psum / dso - pdif / dso);
  psi_result.dzdm[jp + INDEX_START_WITH] += (1 / (4 * PI)) * (psum * (dsio + dsim) + pdif * (dsio - dsim));

  //------- dpsi/dni
  double psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
  double pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;

  psi_result.dqdm[jm + INDEX_START_WITH] += (1 / (4 * PI)) * (-psni * dsim + pdni * dsim);
  psi_result.dqdm[jo + INDEX_START_WITH] += (1 / (4 * PI)) * (-psni / dso - pdni / dso);
  psi_result.dqdm[jp + INDEX_START_WITH] += (1 / (4 * PI)) * (psni * (dsio + dsim) + pdni * (dsio - dsim));

  //------- calculate source contribution to psi	for  0-2  half-panel
  dxinv = 1.0 / (x0 - x2);
  psum = x2 * (t2 - apan) - x0 * (theta0 - apan) +
          0.5 * yy * (logr0 - logr22);
  pdif = ((x0 + x2) * psum + rs0 * (theta0 - apan) - rs2 * (t2 - apan) +
          (x2 - x0) * yy) *
          dxinv;

  psx0 = -(theta0 - apan);
  double psx2 = t2 - apan;
  psyy = 0.5 * (logr0 - logr22);

  pdx0 =
      ((x0 + x2) * psx0 + psum + 2.0 * x0 * (theta0 - apan) - pdif) * dxinv;
  double pdx2 =
      ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv;
  pdyy =
      ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (theta0 - t2))) * dxinv;

  double dsp = (points.col(jq) - points.col(jo)).norm();
  double dsip = 1.0 / dsp;

  //------- dpsi/dm
  psi_result.dzdm[jo + INDEX_START_WITH] += (1 / (4 * PI)) * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
  psi_result.dzdm[jp + INDEX_START_WITH] += (1 / (4 * PI)) * (psum / dso - pdif / dso);
  psi_result.dzdm[jq + INDEX_START_WITH] += (1 / (4 * PI)) * (psum * dsip + pdif * dsip);

  //------- dpsi/dni
  psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
  pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;

  psi_result.dqdm[jo + INDEX_START_WITH] += (1 / (4 * PI)) * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
  psi_result.dqdm[jp + INDEX_START_WITH] += (1 / (4 * PI)) * (psni / dso - pdni / dso);
  psi_result.dqdm[jq + INDEX_START_WITH] += (1 / (4 * PI)) * (psni * dsip + pdni * dsip);

  return psi_result;
}

XFoil::PsiResult XFoil::psi_te(int iNode, Vector2d point, Vector2d normal_vector) {
  PsiResult psi_result;
  psi_result.psi = 0;
  psi_result.psi_ni = 0;

  //------ skip null panel

  double apan = apanel[n];

  Vector2d r1 = point - points.col(n);
  Vector2d r2 = point - points.col(1);
  Vector2d s = (points.col(1) - points.col(n)).normalized();

  blData1.param.xz = s.dot(r1);
  blData2.param.xz = s.dot(r2);
  double yy = cross2(s, r1);

  //FIXME normを使うと正しく計算されない。計算精度の問題？
  const double rs1 = r1.dot(r1);
  const double rs2 = r2.dot(r2);

  //------ set reflection flag sgn to avoid branch problems with arctan
  double sgn;
  if (iNode >= 1 && iNode <= n) {
    //------- no problem on airfoil surface
    sgn = 1.0;
  } else {
    //------- make sure arctan falls between  -/+  pi/2
    sgn = sign(1.0, yy);
  }

  //------ set log(r^2) and arctan(x/y), correcting for reflection if any
  double logr12, logr22;
  if (iNode != n && rs1 > 0.0) {
    logr12 = log(rs1);
    blData1.param.tz = atan2(sgn * blData1.param.xz, sgn * yy) + (0.5 - 0.5 * sgn) * PI;
  } else {
    logr12 = 0.0;
    blData1.param.tz = 0.0;
  }
  
  if (iNode != 1 && rs2 > 0.0) {
    logr22 = log(rs2);
    blData2.param.tz = atan2(sgn * blData2.param.xz, sgn * yy) + (0.5 - 0.5 * sgn) * PI;
  } else {
    logr22 = 0.0;
    blData2.param.tz = 0.0;
  }

  double scs, sds;
  if (sharp) {
    scs = 1.0;
    sds = 0.0;
  } else {
    scs = ante / dste;
    sds = aste / dste;
  }

  double x1i = s.dot(normal_vector);
  double x2i = s.dot(normal_vector);
  double yyi = cross2(s, normal_vector);

  double psig = 0.5 * yy * (logr12 - logr22) + blData2.param.xz * (blData2.param.tz - apan) -
         blData1.param.xz * (blData1.param.tz - apan);
  double pgam =
      0.5 * blData1.param.xz * logr12 - 0.5 * blData2.param.xz * logr22 + blData2.param.xz - blData1.param.xz + yy * (blData1.param.tz - blData2.param.tz);

  double psigx1 = -(blData1.param.tz - apan);
  double psigx2 = blData2.param.tz - apan;
  double psigyy = 0.5 * (logr12 - logr22);
  double pgamx1 = 0.5 * logr12;
  double pgamx2 = -0.5 * logr22;
  double pgamyy = blData1.param.tz - blData2.param.tz;

  double psigni = psigx1 * x1i + psigx2 * x2i + psigyy * yyi;
  double pgamni = pgamx1 * x1i + pgamx2 * x2i + pgamyy * yyi;

  //---- TE panel source and vortex strengths
  Vector2d sigte_vector = 0.5 * scs * (gamu.col(0) - gamu.col(n - 1));
  Vector2d gamte_vector = -0.5 * sds * (gamu.col(0) - gamu.col(n - 1));

  sigte = 0.5 * scs * (gam(0, 1) - gam(0, n));
  gamte = -0.5 * sds * (gam(0, 1) - gam(0, n));

  //---- TE panel contribution to psi
  psi_result.psi += (1 / (2 * PI)) * (psig * sigte + pgam * gamte);

  //---- dpsi/dgam
  psi_result.dzdg[n] += -(1 / (2 * PI)) * psig * scs * 0.5;
  psi_result.dzdg[1] += +(1 / (2 * PI)) * psig * scs * 0.5;

  psi_result.dzdg[n] += +(1 / (2 * PI)) * pgam * sds * 0.5;
  psi_result.dzdg[1] += -(1 / (2 * PI)) * pgam * sds * 0.5;

  //---- dpsi/dni
  psi_result.psi_ni += (1 / (2 * PI)) * (psigni * sigte + pgamni * gamte);

  psi_result.qtan += (1 / (2 * PI)) * (psigni * sigte_vector + pgamni * gamte_vector);

  psi_result.dqdg[n] += -(1 / (2 * PI)) * (psigni * 0.5 * scs - pgamni * 0.5 * sds);
  psi_result.dqdg[1] += +(1 / (2 * PI)) * (psigni * 0.5 * scs - pgamni * 0.5 * sds);

  return psi_result;
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
XFoil::PsiResult XFoil::pswlin(int i, Vector2d point, Vector2d normal_vector) {
  PsiResult psi_result;
  double g1, g2, t1, t2;
  int io, jo;

  io = i;
  psi_result.psi = 0.0;
  psi_result.psi_ni = 0.0;

  for (jo = n + 1; jo <= n + nw - 1; jo++) {
    int jp = jo + 1;
    int jm = jo - 1;
    int jq = jp + 1;
    if (jo == n + 1) {
      jm = jo;
    } else {
      if (jo == n + nw - 1) jq = jp;
    }
    const double dso = (points.col(jo) - points.col(jp)).norm() ;
    const double dsio = 1.0 / dso;

    const double apan = apanel[jo];
    
    const Vector2d r1 = point - points.col(jo);
    const Vector2d r2 = point - points.col(jp);

    const Vector2d s = (points.col(jp) - points.col(jo)) * dsio;

    blData1.param.xz = s.dot(r1);
    blData2.param.xz = s.dot(r2);
    const double yy = cross2(s, r1);
    const double rs1 = r1.dot(r1);
    const double rs2 = r2.dot(r2);

    double sgn = 1.0;

    if (io >= n + 1 && io <= n + nw) {
      sgn = 1.0;
    } else {
      sgn = sign(1.0, yy);
    }

    if (io != jo && rs1 > 0.0) {
      g1 = log(rs1);
      t1 = atan2(sgn * blData1.param.xz, sgn * yy) - (0.5 - 0.5 * sgn) * PI;
    } else {
      g1 = 0.0;
      t1 = 0.0;
    }

    if (io != jp && rs2 > 0.0) {
      g2 = log(rs2);
      t2 = atan2(sgn * blData2.param.xz, sgn * yy) - (0.5 - 0.5 * sgn) * PI;
    } else {
      g2 = 0.0;
      t2 = 0.0;
    }
    const double x1i = s.dot(normal_vector);
    const double x2i = s.dot(normal_vector);
    const double yyi = cross2(s, normal_vector);
    //------- set up midpoint quantities
    const double x0 = 0.5 * (blData1.param.xz + blData2.param.xz);
    const double rs0 = x0 * x0 + yy * yy;
    const double g0 = log(rs0);
    const double t0 = atan2(sgn * x0, sgn * yy) - (0.5 - 0.5 * sgn) * PI;

    //------- calculate source contribution to psi	for  1-0  half-panel
    double dxinv = 1.0 / (blData1.param.xz - x0);
    double psum = x0 * (t0 - apan) - blData1.param.xz * (t1 - apan) + 0.5 * yy * (g1 - g0);
    double pdif = ((blData1.param.xz + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan) +
            (x0 - blData1.param.xz) * yy) *
           dxinv;

    double psx1 = -(t1 - apan);
    double psx0 = t0 - apan;
    double psyy = 0.5 * (g1 - g0);

    double pdx1 = ((blData1.param.xz + x0) * psx1 + psum + 2.0 * blData1.param.xz * (t1 - apan) - pdif) * dxinv;
    double pdx0 = ((blData1.param.xz + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
    double pdyy = ((blData1.param.xz + x0) * psyy + 2.0 * (x0 - blData1.param.xz + yy * (t1 - t0))) * dxinv;

    const double dsm = (points.col(jp) - points.col(jm)).norm();
    const double dsim = 1.0 / dsm;

    //------- dpsi/dm
    psi_result.dzdm[jm] += (1 / (4 * PI)) * (-psum * dsim + pdif * dsim);
    psi_result.dzdm[jo] += (1 / (4 * PI)) * (-psum / dso - pdif / dso);
    psi_result.dzdm[jp] += (1 / (4 * PI)) * (psum * (dsio + dsim) + pdif * (dsio - dsim));

    //------- dpsi/dni
    double psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
    double pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;

    psi_result.dqdm[jm] += (1 / (4 * PI)) * (-psni * dsim + pdni * dsim);
    psi_result.dqdm[jo] += (1 / (4 * PI)) * (-psni / dso - pdni / dso);
    psi_result.dqdm[jp] += (1 / (4 * PI)) * (psni * (dsio + dsim) + pdni * (dsio - dsim));

    //------- calculate source contribution to psi	for  0-2  half-panel
    dxinv = 1.0 / (x0 - blData2.param.xz);
    psum = blData2.param.xz * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
    pdif = ((x0 + blData2.param.xz) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan) +
            (blData2.param.xz - x0) * yy) *
           dxinv;

    psx0 = -(t0 - apan);
    const double psx2 = t2 - apan;
    psyy = 0.5 * (g0 - g2);

    pdx0 = ((x0 + blData2.param.xz) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv;
    const double pdx2 = ((x0 + blData2.param.xz) * psx2 + psum - 2.0 * blData2.param.xz * (t2 - apan) + pdif) * dxinv;
    pdyy = ((x0 + blData2.param.xz) * psyy + 2.0 * (blData2.param.xz - x0 + yy * (t0 - t2))) * dxinv;

    const double dsp = (points.col(jq) - points.col(jo)).norm();
    const double dsip = 1.0 / dsp;

    //------- dpsi/dm
    psi_result.dzdm[jo] += (1 / (4 * PI)) * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
    psi_result.dzdm[jp] += (1 / (4 * PI)) * (psum / dso - pdif / dso);
    psi_result.dzdm[jq] += (1 / (4 * PI)) * (psum * dsip + pdif * dsip);

    //------- dpsi/dni
    psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
    pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;

    psi_result.dqdm[jo] += (1 / (4 * PI)) * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
    psi_result.dqdm[jp] += (1 / (4 * PI)) * (psni / dso - pdni / dso);
    psi_result.dqdm[jq] += (1 / (4 * PI)) * (psni * dsip + pdni * dsip);
  }

  return psi_result;
}

/** -----------------------------------------------------
 * 	   calculates source panel influence coefficient
 * 	   matrix for current airfoil and wake geometry.
 * ------------------------------------------------------ */
bool XFoil::qdcalc() {
  double sum;
  VectorXd gamu_temp(n + 1);

  // TRACE("calculating source influence matrix ...\n");
  writeString("   Calculating source influence matrix ...\n");

  if (!ladij) {
    //----- calculate source influence matrix for airfoil surface if it doesn't
    // exist
    for (int j = 0; j < n; j++) {
      //------- multiply each dpsi/sig vector by inverse of factored dpsi/dgam
      // matrix
      for (int iu = 0; iu <= n; iu++) {
        gamu_temp[iu] = bij(iu, j);
      }
      gamu_temp = psi_gamma_lu.solve(gamu_temp);
      for (int iu = 0; iu <= n; iu++) {
        bij(iu, j) = gamu_temp[iu];
      }

      //------- store resulting dgam/dsig = dqtan/dsig vector
      for (int i = 0; i < n; i++) {
        dij[i + INDEX_START_WITH][j + INDEX_START_WITH] = bij(i, j);
      }
    }
    ladij = true;
  }

  //---- set up coefficient matrix of dpsi/dm on airfoil surface
  for (int i = 0; i < n; i++) {
    PsiResult psi_result = pswlin(i, points.col(i + INDEX_START_WITH), normal_vectors.col(i + INDEX_START_WITH));
    for (int j = n; j < n + nw; j++) {
      bij(i, j) = -psi_result.dzdm[j + INDEX_START_WITH];
    }
  }

  //---- set up kutta condition (no direct source influence)
  for (int j = n; j < n + nw; j++) bij(n, j) = 0.0;

  //---- multiply by inverse of factored dpsi/dgam matrix
  for (int j = n; j < n + nw; j++) {
    for (int iu = 0; iu <= n; iu++) {
      gamu_temp[iu] = bij(iu, j);
    }
    gamu_temp = psi_gamma_lu.solve(gamu_temp);

    bij.col(j).head(n + 1) = gamu_temp;
  }
  //---- set the source influence matrix for the wake sources
  for (int i = 0; i < n; i++) {
    for (int j = n; j < n + nw; j++) {
      dij[i + INDEX_START_WITH][j + INDEX_START_WITH] = bij(i, j);
    }
  }

  //**** now we need to calculate the influence of sources on the wake
  // velocities

  //---- calculate dqtan/dgam and dqtan/dsig at the wake points

  for (int i = n; i < n + nw; i++) {
    int iw = i - n;
    //------ airfoil contribution at wake panel node
    PsiResult psi_result = psilin(i + INDEX_START_WITH, points.col(i + INDEX_START_WITH), normal_vectors.col(i + INDEX_START_WITH), true);
    for (int j = 0; j < n; j++) {
      cij[iw + INDEX_START_WITH][j + INDEX_START_WITH] = psi_result.dqdg[j + INDEX_START_WITH];
    }
    for (int j = 0; j < n; j++) {
      dij[i + INDEX_START_WITH][j + INDEX_START_WITH] = psi_result.dqdm[j + INDEX_START_WITH];
    }
    //------ wake contribution
    psi_result = pswlin(i + INDEX_START_WITH, points.col(i + INDEX_START_WITH), normal_vectors.col(i + INDEX_START_WITH));
    for (int j = n; j < n + nw; j++) {
      dij[i + INDEX_START_WITH][j + INDEX_START_WITH] = psi_result.dqdm[j + INDEX_START_WITH];
    }
  }

  //---- add on effect of all sources on airfoil vorticity which effects wake
  // qtan
  for (int i = n; i < n + nw; i++) {
    int iw = i - n;

    //------ airfoil surface source contribution first
    for (int j = 0; j < n; j++) {
      sum = 0.0;
      for (int k = 0; k < n; k++) sum = sum + cij[iw + INDEX_START_WITH][k + INDEX_START_WITH] * dij[k + INDEX_START_WITH][j + INDEX_START_WITH];
      dij[i + INDEX_START_WITH][j + INDEX_START_WITH] = dij[i + INDEX_START_WITH][j + INDEX_START_WITH] + sum;
    }

    //------ wake source contribution next
    for (int j = n; j < n + nw; j++) {
      sum = 0.0;
      for (int k = 0; k < n; k++) sum = sum + cij[iw + INDEX_START_WITH][k + INDEX_START_WITH] * bij(k, j);
      dij[i + INDEX_START_WITH][j + INDEX_START_WITH] = dij[i + INDEX_START_WITH][j + INDEX_START_WITH] + sum;
    }
  }

  //---- make sure first wake point has same velocity as trailing edge
  for (int j = 0; j < n + nw; j++) {
    dij[n + 1][j + INDEX_START_WITH] = dij[n][j + INDEX_START_WITH];
  }

  lwdij = true;
  return true;
}

/** -------------------------------------------------------
 *     sets inviscid panel tangential velocity for
 *      current alpha.
 * -------------------------------------------------------- */
bool XFoil::qiset() {
  Matrix2d rotateMatrix = Matrix2d {
    {cos(alfa), sin(alfa)},
    {-sin(alfa), cos(alfa)}
  };
  
  for (int i = 1; i <= n + nw; i++) {
    qinv[i] = rotateMatrix.row(0).dot(qinvu.col(i));
    qinv_a[i] = rotateMatrix.row(1).dot(qinvu.col(i));
  }

  return true;
}

/** -------------------------------------------------------------
 *     sets panel viscous tangential velocity from viscous ue
 * -------------------------------------------------------------- */
bool XFoil::qvfue() {
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 2; ibl <= nbl.get(is); ibl++) {
      int i = ipan.get(is)[ibl];
      qvis[i] = vti.get(is)[ibl] * uedg.get(is)[ibl];
    }
  }

  return true;
}

/** ---------------------------------------------------------------
 *      sets inviscid tangential velocity for alpha = 0, 90
 *      on wake due to freestream and airfoil surface vorticity.
 * --------------------------------------------------------------- */
bool XFoil::qwcalc() {

  //---- first wake point (same as te)
  qinvu.col(n + 1) = qinvu.col(n);

  //---- rest of wake
  for (int i = n + 2; i <= n + nw; i++) {
    qinvu.col(i) = psilin(i, points.col(i), normal_vectors.col(i), false).qtan;
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
  //int i, ibl, iv, iw, j, js = 0, jv, jbl, is = 0;
  int ile1 = 0, ile2 = 0, ite1 = 0, ite2 = 0, jvte1 = 0, jvte2 = 0;
  double u1_m[2 * IVX + 1], u2_m[2 * IVX + 1];
  double d1_m[2 * IVX + 1], d2_m[2 * IVX + 1];
  double ule1_m[2 * IVX + 1], ule2_m[2 * IVX + 1];
  double ute1_m[2 * IVX + 1], ute2_m[2 * IVX + 1];

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

    mrchue();
    lblini = true;
  }

  //---- march bl with current ue and ds to establish transition
  mrchdu();

  SidePair<VectorXd> usav = uedg;

  ueset();

  for (int is = 1; is <= 2; is++) {
    for (int ibl = 2; ibl <= nbl.get(is); ibl++) {
      double temp = usav.get(is)[ibl];
      usav.get(is)[ibl] = uedg.get(is)[ibl];
      uedg.get(is)[ibl] = temp;
    }
  }
  ile1 = ipan.top[2];
  ile2 = ipan.bottom[2];
  ite1 = ipan.top[iblte.top];
  ite2 = ipan.bottom[iblte.bottom];

  jvte1 = isys.top[iblte.top];
  jvte2 = isys.bottom[iblte.bottom];

  dule1 = uedg.top[2] - usav.top[2];
  dule2 = uedg.bottom[2] - usav.bottom[2];

  //---- set le and te ue sensitivities wrt all m values
  for (int js = 1; js <= 2; js++) {
    for (int jbl = 2; jbl <= nbl.get(js); jbl++) {
      int j = ipan.get(js)[jbl];
      int jv = isys.get(js)[jbl];
      ule1_m[jv] = -vti.top[2] * vti.get(js)[jbl] * dij[ile1][j];
      ule2_m[jv] = -vti.bottom[2] * vti.get(js)[jbl] * dij[ile2][j];
      ute1_m[jv] = -vti.top[iblte.top] * vti.get(js)[jbl] * dij[ite1][j];
      ute2_m[jv] = -vti.bottom[iblte.bottom] * vti.get(js)[jbl] * dij[ite2][j];
    }
  }

  ule1_a = uinv_a.top[2];
  ule2_a = uinv_a.bottom[2];

  writeString(" \n");

  //*** go over each boundary layer/wake
  for (int is = 1; is <= 2; is++) {
    //---- there is no station "1" at similarity, so zero everything out
    for (int js = 1; js <= 2; js++) {
      for (int jbl = 2; jbl <= nbl.get(js); jbl++) {
        int jv = isys.get(js)[jbl];
        u1_m[jv] = 0.0;
        d1_m[jv] = 0.0;
      }
    }
    double u1_a = 0.0;
    double d1_a = 0.0;

    double due1 = 0.0;
    double dds1 = 0.0;

    //---- similarity station pressure gradient parameter  x/u du/dx
    bule = 1.0;

    //---- set forced transition arc length position
    xiforc = xifset(is);

    tran = false;
    turb = false;

    //**** sweep downstream setting up bl equation linearizations
    for (int ibl = 2; ibl <= nbl.get(is); ibl++) {
      int iv = isys.get(is)[ibl];

      simi = (ibl == 2);
      wake = (ibl > iblte.get(is));
      tran = (ibl == itran.get(is));
      turb = (ibl > itran.get(is));

      int i = ipan.get(is)[ibl];

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
        dswaki = wgap[iw];
      } else
        dswaki = 0.0;

      //---- set derivatives of dsi (= d2)
      d2_m2 = 1.0 / uei;
      d2_u2 = -dsi / uei;

      for (int js = 1; js <= 2; js++) {
        for (int jbl = 2; jbl <= nbl.get(js); jbl++) {
          int j = ipan.get(js)[jbl];
          int jv = isys.get(js)[jbl];
          u2_m[jv] = -vti.get(is)[ibl] * vti.get(js)[jbl] * dij[i][j];
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

      blprv(xsi, ami, cti, thi, dsi, dswaki, uei);  // cti
      blkin();

      //---- check for transition and set tran, xt, etc. if found
      if (tran) {
        trchek();
        ami = blData2.param.amplz;
      }

      if (ibl == itran.get(is) && !tran) {
        // TRACE("setbl: xtr???  n1=%d n2=%d: \n", ampl1, ampl2);

        ss << "setbl: xtr???  n1=" << blData1.param.amplz << " n2=" << blData2.param.amplz << ":\n";
        writeString(ss.str());
        ss.str("");
      }

      //---- assemble 10x4 linearized system for dctau, dth, dds, due, dxi
      //	   at the previous "1" station and the current "2" station

      if (ibl == iblte.get(is) + 1) {
        //----- define quantities at start of wake, adding te base thickness to
        // dstar
        tte = thet.top[iblte.top] + thet.bottom[iblte.bottom];
        dte = dstr.top[iblte.top] + dstr.bottom[iblte.bottom] + ante;
        cte = (ctau.top[iblte.top] * thet.top[iblte.top] +
               ctau.bottom[iblte.bottom] * thet.bottom[iblte.bottom]) /
              tte;
        tesys(cte, tte, dte);

        tte_tte1 = 1.0;
        tte_tte2 = 1.0;
        dte_mte1 = 1.0 / uedg.top[iblte.top];
        dte_ute1 = -dstr.top[iblte.top] / uedg.top[iblte.top];
        dte_mte2 = 1.0 / uedg.bottom[iblte.bottom];
        dte_ute2 = -dstr.bottom[iblte.bottom] / uedg.bottom[iblte.bottom];
        cte_cte1 = thet.top[iblte.top] / tte;
        cte_cte2 = thet.bottom[iblte.bottom] / tte;
        cte_tte1 = (ctau.top[iblte.top] - cte) / tte;
        cte_tte2 = (ctau.bottom[iblte.bottom] - cte) / tte;

        //----- re-define d1 sensitivities wrt m since d1 depends on both te ds
        // values
        for (int js = 1; js <= 2; js++) {
          for (int jbl = 2; jbl <= nbl.get(js); jbl++) {            
            int jv = isys.get(js)[jbl];
            d1_m[jv] = dte_ute1 * ute1_m[jv] + dte_ute2 * ute2_m[jv];
          }
        }
        d1_m[jvte1] = d1_m[jvte1] + dte_mte1;
        d1_m[jvte2] = d1_m[jvte2] + dte_mte2;

        //----- "forced" changes from  uedg --- usav=uinv+dij*mass	mismatch
        due1 = 0.0;
        dds1 = dte_ute1 * (uedg.top[iblte.top] - usav.top[iblte.top]) +
               dte_ute2 * (uedg.bottom[iblte.bottom] - usav.bottom[iblte.bottom]);
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
        vm[0][jv][iv] = vs1(0, 2) * d1_m[jv] + vs1(0, 3) * u1_m[jv] +
                        vs2(0, 2) * d2_m[jv] + vs2(0, 3) * u2_m[jv] +
                        (vs1(0, 4) + vs2(0, 4) + vsx[0]) *
                            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
      }

      vb[0][0][iv] = vs1(0, 0);
      vb[0][1][iv] = vs1(0, 1);

      va[0][0][iv] = vs2(0, 0);
      va[0][1][iv] = vs2(0, 1);

      if (lalfa)
        vdel[0][1][iv] = vsr[0] * re_clmr + vsm[0] * msq_clmr;
      else
        vdel[0][1][iv] = (vs1(0, 3) * u1_a + vs1(0, 2) * d1_a) +
                         (vs2(0, 3) * u2_a + vs2(0, 2) * d2_a) +
                         (vs1(0, 4) + vs2(0, 4) + vsx[0]) *
                             (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

      vdel[0][0][iv] = vsrez[0] + (vs1(0, 3) * due1 + vs1(0, 2) * dds1) +
                       (vs2(0, 3) * due2 + vs2(0, 2) * dds2) +
                       (vs1(0, 4) + vs2(0, 4) + vsx[0]) *
                           (xi_ule1 * dule1 + xi_ule2 * dule2);

      for (int jv = 1; jv <= nsys; jv++) {
        vm[1][jv][iv] = vs1(1, 2) * d1_m[jv] + vs1(1, 3) * u1_m[jv] +
                        vs2(1, 2) * d2_m[jv] + vs2(1, 3) * u2_m[jv] +
                        (vs1(1, 4) + vs2(1, 4) + vsx[1]) *
                            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
      }
      vb[1][0][iv] = vs1(1, 0);
      vb[1][1][iv] = vs1(1, 1);

      va[1][0][iv] = vs2(1, 0);
      va[1][1][iv] = vs2(1, 1);

      if (lalfa)
        vdel[1][1][iv] = vsr[1] * re_clmr + vsm[1] * msq_clmr;
      else
        vdel[1][1][iv] = (vs1(1, 3) * u1_a + vs1(1, 2) * d1_a) +
                         (vs2(1, 3) * u2_a + vs2(1, 2) * d2_a) +
                         (vs1(1, 4) + vs2(1, 4) + vsx[1]) *
                             (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

      vdel[1][0][iv] = vsrez[1] + (vs1(1, 3) * due1 + vs1(1, 2) * dds1) +
                       (vs2(1, 3) * due2 + vs2(1, 2) * dds2) +
                       (vs1(1, 4) + vs2(1, 4) + vsx[1]) *
                           (xi_ule1 * dule1 + xi_ule2 * dule2);

      // memory overlap problem
      for (int jv = 1; jv <= nsys; jv++) {
        vm[2][jv][iv] = vs1(2, 2) * d1_m[jv] + vs1(2, 3) * u1_m[jv] +
                        vs2(2, 2) * d2_m[jv] + vs2(2, 3) * u2_m[jv] +
                        (vs1(2, 4) + vs2(2, 4) + vsx[2]) *
                            (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
      }

      vb[2][0][iv] = vs1(2, 0);
      vb[2][1][iv] = vs1(2, 1);

      va[2][0][iv] = vs2(2, 0);
      va[2][1][iv] = vs2(2, 1);

      if (lalfa)
        vdel[2][1][iv] = vsr[2] * re_clmr + vsm[2] * msq_clmr;
      else
        vdel[2][1][iv] = (vs1(2, 3) * u1_a + vs1(2, 2) * d1_a) +
                         (vs2(2, 3) * u2_a + vs2(2, 2) * d2_a) +
                         (vs1(2, 4) + vs2(2, 4) + vsx[2]) *
                             (xi_ule1 * ule1_a + xi_ule2 * ule2_a);

      vdel[2][0][iv] = vsrez[2] + (vs1(2, 3) * due1 + vs1(2, 2) * dds1) +
                       (vs2(2, 3) * due2 + vs2(2, 2) * dds2) +
                       (vs1(2, 4) + vs2(2, 4) + vsx[2]) *
                           (xi_ule1 * dule1 + xi_ule2 * dule2);

      if (ibl == iblte.get(is) + 1) {
        //----- redefine coefficients for tte, dte, etc
        vz[0][0] = vs1(0, 0) * cte_cte1;
        vz[0][1] = vs1(0, 0) * cte_tte1 + vs1(0, 1) * tte_tte1;
        vb[0][0][iv] = vs1(0, 0) * cte_cte2;
        vb[0][1][iv] = vs1(0, 0) * cte_tte2 + vs1(0, 1) * tte_tte2;

        vz[1][0] = vs1(1, 0) * cte_cte1;
        vz[1][1] = vs1(1, 0) * cte_tte1 + vs1(1, 1) * tte_tte1;
        vb[1][0][iv] = vs1(1, 0) * cte_cte2;
        vb[1][1][iv] = vs1(1, 0) * cte_tte2 + vs1(1, 1) * tte_tte2;

        vz[2][0] = vs1(2, 0) * cte_cte1;
        vz[2][1] = vs1(2, 0) * cte_tte1 + vs1(2, 1) * tte_tte1;
        vb[2][0][iv] = vs1(2, 0) * cte_cte2;
        vb[2][1][iv] = vs1(2, 0) * cte_tte2 + vs1(2, 1) * tte_tte2;
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
        blvar(blData2, 3);
        blmid(3);
      }

      for (int js = 1; js <= 2; js++) {
        for (int jbl = 2; jbl <= nbl.get(js); jbl++) {
          int jv = isys.get(js)[jbl];
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
  int nex;
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
    writeString("setexp: cannot fill array.  n too small\n");
    return false;
  } else {
    if (nex == 2)
      ratio = -ccc / bbb + 1.0;
    else
      ratio = (-bbb + sqrt(disc)) / (2.0 * aaa) + 1.0;
  }

  //-- newton iteration for actual geometric ratio
  for (int iter = 1; iter <= 100; iter++) {
    sigman = (pow(ratio, (double)nex) - 1.0) / (ratio - 1.0);
    res = pow(sigman, rni) - pow(sigma, rni);
    dresdr = rni * pow(sigman, rni) *
             (rnex * pow(ratio, (double)(nex - 1)) - sigman) /
             (pow(ratio, (double)nex) - 1.0);

    dratio = -res / dresdr;
    ratio = ratio + dratio;

    if (fabs(dratio) < 1.0e-5) {
      break;
    }
  }

  spline_length[1] = 0.0;
  ds = ds1;
  for (int i = 2; i <= nn; i++) {
    spline_length[i] = spline_length[i - 1] + ds;
    ds = ds * ratio;
  }
  return true;
}

bool XFoil::setMach() {
  minf_cl = getActualMach(1.0, mach_type);
  reinf_cl = getActualReynolds(1.0, reynolds_type);
  comset();
  cpi = cpcalc(n, qinv, qinf, minf);
  if (lvisc) {
    cpv = cpcalc(n + nw, qvis, qinf, minf);
  }
  clcalc(cmref);
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
 *      Converges to specified alpha.
 */
bool XFoil::specal() {
  double minf_clm;
  double clm;

  //---- calculate surface vorticity distributions for alpha = 0, 90 degrees
  if (!lgamu || !lqaij) ggcalc();

  Matrix2d rotateMatrix = Matrix2d {
    {cos(alfa), sin(alfa)},
    {-sin(alfa), cos(alfa)}
  };

  //---- superimpose suitably weighted  alpha = 0, 90  distributions
  for (int i = 1; i <= n; i++) {
    gam(0, i) = rotateMatrix.row(0).dot(gamu.col(i - INDEX_START_WITH));
    gam(1, i) = rotateMatrix.row(1).dot(gamu.col(i - INDEX_START_WITH));
  }

  tecalc();
  qiset();

  //---- set initial guess for the newton variable clm
  clm = 1.0;

  //---- set corresponding  m(clm), re(clm)
  minf_clm = getActualMach(clm, mach_type);

  comset();

  //---- set corresponding cl(m)
  clcalc(cmref);
  //---- iterate on clm
  bool bConv = false;
  for (int itcl = 1; itcl <= 20; itcl++) {
    const double msq_clm = 2.0 * minf * minf_clm;
    const double dclm = (cl - clm) / (1.0 - cl_msq * msq_clm);

    const double clm1 = clm;
    rlx = 1.0;

    //------ under-relaxation loop to avoid driving m(cl) above 1
    for (int irlx = 1; irlx <= 12; irlx++) {
      clm = clm1 + rlx * dclm;

      //-------- set new freestream mach m(clm)
      minf_clm = getActualMach(clm, mach_type);

      //-------- if mach is ok, go do next newton iteration
      //FIXME double型の==比較
      if (mach_type == MachType::CONSTANT || minf == 0.0 || minf_clm != 0.0) break;

      rlx = 0.5 * rlx;
    }

    //------ set new cl(m)
    comset();
    clcalc(cmref);

    if (fabs(dclm) <= 1.0e-6) {
      bConv = true;
      break;
    }
  }
  if (!bConv) {
    writeString("Specal:  MInf convergence failed\n");
    return false;
  }

  //---- set final mach, cl, cp distributions, and hinge moment
  minf_cl = getActualMach(cl, mach_type);
  reinf_cl = getActualReynolds(cl, reynolds_type);
  comset();
  clcalc(cmref);

  cpi = cpcalc(n, qinv, qinf, minf);
  if (lvisc) {
    cpv = cpcalc(n + nw, qvis, qinf, minf);
    cpi = cpcalc(n + nw, qinv, qinf, minf);
  } else
    cpi = cpcalc(n, qinv, qinf, minf);

  for (int i = 1; i <= n; i++) {
    qgamm[i] = gam(0, i);
  }

  return true;
}

bool XFoil::speccl() {
  //-----------------------------------------
  //     converges to specified inviscid cl.
  //-----------------------------------------

  //---- calculate surface vorticity distributions for alpha = 0, 90 degrees
  if (!lgamu || !lqaij) ggcalc();

  //---- set freestream mach from specified cl -- mach will be held fixed
  minf_cl = getActualMach(clspec, mach_type);
  reinf_cl = getActualReynolds(clspec, reynolds_type);
  comset();

  Matrix2d rotateMatrix = Matrix2d {
    {cos(alfa), sin(alfa)},
    {-sin(alfa), cos(alfa)}
  };

  //---- superimpose suitably weighted  alpha = 0, 90  distributions
  for (int i = 1; i <= n; i++) {
    gam(0, i) = rotateMatrix.row(0).dot(gamu.col(i - INDEX_START_WITH));
    gam(1, i) = rotateMatrix.row(1).dot(gamu.col(i - INDEX_START_WITH));
  }

  //---- get corresponding cl, cl_alpha, cl_mach
  clcalc(cmref);

  //---- newton loop for alpha to get specified inviscid cl
  bool bConv = false;
  for (int ital = 1; ital <= 20; ital++) {
    const double dalfa = (clspec - cl) / cl_alf;
    rlx = 1.0;

    alfa = alfa + rlx * dalfa;

    //------ set new surface speed distribution
    Matrix2d rotateMatrix = Matrix2d {
      {cos(alfa), sin(alfa)},
      {-sin(alfa), cos(alfa)}
    };
    for (int i = 1; i <= n; i++) {
      gam(0, i) = rotateMatrix.row(0).dot(gamu.col(i - INDEX_START_WITH));
      gam(1, i) = rotateMatrix.row(1).dot(gamu.col(i - INDEX_START_WITH));
    }

    //------ set new cl(alpha)
    clcalc(cmref);

    if (fabs(dalfa) <= 1.0e-6) {
      bConv = true;
      break;
    }
  }
  if (!bConv) {
    writeString("Speccl:  cl convergence failed");
    return false;
  }

  //---- set final surface speed and cp distributions
  tecalc();
  qiset();

  if (lvisc) {
    cpv = cpcalc(n + nw, qvis, qinf, minf);
    cpi = cpcalc(n + nw, qinv, qinf, minf);

  } else {
    cpi = cpcalc(n, qinv, qinf, minf);
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
    if (gam(0, i) >= 0.0 && gam(0, i + 1) < 0.0) {
      bFound = true;
      break;
    }
  }

  if (!bFound) {
    writeString("stfind: Stagnation point not found. Continuing ...\n");
    i = n / 2;
  }

  ist = i;
  dgam = gam(0, i + 1) - gam(0, i);
  ds = spline_length[i + 1] - spline_length[i];

  //---- evaluate so as to minimize roundoff for very small gam[i] or gam[i+1]
  if (gam(0, i) < -gam(0, i + 1))
    sst = spline_length[i] - ds * (gam(0, i) / dgam);
  else
    sst = spline_length[i + 1] - ds * (gam(0, i + 1) / dgam);

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
  int istold;
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

      itran.top = itran.top + idif;
      itran.bottom = itran.bottom - idif;

      //---- move top side bl variables downstream
      for (int ibl = nbl.top; ibl >= idif + 2; ibl--) {
        ctau.top[ibl] = ctau.top[ibl - idif];
        thet.top[ibl] = thet.top[ibl - idif];
        dstr.top[ibl] = dstr.top[ibl - idif];
        uedg.top[ibl] = uedg.top[ibl - idif];
      }

      //---- set bl variables between old and new stagnation point
      const double dudx = uedg.top[idif + 2] / xssi.top[idif + 2];
      for (int ibl = idif + 1; ibl >= 2; ibl--) {
        ctau.top[ibl] = ctau.top[idif + 2];
        thet.top[ibl] = thet.top[idif + 2];
        dstr.top[ibl] = dstr.top[idif + 2];
        uedg.top[ibl] = dudx * xssi.top[ibl];
      }

      //---- move bottom side bl variables upstream
      for (int ibl = 2; ibl <= nbl.bottom; ibl++) {
        ctau.bottom[ibl] = ctau.bottom[ibl + idif];
        thet.bottom[ibl] = thet.bottom[ibl + idif];
        dstr.bottom[ibl] = dstr.bottom[ibl + idif];
        uedg.bottom[ibl] = uedg.bottom[ibl + idif];
      }
    } else {
      //---- increase in number of points on bottom side (is=2)
      int idif = istold - ist;

      itran.top = itran.top - idif;
      itran.bottom = itran.bottom + idif;

      //---- move bottom side bl variables downstream
      for (int ibl = nbl.bottom; ibl >= idif + 2; ibl--) {
        ctau.bottom[ibl] = ctau.bottom[ibl - idif];
        thet.bottom[ibl] = thet.bottom[ibl - idif];
        dstr.bottom[ibl] = dstr.bottom[ibl - idif];
        uedg.bottom[ibl] = uedg.bottom[ibl - idif];
      }

      //---- set bl variables between old and new stagnation point
      const double dudx = uedg.bottom[idif + 2] / xssi.bottom[idif + 2];
      for (int ibl = idif + 1; ibl >= 2; ibl--) {
        ctau.bottom[ibl] = ctau.bottom[idif + 2];
        thet.bottom[ibl] = thet.bottom[idif + 2];
        dstr.bottom[ibl] = dstr.bottom[idif + 2];
        uedg.bottom[ibl] = dudx * xssi.bottom[ibl];
      }

      //---- move top side bl variables upstream
      for (int ibl = 2; ibl <= nbl.top; ibl++) {
        ctau.top[ibl] = ctau.top[ibl + idif];
        thet.top[ibl] = thet.top[ibl + idif];
        dstr.top[ibl] = dstr.top[ibl + idif];
        uedg.top[ibl] = uedg.top[ibl + idif];
      }
    }
  }

  //-- set new mass array since ue has been tweaked
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 2; ibl <= nbl.get(is); ibl++)
      mass.get(is)[ibl] = dstr.get(is)[ibl] * uedg.get(is)[ibl];
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
  Vector2d point_te = points.col(1) - points.col(n);
  
  Vector2d dpoint_ds_te = 0.5 * (-dpoints_ds.col(1) + dpoints_ds.col(n));

  //---- normal and streamwise projected TE gap areas
  ante = cross2(dpoint_ds_te, point_te);
  aste = point_te.dot(dpoint_ds_te);

  //---- total TE gap area
  dste = point_te.norm();

  sharp = dste < 0.0001 * chord;

  if (sharp) {
    scs = 1.0;
    sds = 0.0;
  } else {
    scs = ante / dste;
    sds = aste / dste;
  }

  //---- TE panel source and vorticity strengths
  sigte = 0.5 * (gam(0, 1) - gam(0, n)) * scs;
  gamte = -.5 * (gam(0, 1) - gam(0, n)) * sds;

  return true;
}

bool XFoil::tesys(double cte, double tte, double dte) {
  //--------------------------------------------------------
  //	   sets up "dummy" bl system between airfoil te point
  //	   and first wake point infinitesimally behind te.
  //--------------------------------------------------------

  vsrez = Vector<double, 4>::Zero();
  vsm = Vector<double, 4>::Zero();
  vsr = Vector<double, 4>::Zero();
  vsx = Vector<double, 4>::Zero();
  vs1 = Matrix<double, 4, 5>::Zero();
  vs2 = Matrix<double, 4, 5>::Zero();

  blvar(blData2, 3);

  vs1(0, 0) = -1.0;
  vs2(0, 0) = 1.0;
  vsrez[0] = cte - blData2.param.sz;

  vs1(1, 1) = -1.0;
  vs2(1, 1) = 1.0;
  vsrez[1] = tte - blData2.param.tz;

  vs1(2, 2) = -1.0;
  vs2(2, 2) = 1.0;
  vsrez[2] = dte - blData2.param.dz - blData2.param.dwz;

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
  double ut_a2, hkt, hkt_tt, hkt_dt, hkt_ut, hkt_ms, rtt_tt, rtt_ut, rtt_ms,
      rtt, rtt_re;
  double daeps = 0.00005;

  amplt_a2 = 0.0;
  xt_a2 = dt_a2 = tt_a2 = 0.0;
  ut_a2 = hkt_tt = hkt_dt = hkt_ut = hkt_ms = rtt_tt = rtt_ut = rtt_ms = rtt_re = 0.0;

  //---- save variables and sensitivities at ibl ("2") for future restoration
  saveblData(2);

  //---- calculate average amplification rate ax over x1..x2 interval
  AxResult ax_result = axset(blData1.hkz.scalar, blData1.param.tz, blData1.rtz.scalar, blData1.param.amplz, blData2.hkz.scalar, blData2.param.tz, blData2.rtz.scalar, blData2.param.amplz, amcrit);

  //---- set initial guess for iterate n2 (ampl2) at x2
  blData2.param.amplz = blData1.param.amplz + ax_result.ax * (blData2.param.xz - blData1.param.xz);

  //---- solve implicit system for amplification ampl2
  for (int itam = 0; itam < 30; itam++) {
    //---- define weighting factors wf1,wf2 for defining "t" quantities from 1,2
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
      sfa = (amplt - blData1.param.amplz) / (blData2.param.amplz - blData1.param.amplz);
      sfa_a1 = (sfa - 1.0) / (blData2.param.amplz - blData1.param.amplz);
      sfa_a2 = (-sfa) / (blData2.param.amplz - blData1.param.amplz);
    }

    if (xiforc < blData2.param.xz) {
      sfx = (xiforc - blData1.param.xz) / (blData2.param.xz - blData1.param.xz);
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
    blkin();

    

    hkt = blData2.hkz.scalar;
    hkt_tt = blData2.hkz.t();
    hkt_dt = blData2.hkz.d();
    hkt_ut = blData2.hkz.u();
    hkt_ms = blData2.hkz.ms();

    rtt = blData2.rtz.scalar;
    rtt_tt = blData2.rtz.t();
    rtt_ut = blData2.rtz.u();
    rtt_ms = blData2.rtz.ms();
    rtt_re = blData2.rtz.re();

    //---- restore clobbered "2" variables, except for ampl2
    amsave = blData2.param.amplz;

    restoreblData(2);

    blData2.param.amplz = amsave;

    //---- calculate amplification rate ax over current x1-xt interval
    ax_result = axset(blData1.hkz.scalar, blData1.param.tz, blData1.rtz.scalar, blData1.param.amplz, hkt, tt, rtt, amplt, amcrit);

    //---- punch out early if there is no amplification here
    if (ax_result.ax <= 0.0) {
      break;
    }

    //---- set sensitivity of ax(a2)
    ax_result.ax_a2 = (ax_result.ax_hk2 * hkt_tt + ax_result.ax_t2 + ax_result.ax_rt2 * rtt_tt) * tt_a2 +
            (ax_result.ax_hk2 * hkt_dt) * dt_a2 +
            (ax_result.ax_hk2 * hkt_ut + ax_result.ax_rt2 * rtt_ut) * ut_a2 + ax_result.ax_a2 * amplt_a2;

    //---- residual for implicit ampl2 definition (amplification equation)
    res = blData2.param.amplz - blData1.param.amplz - ax_result.ax * (blData2.param.xz - blData1.param.xz);
    res_a2 = 1.0 -  ax_result.ax_a2 * (blData2.param.xz - blData1.param.xz);

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
      break;
    }

    if ((blData2.param.amplz > amcrit && blData2.param.amplz + rlx * da2 < amcrit) ||
        (blData2.param.amplz < amcrit && blData2.param.amplz + rlx * da2 > amcrit)) {
      //------ limited newton step so ampl2 doesn't step across amcrit either
      // way
      blData2.param.amplz = amcrit;
    }
    else {
      //------ regular newton step
      blData2.param.amplz = blData2.param.amplz + rlx * da2;
    }
  }

  // TRACE("trchek2 - n2 convergence failed\n");
  writeString("trchek2 - n2 convergence failed\n");
  if (s_bCancel) return false;

  //---- test for free or forced transition
  trfree = (blData2.param.amplz >= amcrit);
  trforc = (xiforc > blData1.param.xz) && (xiforc <= blData2.param.xz);

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

  //---- set sensitivities of ax( t1 d1 u1 a1 t2 d2 u2 a2 ms re )
  double ax_t1 = ax_result.ax_hk1 * blData1.hkz.t() + ax_result.ax_t1 + ax_result.ax_rt1 * blData1.rtz.t() +
          (ax_result.ax_hk2 * hkt_tt + ax_result.ax_t2 + ax_result.ax_rt2 * rtt_tt) * tt_t1;
  double ax_d1 = ax_result.ax_hk1 * blData1.hkz.d() + (ax_result.ax_hk2 * hkt_dt) * dt_d1;
  double ax_u1 = ax_result.ax_hk1 * blData1.hkz.u() + ax_result.ax_rt1 * blData1.rtz.u() +
          (ax_result.ax_hk2 * hkt_ut + ax_result.ax_rt2 * rtt_ut) * ut_u1;
  double ax_a1 = ax_result.ax_a1 +
          (ax_result.ax_hk2 * hkt_tt + ax_result.ax_t2 + ax_result.ax_rt2 * rtt_tt) *
              tt_a1
          + (ax_result.ax_hk2 * hkt_dt) * dt_a1 +
          (ax_result.ax_hk2 * hkt_ut + ax_result.ax_rt2 * rtt_ut) * ut_a1;
  double ax_x1 = (ax_result.ax_hk2 * hkt_tt + ax_result.ax_t2 + ax_result.ax_rt2 * rtt_tt) * tt_x1 +
          (ax_result.ax_hk2 * hkt_dt) * dt_x1 +
          (ax_result.ax_hk2 * hkt_ut + ax_result.ax_rt2 * rtt_ut) * ut_x1;

  double ax_t2 = (ax_result.ax_hk2 * hkt_tt + ax_result.ax_t2 + ax_result.ax_rt2 * rtt_tt) * tt_t2;
  double ax_d2 = (ax_result.ax_hk2 * hkt_dt) * dt_d2;
  double ax_u2 = (ax_result.ax_hk2 * hkt_ut + ax_result.ax_rt2 * rtt_ut) * ut_u2;
  double ax_a2 = ax_result.ax_a2 * amplt_a2 +
          (ax_result.ax_hk2 * hkt_tt + ax_result.ax_t2 + ax_result.ax_rt2 * rtt_tt) *
              tt_a2 
          + (ax_result.ax_hk2 * hkt_dt) * dt_a2 +
          (ax_result.ax_hk2 * hkt_ut + ax_result.ax_rt2 * rtt_ut) * ut_a2;
  double ax_x2 = (ax_result.ax_hk2 * hkt_tt + ax_result.ax_t2 + ax_result.ax_rt2 * rtt_tt) * tt_x2 +
          (ax_result.ax_hk2 * hkt_dt) * dt_x2 +
          (ax_result.ax_hk2 * hkt_ut + ax_result.ax_rt2 * rtt_ut) * ut_x2;

  double ax_ms = ax_result.ax_hk2 * hkt_ms + ax_result.ax_rt2 * rtt_ms + ax_result.ax_hk1 * blData1.hkz.ms() + ax_result.ax_rt1 * blData1.rtz.ms();
  double ax_re = ax_result.ax_rt2 * rtt_re + ax_result.ax_rt1 * blData1.rtz.re();

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
  double wf = (xt - blData1.param.xz) / (blData2.param.xz - blData1.param.xz);
  double wf_xt = 1.0 / (blData2.param.xz - blData1.param.xz);

  double wf_a1 = wf_xt * xt_a1;
  double wf_x1 = wf_xt * xt_x1 + (wf - 1.0) / (blData2.param.xz - blData1.param.xz);
  double wf_x2 = wf_xt * xt_x2 - wf / (blData2.param.xz - blData1.param.xz);
  double wf_t1 = wf_xt * xt_t1;
  double wf_t2 = wf_xt * xt_t2;
  double wf_d1 = wf_xt * xt_d1;
  double wf_d2 = wf_xt * xt_d2;
  double wf_u1 = wf_xt * xt_u1;
  double wf_u2 = wf_xt * xt_u2;
  double wf_ms = wf_xt * xt_ms;
  double wf_re = wf_xt * xt_re;
  double wf_xf = wf_xt * xt_xf;

  //**** first,  do laminar part between x1 and xt

  //-----interpolate primary variables to transition point
  tt = blData1.param.tz * (1 - wf) + blData2.param.tz * wf;
  tt_a1 = (blData2.param.tz - blData1.param.tz) * wf_a1;
  tt_x1 = (blData2.param.tz - blData1.param.tz) * wf_x1;
  tt_x2 = (blData2.param.tz - blData1.param.tz) * wf_x2;
  tt_t1 = (blData2.param.tz - blData1.param.tz) * wf_t1 + (1 - wf);
  tt_t2 = (blData2.param.tz - blData1.param.tz) * wf_t2 + wf;
  tt_d1 = (blData2.param.tz - blData1.param.tz) * wf_d1;
  tt_d2 = (blData2.param.tz - blData1.param.tz) * wf_d2 ;
  tt_u1 = (blData2.param.tz - blData1.param.tz) * wf_u1;
  tt_u2 = (blData2.param.tz - blData1.param.tz) * wf_u2;
  tt_ms = (blData2.param.tz - blData1.param.tz) * wf_ms;
  tt_re = (blData2.param.tz - blData1.param.tz) * wf_re;
  tt_xf = (blData2.param.tz - blData1.param.tz) * wf_xf;

  dt = blData1.param.dz * (1 - wf) + blData2.param.dz * wf;
  dt_a1 = (blData2.param.dz - blData1.param.dz) * wf_a1;
  dt_x1 = (blData2.param.dz - blData1.param.dz) * wf_x1;
  dt_x2 = (blData2.param.dz - blData1.param.dz) * wf_x2;
  dt_t1 = (blData2.param.dz - blData1.param.dz) * wf_t1;
  dt_t2 = (blData2.param.dz - blData1.param.dz) * wf_t2;
  dt_d1 = (blData2.param.dz - blData1.param.dz) * wf_d1 + (1 - wf);
  dt_d2 = (blData2.param.dz - blData1.param.dz) * wf_d2 + wf;
  dt_u1 = (blData2.param.dz - blData1.param.dz) * wf_u1;
  dt_u2 = (blData2.param.dz - blData1.param.dz) * wf_u2;
  dt_ms = (blData2.param.dz - blData1.param.dz) * wf_ms;
  dt_re = (blData2.param.dz - blData1.param.dz) * wf_re;
  dt_xf = (blData2.param.dz - blData1.param.dz) * wf_xf;

  ut = blData1.param.uz * (1 - wf) + blData2.param.uz * wf;
  ut_a1 = (blData2.param.uz - blData1.param.uz) * wf_a1;
  ut_x1 = (blData2.param.uz - blData1.param.uz) * wf_x1;
  ut_x2 = (blData2.param.uz - blData1.param.uz) * wf_x2;
  ut_t1 = (blData2.param.uz - blData1.param.uz) * wf_t1;
  ut_t2 = (blData2.param.uz - blData1.param.uz) * wf_t2;
  ut_d1 = (blData2.param.uz - blData1.param.uz) * wf_d1;
  ut_d2 = (blData2.param.uz - blData1.param.uz) * wf_d2;
  ut_u1 = (blData2.param.uz - blData1.param.uz) * wf_u1 + (1 - wf);
  ut_u2 = (blData2.param.uz - blData1.param.uz) * wf_u2 + wf;
  ut_ms = (blData2.param.uz - blData1.param.uz) * wf_ms;
  ut_re = (blData2.param.uz - blData1.param.uz) * wf_re;
  ut_xf = (blData2.param.uz - blData1.param.uz) * wf_xf;

  //---- set primary "t" variables at xt  (really placed into "2" variables)
  blData2.param.xz = xt;
  blData2.param.tz = tt;
  blData2.param.dz = dt;
  blData2.param.uz = ut;

  blData2.param.amplz = amcrit;
  blData2.param.sz = 0.0;

  //---- calculate laminar secondary "t" variables
  blkin();
  blvar(blData2, 1);

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
  for (int k = 1; k < 3; k++) {
    blrez[k] = vsrez[k];
    blm[k] = vsm[k] + vs2(k, 1) * tt_ms + vs2(k, 2) * dt_ms +
             vs2(k, 3) * ut_ms + vs2(k, 4) * xt_ms;
    blr[k] = vsr[k] + vs2(k, 1) * tt_re + vs2(k, 2) * dt_re +
             vs2(k, 3) * ut_re + vs2(k, 4) * xt_re;
    blx[k] = vsx[k] + vs2(k, 1) * tt_xf + vs2(k, 2) * dt_xf +
             vs2(k, 3) * ut_xf + vs2(k, 4) * xt_xf;

    bl1(k, 0) = vs1(k, 0) + vs2(k, 1) * tt_a1 + vs2(k, 2) * dt_a1 +
                vs2(k, 3) * ut_a1 + vs2(k, 4) * xt_a1;
    bl1(k, 1) = vs1(k, 1) + vs2(k, 1) * tt_t1 + vs2(k, 2) * dt_t1 +
                vs2(k, 3) * ut_t1 + vs2(k, 4) * xt_t1;
    bl1(k, 2) = vs1(k, 2) + vs2(k, 1) * tt_d1 + vs2(k, 2) * dt_d1 +
                vs2(k, 3) * ut_d1 + vs2(k, 4) * xt_d1;
    bl1(k, 3) = vs1(k, 3) + vs2(k, 1) * tt_u1 + vs2(k, 2) * dt_u1 +
                vs2(k, 3) * ut_u1 + vs2(k, 4) * xt_u1;
    bl1(k, 4) = vs1(k, 4) + vs2(k, 1) * tt_x1 + vs2(k, 2) * dt_x1 +
                vs2(k, 3) * ut_x1 + vs2(k, 4) * xt_x1;

    bl2(k, 0) = 0.0;
    bl2(k, 1) = vs2(k, 1) * tt_t2 + vs2(k, 2) * dt_t2 + vs2(k, 3) * ut_t2 +
                vs2(k, 4) * xt_t2;
    bl2(k, 2) = vs2(k, 1) * tt_d2 + vs2(k, 2) * dt_d2 + vs2(k, 3) * ut_d2 +
                vs2(k, 4) * xt_d2;
    bl2(k, 3) = vs2(k, 1) * tt_u2 + vs2(k, 2) * dt_u2 + vs2(k, 3) * ut_u2 +
                vs2(k, 4) * xt_u2;
    bl2(k, 4) = vs2(k, 1) * tt_x2 + vs2(k, 2) * dt_x2 + vs2(k, 3) * ut_x2 +
                vs2(k, 4) * xt_x2;
  }

  //**** second, set up turbulent part between xt and x2  ****

  //---- calculate equilibrium shear coefficient cqt at transition point
  blvar(blData2, 2);

  //---- set initial shear coefficient value st at transition point
  //-    ( note that cq2, cq2_t2, etc. are really "cqt", "cqt_tt", etc.)

  ctr = 1.8 * exp(-3.3 / (blData2.hkz.scalar - 1.0));
  ctr_hk2 = ctr * 3.3 / (blData2.hkz.scalar - 1.0) / (blData2.hkz.scalar - 1.0);

  st = ctr * blData2.cqz.scalar;
  st_tt = ctr * blData2.cqz.t() + blData2.cqz.scalar * ctr_hk2 * blData2.hkz.t();
  st_dt = ctr * blData2.cqz.d() + blData2.cqz.scalar * ctr_hk2 * blData2.hkz.d();
  st_ut = ctr * blData2.cqz.u() + blData2.cqz.scalar * ctr_hk2 * blData2.hkz.u();
  st_ms = ctr * blData2.cqz.ms() + blData2.cqz.scalar * ctr_hk2 * blData2.hkz.ms();
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
  blvar(blData2, 2);

  stepbl();
  restoreblData(2);

  //---- calculate xt-x2 midpoint cfm value
  blmid(2);

  //---- set up newton system for dct, dth, dds, due, dxi  at  xt and x2
  bldif(2);

  //---- convert sensitivities wrt "t" variables into sensitivities
  //-    wrt "1" and "2" variables as done before for the laminar part
  Matrix<double, 5, 5> bt1_right = Matrix<double, 5, 5> {
    {st_a1, st_t1, st_d1, st_u1, st_x1},
    {tt_a1, tt_t1, tt_d1, tt_u1, tt_x1},
    {dt_a1, dt_t1, dt_d1, dt_u1, dt_x1},
    {ut_a1, ut_t1, ut_d1, ut_u1, ut_x1},
    {xt_a1, xt_t1, xt_d1, xt_u1, xt_x1}
  };

    Matrix<double, 5, 5> bt2_right = Matrix<double, 5, 5> {
    {0, st_t2, st_d2, st_u2, st_x2},
    {0, tt_t2, tt_d2, tt_u2, tt_x2},
    {0, dt_t2, dt_d2, dt_u2, dt_x2},
    {0, ut_t2, ut_d2, ut_u2, ut_x2},
    {0, xt_t2, xt_d2, xt_u2, xt_x2}
  };
  bt1.block(0, 0, 3, 5) = vs1.block(0, 0, 3, 5) * bt1_right;
  bt2.block(0, 0, 3, 5) = vs1.block(0, 0, 3, 5) * bt2_right;
  bt2 += vs2;
  for (int k = 0; k < 3; k++) {
    btrez[k] = vsrez[k];
    btm[k] = vsm[k] + vs1(k, 0) * st_ms + vs1(k, 1) * tt_ms +
             vs1(k, 2) * dt_ms + vs1(k, 3) * ut_ms + vs1(k, 4) * xt_ms;
    btr[k] = vsr[k] + vs1(k, 0) * st_re + vs1(k, 1) * tt_re +
             vs1(k, 2) * dt_re + vs1(k, 3) * ut_re + vs1(k, 4) * xt_re;
    btx[k] = vsx[k] + vs1(k, 0) * st_xf + vs1(k, 1) * tt_xf +
             vs1(k, 2) * dt_xf + vs1(k, 3) * ut_xf + vs1(k, 4) * xt_xf;

  }

  //---- add up laminar and turbulent parts to get final system
  //-    in terms of honest-to-god "1" and "2" variables.
  vsrez[0] = btrez[0];
  vsrez[1] = blrez[1] + btrez[1];
  vsrez[2] = blrez[2] + btrez[2];
  vsm[0] = btm[0];
  vsm[1] = blm[1] + btm[1];
  vsm[2] = blm[2] + btm[2];
  vsr[0] = btr[0];
  vsr[1] = blr[1] + btr[1];
  vsr[2] = blr[2] + btr[2];
  vsx[0] = btx[0];
  vsx[1] = blx[1] + btx[1];
  vsx[2] = blx[2] + btx[2];
  vs1.row(0) = bt1.row(0);
  vs2.row(0) = bt2.row(0);
  vs1.middleRows(1, 2) = bl1.middleRows(1, 2) + bt1.middleRows(1, 2);
  vs2.middleRows(1, 2) = bl2.middleRows(1, 2) + bt2.middleRows(1, 2);

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
    for (int ibl = 2; ibl <= nbl.get(is); ibl++) {
      double dui = 0.0;
      for (int js = 1; js <= 2; js++) {
        for (int jbl = 2; jbl <= nbl.get(js); jbl++) {
          double ue_m = -vti.get(is)[ibl] * vti.get(js)[jbl] * dij[ipan.get(is)[ibl]][ipan.get(js)[jbl]];
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
    uinv.get(is)[1] = 0.0;
    uinv_a.get(is)[1] = 0.0;
    for (int ibl = 2; ibl <= nbl.get(is); ibl++) {
      int i = ipan.get(is)[ibl];
      uinv.get(is)[ibl] = vti.get(is)[ibl] * qinv[i];
      uinv_a.get(is)[ibl] = vti.get(is)[ibl] * qinv_a[i];
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

  //int i = 0, is = 0, iv, iw, j, js, jv, ibl, jbl, kbl = 0;
  
  SidePair<VectorXd> unew, u_ac;
  unew.top = VectorXd::Zero(IVX);
  unew.bottom = VectorXd::Zero(IVX);
  u_ac.top = VectorXd::Zero(IVX);
  u_ac.bottom = VectorXd::Zero(IVX);

  double qnew[IQX], q_ac[IQX];
  memset(qnew, 0, IQX * sizeof(double));
  memset(q_ac, 0, IQX * sizeof(double));
  double dalmax = 0.0, dalmin = 0.0, dclmax = 0.0, dclmin = 0.0;
  double dac = 0.0, dhi = 0.0, dlo = 0.0, dctau, dthet, dmass, duedg, ddstr;
  double dn1, dn2, dn3, dn4, rdn1, rdn2, rdn3, rdn4;
  double dswaki, hklim, msq, dsw;
  double dui, dui_ac, ue_m , uinv_ac,
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
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 2; ibl <= nbl.get(is); ibl++) {
      int i = ipan.get(is)[ibl];
      dui = 0.0;
      dui_ac = 0.0;
      for (int js = 1; js <= 2; js++) {
        for (int jbl = 2; jbl <= nbl.get(js); jbl++) {
          int j = ipan.get(js)[jbl];
          int jv = isys.get(js)[jbl];
          ue_m = -vti.get(is)[ibl] * vti.get(js)[jbl] * dij[i][j];
          dui = dui + ue_m * (mass.get(js)[jbl] + vdel[2][0][jv]);
          dui_ac = dui_ac + ue_m * (-vdel[2][1][jv]);
        }
      }

      //------- uinv depends on "ac" only if "ac" is alpha
      if (lalfa)
        uinv_ac = 0.0;
      else
        uinv_ac = uinv_a.get(is)[ibl];

      unew.get(is)[ibl] = uinv.get(is)[ibl] + dui;
      u_ac.get(is)[ibl] = uinv_ac + dui_ac;
    }
  }

  //--- set new qtan from new ue with appropriate sign change

  for (int is = 1; is <= 2; is++) {
    for (int ibl = 2; ibl <= iblte.get(is); ibl++) {
      int i = ipan.get(is)[ibl];
      qnew[i] = vti.get(is)[ibl] * unew.get(is)[ibl];
      q_ac[i] = vti.get(is)[ibl] * u_ac.get(is)[ibl];
    }
  }

  //--- calculate new cl from this new qtan
  beta = sqrt(1.0 - minf * minf);
  beta_msq = -0.5 / beta;

  bfac = 0.5 * minf * minf / (1.0 + beta);
  bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;

  clnew = 0.0;
  cl_a = 0.0;
  cl_ms = 0.0;
  cl_ac = 0.0;

  cginc = 1.0 - (qnew[1] / qinf) * (qnew[1] / qinf);
  cpg1 = cginc / (beta + bfac * cginc);
  cpg1_ms = -cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

  cpi_q = -2.0 * qnew[1] / qinf / qinf;
  cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
  cpg1_ac = cpc_cpi * cpi_q * q_ac[1];

  for (int i = 1; i <= n; i++) {
    int ip = i + 1;
    if (i == n) ip = 1;

    cginc = 1.0 - (qnew[ip] / qinf) * (qnew[ip] / qinf);
    const double cpg2 = cginc / (beta + bfac * cginc);
    const double cpg2_ms = -cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

    cpi_q = -2.0 * qnew[ip] / qinf / qinf;
    cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
    const double cpg2_ac = cpc_cpi * cpi_q * q_ac[ip];

    Matrix2d rotateMatrix = Matrix2d {
      {cos(alfa), sin(alfa)},
      {-sin(alfa), cos(alfa)}
    };

    Vector2d dpoint = rotateMatrix * (points.col(ip) - points.col(i));

    const double ag = 0.5 * (cpg2 + cpg1);
    const double ag_ms = 0.5 * (cpg2_ms + cpg1_ms);
    const double ag_ac = 0.5 * (cpg2_ac + cpg1_ac);

    clnew = clnew + dpoint.x() * ag;
    cl_a = cl_a + dpoint.y() * ag;
    cl_ms = cl_ms + dpoint.x() * ag_ms;
    cl_ac = cl_ac + dpoint.x() * ag_ac;

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

  for (int is = 1; is <= 2; is++) {
    for (int ibl = 2; ibl <= nbl.get(is); ibl++) {
      int iv = isys.get(is)[ibl];
      //------- set changes without underrelaxation
      dctau = vdel[0][0][iv] - dac * vdel[0][1][iv];
      dthet = vdel[1][0][iv] - dac * vdel[1][1][iv];
      dmass = vdel[2][0][iv] - dac * vdel[2][1][iv];
      duedg = unew.get(is)[ibl] + dac * u_ac.get(is)[ibl] - uedg.get(is)[ibl];
      ddstr = (dmass - dstr.get(is)[ibl] * duedg) / uedg.get(is)[ibl];
      //------- normalize changes
      if (ibl < itran.get(is))
        dn1 = dctau / 10.0;
      else
        dn1 = dctau / ctau.get(is)[ibl];
      dn2 = dthet / thet.get(is)[ibl];
      dn3 = ddstr / dstr.get(is)[ibl];
      dn4 = fabs(duedg) / 0.25;
      //------- accumulate for rms change
      rmsbl = rmsbl + dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;
      //------- see if ctau needs underrelaxation
      rdn1 = rlx * dn1;
      if (fabs(dn1) > fabs(rmxbl)) {
        rmxbl = dn1;
        if (ibl < itran.get(is)) vmxbl = "n";
        if (ibl >= itran.get(is)) vmxbl = "c";
      }
      if (rdn1 > dhi) rlx = dhi / dn1;
      if (rdn1 < dlo) rlx = dlo / dn1;
      //------- see if theta needs underrelaxation
      rdn2 = rlx * dn2;
      if (fabs(dn2) > fabs(rmxbl)) {
        rmxbl = dn2;
        vmxbl = "t";
      }
      if (rdn2 > dhi) rlx = dhi / dn2;
      if (rdn2 < dlo) rlx = dlo / dn2;
      //------- see if dstar needs underrelaxation
      rdn3 = rlx * dn3;
      if (fabs(dn3) > fabs(rmxbl)) {
        rmxbl = dn3;
        vmxbl = "d";
      }
      if (rdn3 > dhi) rlx = dhi / dn3;
      if (rdn3 < dlo) rlx = dlo / dn3;

      //------- see if ue needs underrelaxation
      rdn4 = rlx * dn4;
      if (fabs(dn4) > fabs(rmxbl)) {
        rmxbl = duedg;
        vmxbl = "u";
      }
      if (rdn4 > dhi) rlx = dhi / dn4;
      if (rdn4 < dlo) rlx = dlo / dn4;
    }
  }

  //--- set true rms change
  rmsbl = sqrt(rmsbl / (4.0 * double(nbl.top + nbl.bottom)));

  if (lalfa) {
    //---- set underrelaxed change in reynolds number from change in lift
    cl = cl + rlx * dac;
  } else {
    //---- set underrelaxed change in alpha
    alfa = alfa + rlx * dac;
  }

  //--- update bl variables with underrelaxed changes
  for (int is = 1; is <= 2; is++) {
    for (int ibl = 2; ibl <= nbl.get(is); ibl++) {
      int iv = isys.get(is)[ibl];

      dctau = vdel[0][0][iv] - dac * vdel[0][1][iv];
      dthet = vdel[1][0][iv] - dac * vdel[1][1][iv];
      dmass = vdel[2][0][iv] - dac * vdel[2][1][iv];
      duedg = unew.get(is)[ibl] + dac * u_ac.get(is)[ibl] - uedg.get(is)[ibl];
      ddstr = (dmass - dstr.get(is)[ibl] * duedg) / uedg.get(is)[ibl];

      ctau.get(is)[ibl] = ctau.get(is)[ibl] + rlx * dctau;
      thet.get(is)[ibl] = thet.get(is)[ibl] + rlx * dthet;
      dstr.get(is)[ibl] = dstr.get(is)[ibl] + rlx * ddstr;
      uedg.get(is)[ibl] = uedg.get(is)[ibl] + rlx * duedg;

      if (ibl > iblte.get(is)) {
        dswaki = wgap[ibl - iblte.get(is)];
      } else
        dswaki = 0.0;
      //------- eliminate absurd transients
      if (ibl >= itran.get(is)) ctau.get(is)[ibl] = std::min(ctau.get(is)[ibl], 0.25);

      if (ibl <= iblte.get(is))
        hklim = 1.02;
      else
        hklim = 1.00005;

      msq = uedg.get(is)[ibl] * uedg.get(is)[ibl] * hstinv /
            (gamm1 * (1.0 - 0.5 * uedg.get(is)[ibl] * uedg.get(is)[ibl] * hstinv));
      dsw = dstr.get(is)[ibl] - dswaki;
      dslim(dsw, thet.get(is)[ibl], msq, hklim);
      dstr.get(is)[ibl] = dsw + dswaki;

      //------- set new mass defect (nonlinear update)
      mass.get(is)[ibl] = dstr.get(is)[ibl] * uedg.get(is)[ibl];
    }
  }

  //--- equate upper wake arrays to lower wake arrays
  for (int kbl = 1; kbl <= nbl.bottom - iblte.bottom; kbl++) {
    ctau.top[iblte.top + kbl] = ctau.bottom[iblte.bottom + kbl];
    thet.top[iblte.top + kbl] = thet.bottom[iblte.bottom + kbl];
    dstr.top[iblte.top + kbl] = dstr.bottom[iblte.bottom + kbl];
    uedg.top[iblte.top + kbl] = uedg.bottom[iblte.bottom + kbl];
    ctq.top[iblte.top + kbl] = ctq.bottom[iblte.bottom + kbl];
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
    for (int ibl = 1; ibl <= nbl.top; ibl++) {
      uedg.top[ibl] = uinv.top[ibl];
    }
    for (int ibl = 1; ibl <= nbl.bottom; ibl++) {
      uedg.bottom[ibl] = uinv.bottom[ibl];
    }
  }

  if (lvconv) {
    //	----- set correct cl if converged point exists
    qvfue();
    
    if (lvisc) {
      cpv = cpcalc(n + nw, qvis, qinf, minf);
      cpi = cpcalc(n + nw, qinv, qinf, minf);
    } else
      cpi = cpcalc(n, qinv, qinf, minf);

    gamqv();
    clcalc(cmref);
    cdcalc();
  }

  //	---- set up source influence matrix if it doesn't exist
  if (!lwdij || !ladij) qdcalc();

  return true;
}

bool XFoil::ViscalEnd() {

  cpi = cpcalc(n + nw, qinv, qinf, minf);
  cpv = cpcalc(n + nw, qvis, qinf, minf);

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
    minf_cl = getActualMach(cl, mach_type);
    reinf_cl = getActualReynolds(cl, reynolds_type);
    comset();
  } else {  //	------- set new inviscid speeds qinv and uinv for new alpha
    qiset();
    uicalc();
  }

  qvfue();   //	------ calculate edge velocities qvis(.) from uedg(..)
  gamqv();   //	------ set gam distribution from qvis
  stmove();  //	------ relocate stagnation point

  //	------ set updated cl,cd
  clcalc(cmref);
  cdcalc();

  if (rmsbl < eps1) {
    lvconv = true;
    avisc = alfa;
    mvisc = minf;
    writeString("----------CONVERGED----------\n\n");
  }

  return true;
}

bool XFoil::xicalc() {
  //-------------------------------------------------------------
  //     sets bl arc length array on each airfoil side and wake
  //-------------------------------------------------------------

  xssi.top[1] = 0.0;
  for (int ibl = 2; ibl <= iblte.top; ibl++) {
    xssi.top[ibl] = sst - spline_length[ipan.top[ibl]];
  }

  xssi.bottom[1] = 0.0;
  for (int ibl = 2; ibl <= iblte.bottom; ibl++) {
    xssi.bottom[ibl] = spline_length[ipan.bottom[ibl]] - sst;
  }

  xssi.bottom[iblte.bottom + 1] = xssi.bottom[iblte.bottom];
  for (int ibl = iblte.bottom + 2; ibl <= nbl.bottom; ibl++) {
    xssi.bottom[ibl] =
        xssi.bottom[ibl - 1] + (points.col(ipan.bottom[ibl]) - points.col(ipan.bottom[ibl] - 1)).norm();
  }

  //---- trailing edge flap length to te gap ratio
  const double telrat = 2.50;

  //---- set up parameters for te flap cubics

  const double crosp = cross2(dpoints_ds.col(n).normalized(), dpoints_ds.col(1).normalized());
  double dwdxte = crosp / sqrt(1.0 - crosp * crosp);

  //---- limit cubic to avoid absurd te gap widths
  dwdxte = std::max(dwdxte, -3.0 / telrat);
  dwdxte = std::min(dwdxte, 3.0 / telrat);

  const double aa = 3.0 + telrat * dwdxte;
  const double bb = -2.0 - telrat * dwdxte;

  if (sharp) {
    for (int iw = 1; iw <= nw; iw++) wgap[iw] = 0.0;
  }

  else {
    //----- set te flap (wake gap) array
    for (int iw = 1; iw <= nw; iw++) {
      const double zn = 1.0 - (xssi.bottom[iblte.bottom + iw] - xssi.bottom[iblte.bottom]) / (telrat * ante);
      wgap[iw] = 0.0;
      if (zn >= 0.0) wgap[iw] = ante * (aa + bb * zn) * zn * zn;
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

  if (xstrip.get(is) >= 1.0) {
    return xssi.get(is)[iblte.get(is)];
  }

  Vector2d point_chord = point_te - point_le;

  //---- calculate chord-based x/c, y/c
  for (int i = 1; i <= n; i++) {
    w1[i] = (points.col(i) - point_le).dot(point_chord.normalized());
    w2[i] = cross2(points.col(i) - point_le, point_chord.normalized());
  }

  w3 = spline::splind(w1, spline_length, n);
  w4 = spline::splind(w2, spline_length, n);

  if (is == 1) {
    //----- set approximate arc length of forced transition point for sinvrt
    str = sle + (spline_length[1] - sle) * xstrip.top;

    //----- calculate actual arc length
    str = spline::sinvrt(str, xstrip.top, w1, w3, spline_length, n);

    //----- set bl coordinate value
    xiforc = std::min((sst - str), xssi.get(is)[iblte.get(is)]);
  } else {
    //----- same for bottom side

    str = sle + (spline_length[n] - sle) * xstrip.bottom;
    str = spline::sinvrt(str, xstrip.bottom, w1, w3, spline_length, n);
    xiforc = std::min((str - sst), xssi.get(is)[iblte.get(is)]);
  }

  if (xiforc < 0.0) {
    ss << " ***  stagnation point is past trip on side " << is << "\n";
    writeString(ss.str());

    xiforc = xssi.get(is)[iblte.get(is)];
  }

  return xiforc;
}

bool XFoil::xyWake() {
  //-----------------------------------------------------
  //     sets wake coordinate array for current surface
  //     vorticity and/or mass source distributions.
  //-----------------------------------------------------
  double ds1, sx, sy, smod;
  //
  writeString("   Calculating wake trajectory ...\n");
  //
  //--- number of wake points
  nw = n / 8 + 2;
  if (nw > IWX) {
    writeString(
        " XYWake: array size (IWX) too small.\n  Last wake point index reduced.");
    nw = IWX;
  }

  ds1 = 0.5 * (spline_length[2] - spline_length[1] + spline_length[n] - spline_length[n - 1]);
  setexp(snew.data() + n, ds1, waklen * chord, nw);

  point_te = 0.5 * (points.col(1) + points.col(n));

  //-- set first wake point a tiny distance behind te
  sx = 0.5 * (dpoints_ds.col(n).y() - dpoints_ds.col(1).y());
  sy = 0.5 * (dpoints_ds.col(1).x() - dpoints_ds.col(n).x());
  smod = sqrt(sx * sx + sy * sy);
  normal_vectors.col(n + 1).x() = sx / smod;
  normal_vectors.col(n + 1).y() = sy / smod;
  points.col(n + 1).x() = point_te.x() - 0.0001 * normal_vectors.col(n + 1).y();
  points.col(n + 1).y() = point_te.y() + 0.0001 * normal_vectors.col(n + 1).x();
  spline_length[n + 1] = spline_length[n];

  //---- calculate streamfunction gradient components at first point
  Vector2d psi = {
    psilin(n + 1, points.col(n + 1), {1.0, 0.0}, false).psi_ni,
    psilin(n + 1, points.col(n + 1), {0.0, 1.0}, false).psi_ni
  };

  //---- set unit vector normal to wake at first point
  normal_vectors.col(n + 2) = -psi.normalized();

  //---- set angle of wake panel normal
  apanel[n + 1] = atan2(psi.y(), psi.x());

  //---- set rest of wake points
  for (int i = n + 2; i <= n + nw; i++) {
    const double ds = snew[i] - snew[i - 1];

    //------ set new point ds downstream of last point
    points.col(i).x() = points.col(i - 1).x() - ds * normal_vectors.col(i).y();
    points.col(i).y() = points.col(i - 1).y() + ds * normal_vectors.col(i).x();
    spline_length[i] = spline_length[i - 1] + ds;

    if (i != n + nw) {
      //---- calculate streamfunction gradient components at first point
      Vector2d psi = {
        psilin(i, points.col(i), {1.0, 0.0}, false).psi_ni,
        psilin(i, points.col(i), {0.0, 1.0}, false).psi_ni
      };

      //---- set unit vector normal to wake at first point
      normal_vectors.col(i + 1) = -psi.normalized();

      //------- set angle of wake panel normal
      apanel[i] = atan2(psi.y(), psi.x());
    }
  }

  //---- set wake presence flag and corresponding alpha
  lwake = true;
  awake = alfa;

  //---- old source influence matrix is invalid for the new wake geometry
  lwdij = false;

  return true;
}

bool XFoil::isValidFoilAngles(Matrix2Xd points) {
  
  double max_angle = cang(points.middleCols(INDEX_START_WITH, points.cols() - INDEX_START_WITH));
  return max_angle <= angtol;
}

bool XFoil::isValidFoilPointSize(Matrix2Xd points) {
  return points.cols() >= 3 + INDEX_START_WITH;
}
