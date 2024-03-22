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
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*****************************************************************************/

/**
 *@file This class defines the Xfoil object.
 */

#pragma once

/**
*@class XFoil
*@brief  The class which defines the XFoil object.

This is a translation to C++ of the original Fortran code of Mark Drela and
Harold Youngren. See http://raphael.mit.edu/xfoil for more information.
*/

#include <complex>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <iterator>
#include <numeric>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/StdVector"

#include "model/matrix.hpp"
#include "model/spline.hpp"
#include "xfoil_params.h"

using namespace std;
using namespace Eigen;
//------ derived dimensioning limit parameters

struct blData {
 public:
  double xz, uz, tz, dz, sz, amplz, uz_uei, uz_ms, dwz, 
      hz, hz_tz, hz_dz, 
      mz, mz_uz, mz_ms, 
      rz, rz_uz, rz_ms, 
      vz, vz_uz, vz_ms, vz_re,
      hkz, hkz_uz, hkz_tz, hkz_dz, hkz_ms, 
      hsz, hsz_uz, hsz_tz, hsz_dz, hsz_ms, hsz_re, 
      hcz, hcz_uz, hcz_tz, hcz_dz, hcz_ms, 
      rtz, rtz_uz, rtz_tz, rtz_ms, rtz_re, 
      cfz, cfz_uz, cfz_tz, cfz_dz, cfz_ms, cfz_re,
      diz, diz_uz, diz_tz, diz_dz, diz_sz, diz_ms, diz_re,
      usz, usz_uz, usz_tz, usz_dz, usz_ms, usz_re, 
      cqz, cqz_uz, cqz_tz, cqz_dz, cqz_ms, cqz_re, 
      dez, dez_uz, dez_tz, dez_dz, dez_ms;
};

class XFoil {
 public:
  XFoil();
  virtual ~XFoil();

 public:

  bool isValidFoilAngles(Matrix2Xd points);
  bool isValidFoilPointSize(Matrix2Xd points);
  
  bool initialize();
  bool initXFoilGeometry(int fn, const double *fx, const double *fy, double *fnx,
                         double *fny);

  enum class ReynoldsType {
    CONSTANT,
    FIXED_LIFT,
    FIXED_LIFT_AND_DYNAMIC_PRESSURE
  };
  enum class MachType {
    CONSTANT,
    FIXED_LIFT,
    FIXED_LIFT_AND_DYNAMIC_PRESSURE
  };
  ReynoldsType reynolds_type;
  MachType mach_type;
  bool initXFoilAnalysis(double Re, double alpha, double Mach, double NCrit,
                         double XtrTop, double XtrBot, ReynoldsType reType, MachType maType,
                         bool bViscous, std::stringstream &outStream);


  bool clcalc(double xref, double yref);

  void writeString(std::string str, bool bFullReport = false);
  bool specal();
  bool speccl();
  bool viscal();
  bool ViscalEnd();
  bool ViscousIter();
  
  bool abcopy(Matrix2Xd copyFrom);

  bool isBLInitialized() const { return lblini; }
  void setBLInitialized(bool bInitialized) { lblini = bInitialized; }

  double QInf() const { return qinf; }
  void setQInf(double v) { qinf = v; }

  double alpha() const { return alfa; }
  void setAlpha(double aoa) { alfa = aoa; }

  double ClSpec() const { return clspec; }
  void setClSpec(double cl) { clspec = cl; }

  static bool isCancelled() { return s_bCancel; }
  static void setCancel(bool bCancel) { s_bCancel = bCancel; }
  static void setFullReport(bool bFull) { s_bFullReport = bFull; }
  static bool fullReport() { return s_bFullReport; }
  static double VAccel() { return vaccel; }
  static void setVAccel(double accel) { vaccel = accel; }

  bool setMach();

  bool comset();
  bool mrcl(double cls, double &m_cls, double &r_cls);
  bool restoreblData(int icom);
  bool saveblData(int icom);

  bool apcalc();
  class AxResult {
    public:
    double ax;
    double ax_hk1;
    double ax_hk2;
    double ax_t1;
    double ax_t2;
    double ax_rt1;
    double ax_rt2;
    double ax_a1;
    double ax_a2;
  };
  AxResult axset(double hk1, double thet1, double rt1, double a1, double hk2,
    double thet2, double rt2, double a2, double acrit);
  bool bldif(int ityp);
  bool blkin();
  bool blmid(int ityp);
  bool blprv(double xsi, double ami, double cti, double thi, double dsi,
             double dswaki, double uei);
  bool blsolve();
  bool blsys();
  bool blvar(int ityp);
  double cang(Matrix2Xd points);
  bool cdcalc();

  /**
   * @brief 
   * 
   */
  class C_f {
    public:
    double cf;
    double hk;
    double rt;
    /** squared freestream mach number at current cl*/
    double msq;
  };
  C_f cfl(double hk, double rt);
  C_f cft(double hk, double rt, double msq);
  enum class SideType {
    TOP = 1,
    BOTTOM = 2
  };
  template <class T>
  class SidePair {
    public:
    T top;
    T bottom;
    //FIXME deprecated
    T get(int side) {
      if (side == 1) {
        return top;
      }
      else if (side == 2) {
        return bottom;
      }
      throw invalid_argument("invalid side type");
    }
  };
  bool cpcalc(int n, const double q[], double qinf, double minf, double cp[]);
  bool dampl(double hk, double th, double rt, double &ax, double &ax_hk,
             double &ax_th, double &ax_rt);
  bool dil(double hk, double rt, double &di, double &di_hk, double &di_rt);
  bool dilw(double hk, double rt, double &di, double &di_hk, double &di_rt);
  bool dslim(double &dstr, double thet, double msq, double hklim);

  bool gamqv();
  bool Gauss(int nn, double z[][6], double r[5]);
  bool getxyf(Matrix2Xd points, Matrix2Xd dpoints_ds, VectorXd s,
              int n, double &tops, double &bots, double xf, double &yf);
  bool ggcalc();
  bool hct(double hk, double msq, double &hc, double &hc_hk, double &hc_msq);
  bool hkin(double h, double msq, double &hk, double &hk_h, double &hk_msq);
  bool hsl(double hk, double &hs, double &hs_hk, double &hs_rt, double &hs_msq);
  bool hst(double hk, double rt, double msq, double &hs, double &hs_hk,
           double &hs_rt, double &hs_msq);
  bool iblpan();
  bool iblsys();
  bool lefind(double &sle, Matrix2Xd points, Matrix2Xd dpoints_ds, VectorXd s, int n);
  
  bool mrchdu();
  bool mrchue();
  Matrix2Xd ncalc(Matrix2Xd point, VectorXd spline_length, int n);
  class PsiResult {
    public:
    double psi;
    double psi_ni;
    static PsiResult sum(PsiResult a, PsiResult b) {
      PsiResult result;
      result.psi = a.psi + b.psi;
      result.psi_ni = a.psi_ni + b.psi_ni;
      return result;
    }
  };
  PsiResult psilin(int i, Vector2d point, Vector2d normal_vector, bool siglin);
  PsiResult psisig(int iNode, int jNode, Vector2d point, Vector2d normal_vector);
  PsiResult psi_te(int i, Vector2d point, Vector2d normal_vector);
  
  bool pswlin(int i, double xi, double yi, double nxi, double nyi, double &psi,
              double &psi_ni);
  bool qdcalc();
  bool qiset();
  bool qvfue();
  bool qwcalc();
  
  bool setbl();
  bool setexp(double s[], double ds1, double smax, int nn);

  bool stepbl();
  bool stfind();
  bool stmove();
  bool tecalc();
  bool tesys(double cte, double tte, double dte);
  bool trchek();
  bool trdif();
  bool ueset();
  bool uicalc();
  bool update();
  bool xicalc();
  double xifset(int is);
  bool xyWake();
  double aint(double number);
  double atanc(double y, double x, double thold);

  double sign(double a, double b);

 public:
  static double vaccel;
  static bool s_bCancel;
  static bool s_bFullReport;
  std::stringstream *m_pOutStream;

  double clspec;

  int nc1;

  Matrix2Xd normal_vectors;

  double cl, cm, cd, cpi[IZX], cpv[IZX], acrit;
  double xcp;
  double alfa, avisc, awake, reinf1, qinf, mvisc, rmsbl, ante;
  double cpmn;
  double minf, reinf;
  bool lalfa, lvisc, lvconv, lwake;
  double qgamm[IBX + 1];
  double hmom;

  double rmxbl;

  double qtan1, qtan2;
  double amax;  // needed for preprocessing
  double minf1;
  bool lblini, lipan;
  
  int n, ipan[IVX][ISX];
  SidePair<int> iblte, nbl;
  
  Matrix2Xd points; //formerly x,y
  SidePair<double> xstrip;
  
  double qvis[IZX];
  
  double adeg, xcmref, ycmref;
  double tklam;
  Matrix2Xd dpoints_ds; //formerly xp, yp
  VectorXd spline_length;
  double dtor;

  double thet[IVX][ISX], ctau[IVX][ISX], ctq[IVX][ISX], uedg[IVX][ISX];

  double dstr[IVX][ISX];
  
  int itran[ISX];

 public: //private:
  double qf0[IQX + 1], qf1[IQX + 1], qf2[IQX + 1], qf3[IQX + 1];

  double rlx;

  double minf_cl, reinf_cl;
  const double angtol = 40.0; // foil angle tolerance

  blData blsav[3];

  bool lgamu, sharp, lqaij, ladij, lwdij;

  double sccon, gacon, gbcon, gbc0, gbc1, gccon, dlcon, ctcon;

  int nsys;
  double isys[IVX][ISX];
  VectorXd snew;

  double sle;
  Vector2d point_le;
  Vector2d point_te;
  double chord, wgap[IWX], waklen;
  int nw, ist;


  double cl_alf, cl_msq;
  double psio, cosa, sina, gamma, gamm1;
  double tkl_msq, cpstar, qstar;
  double xssi[IVX][ISX], uinv[IVX][ISX], mass[IVX][ISX];

  double vti[IVX][ISX];
  double uinv_a[IVX][ISX];
  double gam[IQX], gam_a[IQX], sig[IZX];
  Matrix2Xd gamu;
  double apanel[IZX], sst, sst_go, sst_gp, gamte, sigte;
  double dste, aste;
  double qinv[IZX], qinvu[IZX][3], qinv_a[IZX];
  double dzdg[IQX], dzdm[IZX], dqdg[IQX], dqdm[IZX];
  FullPivLU<MatrixXd> psi_gamma_lu;
  
  double bij[IQX][IZX], dij[IZX][IZX];
  double cij[IWX][IQX];

  double vs1[5][6], vs2[5][6], vsrez[5], vsr[5], vsm[5], vsx[5];

  bool trforc, simi, tran, turb, wake, trfree;

  double qinfbl, tkbl, tkbl_ms, rstbl, rstbl_ms, hstinv, hstinv_ms;
  double reybl, reybl_ms, reybl_re, gm1bl, hvrat, bule, xiforc, amcrit;

  blData blData1, blData2;

  int imxbl, ismxbl;

  double cfm, cfm_ms, cfm_re, cfm_u1, cfm_t1, cfm_d1, cfm_u2, cfm_t2, cfm_d2;
  double xt, xt_a1, xt_ms, xt_re, xt_xf, xt_x1, xt_t1, xt_d1, xt_u1, xt_x2,
      xt_t2, xt_d2, xt_u2;
  double va[4][3][IZX], vb[4][3][IZX], vdel[4][3][IZX], vm[4][IZX][IZX], vz[4][3];


  /*
  c
  c-    sccon  =  shear coefficient lag constant
  c-    gacon  =  g-beta locus constants...
  c-    gbcon  =  g = gacon * sqrt(1.0 + gbcon*beta)
  c-    gccon  =         + gccon / [h*rtheta*sqrt(cf/2)]   <-- wall term
  c-    dlcon  =  wall/wake dissipation length ratio  lo/l
  c-    ctcon  =  ctau weighting coefficient (implied by g-beta constants)

  c   version     version number of this xfoil implementation
  c
  c   fname       airfoil data filename
  c   pfname[.]   polar append filename
  c   pfnamx[.]   polar append x/c dump filename
  c   oname       default overlay airfoil filename
  c   prefix      default filename prefix
  c   name        airfoil name
  c
  c   ispars      ises domain parameters  [not used in xfoil]
  c
  c   q[..]       generic coefficient matrix
  c   dq[.]       generic matrix righthand side
  c
  c   dzdg[.]     dpsi/dgam
  c   dzdn[.]     dpsi/dn
  c   dzdm[.]     dpsi/dsig
  c
  c   dqdg[.]     dqtan/dgam
  c   dqdm[.]     dqtan/dsig
  c   qtan1       qtan at alpha =  0 deg.
  c   qtan2       qtan at alpha = 90 deg.
  c
  c   z_qinf      dpsi/dqinf
  c   z_alfa      dpsi/dalfa
  c   z_qdof0     dpsi/dqdof0
  c   z_qdof1     dpsi/dqdof1
  c   z_qdof2     dpsi/dqdof2
  c   z_qdof3     dpsi/dqdof3
  c
  c   aij[..]     dpsi/dgam  influence coefficient matrix [factored if lqaij=t]
  c   bij[..]     dgam/dsig  influence coefficient matrix
  c   cij[..]     dqtan/dgam influence coefficient matrix
  c   dij[..]     dqtan/dsig influence coefficient matrix
  c   qinv[.]     tangential velocity due to surface vorticity
  c   qvis[.]     tangential velocity due to surface vorticity & mass sources
  c   qinvu[..]   qinv for alpha = 0, 90 deg.
  c   qinv_a[.]   dqinv/dalpha
  c
  c   x[.],y[.]   airfoil [1<i<n] and wake [n+1<i<n+nw] coordinate arrays
  c   xp[.],yp[.] dx/ds, dy/ds arrays for spline evaluation
  c   s[.]        arc length along airfoil [spline parameter]
  c   sle         value of s at leading edge
  c   xle,yle     leading  edge coordinates
  c   xte,yte     trailing edge coordinates
  c   wgap[.]     thickness of "dead air" region inside wake just behind te
  c   waklen      wake length to chord ratio
  c
  c   gam[.]      surface vortex panel strength array
  c   gamu[.2]    surface vortex panel strength arrays for alpha = 0, 90 deg.
  c   gam_a[.]    dgam/dalfa
  c   sig[.]      surface and wake mass defect array
  c
  c   nx[.],ny[.] normal unit vector components at airfoil and wake coordinates
  c   apanel[.]   surface and wake panel angle array [+ counterclockwise]
  c
  c   sst         s value at stagnation point
  c   sst_go      dsst/dgam[ist]
  c   sst_gp      dsst/dgam[ist+1]
  c
  c   gamte       vortex panel strength across finite-thickness te
  c   sigte       source panel strength across finite-thickness te
  c   gamte_a     dgamte/dalfa
  c   sigte_a     dsigte/dalfa
  c   dste        te panel length
  c   ante,aste   projected te thickness perp.,para. to te bisector
  c   sharp       .true.  if  dste.eq.0.0 ,  .false. otherwise
  c
  c   sspec[.]    normalized arc length around airfoil [qspec coordinate]
  c   xspoc[.]    x/c at sspec points
  c   yspoc[.]    y/c at sspec points
  c   qspec[..]   specified surface velocity for inverse calculations
  c   qspecp[..]  dqspec/dsspec
  c   qgamm[.]    surface velocity for current airfoil geometry
  c   ssple       sspec value at airfoil nose
  c
  c   iq1,iq2     target segment endpoint indices on qspec[s] plot
  c   nsp         number of points in qspec array
  c   nqsp        number qspec arrays
  c   iacqsp      1:  alqsp is prescribed for qspec arrays
  c               2:  clqsp is prescribed for qspec arrays
  c   nc1         number of circle plane points, must be 2**n - 1
  c
  c   nname       number of characters in airfoil name
  c   nprefix     number of characters in default filename prefix
  c
  c   alqsp[.]    alpha,cl,cm corresponding to qspec distributions
  c   clqsp[.]
  c   cmqsp[.]
  c   algam       alpha,cl,cm corresponding to qgamm distribution
  c   clgam
  c   cmgam
  c
  c   qf0[.]      shape function for qspec modification
  c   qf1[.]        "
  c   qf2[.]        "
  c   qf3[.]        "
  c   qdof0       shape function weighting coefficient [inverse dof]
  c   qdof1         "
  c   qdof2         "
  c   qdof3         "
  c   clspec      specified cl
  c   ffilt       circle-plane mapping filter parameter
  c
  c   adeg,alfa   angle of attack in degrees, radians
  c   awake       angle of attack corresponding to wake geometry [radians]
  c   avisc       angle of attack corresponding to bl solution   [radians]
  c   mvisc       mach number corresponding to bl solution
  c   cl,cm       current cl and cm calculated from gam[.] distribution
  c   cd          current cd from bl solution
  c   cl_alf      dcl/dalfa
  c   cl_msq      dcl/d[minf^2]
  c
  c   psio        streamfunction inside airfoil
  c   circ        circulation
  c   cosa,sina   cos[alfa], sin[alfa]
  c   qinf        freestream speed    [defined as 1]
  c   gamma,gamm1 gas constant cp/cv, cp/cv - 1
  c   minf1       freestream mach number at cl=1
  c   minf        freestream mach number at current cl
  c   minf_cl     dminf/dcl
  c   tklam       karman-tsien parameter minf^2 / [1 + sqrt[1-minf^2]]^2
  c   tkl_msq     d[tklam]/d[minf^2]
  c   cpstar      sonic pressure coefficient
  c   qstar       sonic speed
  c
  c   ncpref      number of reference cp vs x/c points
  c   xpref[.]    x/c array corresponding to reference cp data array
  c   cpref[.]    reference cp data array
  c   labref      reference cp data descriptor string
  c
  c   nlref       number of characters in labref string
  c   napol[.]    number of points in each stored polar
  c   npol        number of stored polars
  c   ipact       index of "active" polar being accumulated [0 if none are]
  c   icolp[.]    color for each polar
  c   icolr[.]    color for each reference polar
  c
  c   ndref[..]   number of points in each stored reference polar
  c   npolref     number of stored reference polars
  c
  c   verspol[.]  version number of generating-code for each polar
  c   cpol[...]   cl,cd,and other parameters for each polar
  c   cpolxy[.1.] x,y coordinates of airfoil geometry which generated each polar
  c   cpolxy[.2.]
  c   nxypol[.]   number of x,y points in cpolxy array
  c
  c   pxtr[..]    transition locations for each polar
  c   namepol[.]  airfoil names for each polar
  c   codepol[.]  generating-code names for each polar
  c
  c   nameref[.]  name label of reference polar
  c
  c   pi          3.1415926...
  c   dtor        pi / 180    [degrees to radians conversion factor]
  c
  c   cvpar       curvature attraction parameter for airfoil paneling
  c               0 = uniform panel node spacing around airfoil
  c              ~1 = panel nodes strongly bunched in areas of large curvature
  c   cterat      te panel density / le panel density ratio
  c   ctrrat      local refinement panel density / le panel density ratio
  c   xsref1-2    suction  side local refinement x/c limits
  c   xpref1-2    pressure side local refinement x/c limits
  c
  c   n           number of points on airfoil
  c   nb          number of points in buffer airfoil array
  c   nw          number of points in wake
  c   npan        default/specified number of points on airfoil
  c
  c   ist         stagnation point lies between s[ist], s[ist+1]
  c   itmax       max number of newton iterations
  c   nseqex      max number of unconverged sequence points for early exit
  c
  c   retyp       index giving type of re variation with cl ...
  c            ... 1  re constant
  c            ... 2  re ~ 1/sqrt[cl]    [fixed lift]
  c            ... 3  re ~ 1/cl          [fixed lift and dynamic pressure]
  c
  c   matyp       index giving type of ma variation with cl ...
  c            ... 1  ma constant
  c            ... 2  ma ~ 1/sqrt[cl]    [fixed lift]
  c
  c   aijpiv[.]   pivot index array for lu factoring routine
  c
  c   idev        "device" number for normal screen plotting
  c   idevrp      "device" number for replotting [typically for hardcopy]
  c   ipslu       postscript file specifier
  c   ncolor      number of defined colors in colormap
  c   icols[1]    color indices of top side
  c   icols[2]    color indices of bottom side
  c
  c   nover       number of airfoils overlaid on gdes geometry plot
  c
  c   scrnfr      screen fraction taken up by initial plot window
  c   size        plot width [inches]
  c   plotar      plot aspect ratio
  c   xwind,ywind window size in inches
  c   xpage,ypage plot-page size in inches [for hardcopy]
  c   xmarg,ymarg margin dimensions in inches
  c   pfac        scaling factor for  cp
  c   qfac        scaling factor for  q  [surface speed]
  c   vfac        scaling factor for  cp vectors
  c   ch          character width / plot size  ratio
  c   chg         character width / plot size  ratio for geometry plot
  c   chq         character width / plot size  ratio for qspec[s] plot
  c
  c   xofair      x offset for airfoil in  cp vs x plots
  c   yofair      y offset for airfoil in  cp vs x plots
  c   facair      scale factor for airfoil in  cp vs x plots
  c   xofa        x offset for airfoil in  cp vs x plots in airfoil units
  c   yofa        y offset for airfoil in  cp vs x plots in airfoil units
  c   faca        scale factor for airfoil in  cp vs x plots  in airfoil units
  c   uprwt       u/qinf scale factor for profile plotting
  c   cpmax       max cp  in  cp vs x plots
  c   cpmin       min cp  in  cp vs x plots
  c   cpdel       delta cp  in  cp vs x plots
  c
  c   cpolplf[1,icd]  min cd in cd-cl polar plot
  c   cpolplf[2,icd]  max cd in cd-cl polar plot
  c   cpolplf[3,icd]  delta cd in cd-cl polar plot
  c
  c   xcdwid      width of cd   -cl polar plot
  c   xalwid      width of alpha-cl polar plot
  c   xocwid      width of xtr/c-cl polar plot
  c
  c   ok          user question response
  c   limage      .true. if image airfoil is present
  c   lgamu       .true. if gamu  arrays exist for current airfoil geometry
  c   lqinu       .true. if qinvu arrays exist for current airfoil geometry
  c   lvisc       .true. if viscous option is invoked
  c   lalfa       .true. if alpha is specifed, .false. if cl is specified
  c   lwake       .true. if wake geometry has been calculated
  c   lpacc       .true. if each point calculated is to be saved
  c   lblini      .true. if bl has been initialized
  c   lipan       .true. if bl->panel pointers ipan have been calculated
  c   lqaij       .true. if dpsi/dgam matrix has been computed and factored
  c   ladij       .true. if dq/dsig matrix for the airfoil has been computed
  c   lwdij       .true. if dq/dsig matrix for the wake has been comphd
  c   lqvdes      .true. if viscous ue is to be plotted in qdes routines
  c   lqspec      .true. if qspec has been initialized
  c   lqrefl      .true. if reflected qspec is to be plotted in qdes routines
  c   lvconv      .true. if converged bl solution exists
  c   lcpref      .true. if reference data is to be plotted on cp vs x/c plots
  c   lclock      .true. if source airfoil coordinates are clockwise
  c   lpfile      .true. if polar file is ready to be appended to
  c   lpfilx      .true. if polar dump file is ready to be appended to
  c   lppsho      .true. if cl-cd polar is plotted during point sequence
  c   lbflap      .true. if buffer  airfoil flap parameters are defined
  c   lflap       .true. if current airfoil flap parameters are defined
  c   leiw        .true. if unit circle complex number array is initialized
  c   lscini      .true. if old-airfoil circle-plane arc length s[w] exists
  c   lforef      .true. if cl,cd... data is to be plotted on cp vs x/c plots
  c   lnorm       .true. if input buffer airfoil is to be normalized
  c   lgsame      .true. if current and buffer airfoils are identical
  c
  c   lplcam      .true. if thickness and camber are to be plotted
  c   lqsym       .true. if symmetric qspec will be enforced
  c   lgsym       .true. if symmetric geometry will be enforced
  c   lqgrid      .true. if grid is to overlaid on qspec[s] plot
  c   lggrid      .true. if grid is to overlaid on buffer airfoil geometry plot
  c   lgtick      .true. if node tick marks are to be plotted on buffer airfoil
  c   lqslop      .true. if modified qspec[s] segment is to match slopes
  c   lgslop      .true. if modified geometry segment is to match slopes
  c   lcslop      .true. if modified camber line segment is to match slopes
  c   lqsppl      .true. if current qspec[s] in in plot
  c   lgeopl      .true. if current geometry in in plot
  c   lcpgrd      .true. if grid is to be plotted on cp plots
  c   lblgrd      .true. if grid is to be plotted on bl variable plots
  c   lblsym      .true. if symbols are to be plotted on bl variable plots
  c   lcminp      .true. if min cp is to be written to polar file for cavitation
  c   lhmomp      .true. if hinge moment is to be written to polar file
  c
  c   lpgrid      .true. if polar grid overlay is enabled
  c   lpcdw       .true. if polar cdwave is plotted
  c   lplist      .true. if polar listing lines [at top of plot] are enabled
  c   lplegn      .true. if polar legend is enabled
  c
  c   lplot       .true. if plot page is open
  c   lsym        .true. if symbols are to be plotted in qdes routines
  c   liqset      .true. if inverse target segment is marked off in qdes
  c   lclip       .true. if line-plot clipping is to be performed
  c   lvlab       .true. if label is to be plotted on viscous-variable plots
  c   lcurs       .true. if cursor input is to be used for blowups, etc.
  c   lland       .true. if landscape orientation for postscript is used
  c
  c
  c   xb[.],yb[.] buffer airfoil coordinate arrays
  c   xbp[.]      dxb/dsb
  c   ybp[.]      dyb/dsb
  c   sb[.]       spline parameter for buffer airfoil
  c   snew[.]     new panel endpoint arc length array
  c
  c   xbf,ybf     buffer  airfoil flap hinge coordinates
  c   xof,yof     current airfoil flap hinge coordinates
  c   hmom        moment of flap about hinge point
  c   hfx         x-force of flap on hinge point
  c   hfy         y-force of flap on hinge point
  c
  c~~~~~~~~~~~~~~ properties of current buffer airfoil
  c
  c   xbmin,xbmax  limits of xb array
  c   ybmin,ybmax  limits of yb array
  c   sble         le tangency-point sb location
  c   chordb       chord
  c   areab        area
  c   radble       le radius
  c   angbte       te angle  [rad]
  c
  c   ei11ba       bending inertia about axis 1    x^2 dx dy
  c   ei22ba       bending inertia about axis 2    y^2 dx dy
  c   apx1ba       principal axis 1 angle
  c   apx2ba       principal axis 2 angle
  c
  c   ei11bt       bending inertia about axis 1    x^2 t ds
  c   ei22bt       bending inertia about axis 2    y^2 t ds
  c   apx1bt       principal axis 1 angle
  c   apx2bt       principal axis 2 angle
  c
  c   thickb       max thickness
  c   cambrb       max camber
  c
  c~~~~~~~~~~~~~~
  c
  c   xssi[..]    bl arc length coordinate array on each surface
  c   uedg[..]    bl edge velocity array
  c   uinv[..]    bl edge velocity array without mass defect influence
  c   mass[..]    bl mass defect array  [ = uedg*dstr ]
  c   thet[..]    bl momentum thickness array
  c   dstr[..]    bl displacement thickness array
  c   ctau[..]    sqrt[max shear coefficient] array
  c               [in laminar regions, log of amplification ratio]
  c
  c   ctq[..]     sqrt[equilibrium max shear coefficient] array [  "  ]
  c   vti[..]     +/-1 conversion factor between panel and bl variables
  c   uinv_a[..]  duinv/dalfa array
  c
  c   reinf1      reynolds number  vinf c / ve  for cl=1
  c   reinf       reynolds number for current cl
  c   reinf_cl    dreinf/dcl
  c
  c   acrit       log [critical amplification ratio]
  c   xstrip[.]   transition trip  x/c locations [if xtrip > 0],
  c               transition trip -s/s_side locations [if xtrip < 0],
  c   xoctr[.]    actual transition x/c locations
  c   yoctr[.]    actual transition y/c locations
  c   xssitr[.]   actual transition xi locations
  c
  c   iblte[.]    bl array index at trailing edge
  c   nbl[.]      max bl array index
  c   ipan[..]    panel index corresponding to bl location
  c   isys[..]    bl newton system line number corresponding to bl location
  c   nsys        total number of lines in bl newton system
  c   itran[.]    bl array index of transition interval
  c   tforce[.]   .true. if transition is forced due to transition strip
  c
  c   va,vb[...]  diagonal and off-diagonal blocks in bl newton system
  c   vz[..]      way-off-diagonal block at te station line
  c   vm[...]     mass-influence coefficient vectors in bl newton system
  c   vdel[..]    residual and solution vectors in bl newton system
  c
  c   rmsbl       rms change from bl newton system solution
  c   rmxbl       max change from bl newton system solution
  c   imxbl       location of max change
  c   ismxbl      index of bl side containing max change
  c   vmxbl       character identifying variable with max change
  c   rlx         underrelaxation factor for newton update
  c   vaccel      parameter for accelerating bl newton system solution
  c               [any off-diagonal element < vaccel is not eliminated,
  c                which speeds up each iteration, but may increase
  c                iteration count]
  c                can be set to zero for unadulterated newton method
  c
  c   xoff,yoff   x and y offsets for windowing in qdes,gdes routines
  c   xsf ,ysf    x and y scaling factors for windowing in qdes,gdes routines
  c
  c   xgmin       airfoil grid plot limits
  c   xgmax
  c   ygmin
  c   ygmax
  c   dxyg        airfoil grid-plot annotation increment
  c   gtick       airfoil-plot tick marks size [as fraction of arc length]
  */

  //	complex<double> zcoldw, dzte, chordz, zleold, zc, zc_cn, piq, cn, eiw;

  //----- CIRCLE.INC include file for circle-plane operations
  //   NC         number of circle plane points, must be 2**n + 1
  //   MC         number of Fourier harmonics of P(w) + iQ(w)
  //   MCT        number of Fourier harmonics for which dZC/dCN are calculated
  //
  //   PI         3.1415926
  //   AGTE       trailing edge angle/pi
  //   AG0        angle of airfoil surface at first point
  //   QIM0       Q(w) offset   = Q(0)
  //   QIMOLD     Q(w) offset for old airfoil
  //   DWC        increment of circle-plane coordinate w,  DWC = 2 pi/(NC-1)
  //   WC(.)      circle plane coordinate w for Fourier operations
  //   SC(.)      normalized arc length array s(w)
  //   SCOLD(.)   normalized arc length s(w) of old airfoil
  //   XCOLD(.)   x coordinate x(w) of old airfoil
  //   YCOLD(.)   y coordinate y(w) of old airfoil
  //
  //   DZTE       trailing edge gap specified in the complex plane
  //   CHORDZ     airfoil chord specified in the complex plane
  //   ZLEOLD     leading edge of old airfoil
  //   ZCOLDW(.)  d(x+iy)/dw of old airfoil
  //   ZC(.)      complex airfoil coordinates derived from P(w) + iQ(w)
  //   ZC_CN(..)  sensitivities dZC/dCN for driving geometry constraints
  //   PIQ(.)     complex harmonic function P(w) + iQ(w)
  //   CN(.)      Fourier coefficients of P(w) + iQ(w)
  //   EIW(..)    complex number  exp(inw)  array on the unit circle

  //-----End Specific Inverse MDES-------------------------------
};
