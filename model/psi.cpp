#include <numbers>
#include <algorithm>
#include <cmath>
#include "../XFoil.h"

const int INDEX_START_WITH = 1;

static double cross2(const Vector2d& a, const Vector2d& b) {
  return a[0] * b[1] - a[1] * b[0];
}

PsiResult XFoil::psilin(int iNode, Vector2d point, Vector2d normal_vector, bool siglin) {
  PsiResult psi_result;

  //---- distance tolerance for determining if two points are the same
  const double seps = (spline_length[n - 1] - spline_length[0]) * 0.00001;

  psi_result.psi = 0.0;
  psi_result.psi_ni = 0.0;
  psi_result.qtan = Vector2d::Zero();

  const int panels = n - 1;

  // panel end points
  Matrix2Xd p1 = points.middleCols(INDEX_START_WITH, panels);
  Matrix2Xd p2(2, panels);
  p2.leftCols(panels - 1) = points.middleCols(INDEX_START_WITH + 1, panels - 1);
  p2.col(panels - 1) = points.col(panels + INDEX_START_WITH);

  Matrix2Xd edge = p2 - p1;
  ArrayXd dso = edge.colwise().norm();
  Array<bool, Dynamic, 1> valid_panel = (dso >= 1.0e-7);

  Matrix2Xd s = edge.array().rowwise() / dso.transpose().array();
  Matrix2Xd r1 = point.replicate(1, panels) - p1;
  Matrix2Xd r2 = point.replicate(1, panels) - p2;

  ArrayXd x1 = (s.array() * r1.array()).colwise().sum();
  ArrayXd x2 = (s.array() * r2.array()).colwise().sum();
  ArrayXd yy = s.row(0).array() * r1.row(1).array() - s.row(1).array() * r1.row(0).array();
  ArrayXd rs1 = r1.colwise().squaredNorm();
  ArrayXd rs2 = r2.colwise().squaredNorm();

  ArrayXd sgn = ArrayXd::Ones(panels);
  if (!(iNode >= 1 && iNode <= n)) {
    sgn = (yy >= 0.0).select(ArrayXd::Ones(panels), ArrayXd::Constant(panels, -1.0));
  }

  ArrayXi idx = ArrayXi::LinSpaced(panels, INDEX_START_WITH, INDEX_START_WITH + panels - 1);
  ArrayXi jp_idx = idx + 1;
  Array<bool, Dynamic, 1> mask1 = (rs1 > 0.0) && (idx != iNode);
  Array<bool, Dynamic, 1> mask2 = (rs2 > 0.0) && (jp_idx != iNode);

  ArrayXd angle1 = (sgn * x1).binaryExpr(sgn * yy, [](double a, double b) { return atan2(a, b); });
  ArrayXd angle2 = (sgn * x2).binaryExpr(sgn * yy, [](double a, double b) { return atan2(a, b); });
  ArrayXd t1 = mask1.select(angle1 + (0.5 - 0.5 * sgn) * std::numbers::pi, 0.0);
  ArrayXd t2 = mask2.select(angle2 + (0.5 - 0.5 * sgn) * std::numbers::pi, 0.0);

  ArrayXd logr12 = mask1.select(rs1.log(), 0.0);
  ArrayXd logr22 = mask2.select(rs2.log(), 0.0);

  ArrayXd dx = x1 - x2;
  Array<bool, Dynamic, 1> mask_dx = (dx != 0.0) && valid_panel;
  ArrayXd dxinv = mask_dx.select(1.0 / dx, 0.0);

  ArrayXd psis = 0.5 * x1 * logr12 - 0.5 * x2 * logr22 + x2 - x1 + yy * (t1 - t2);
  psis = valid_panel.select(psis, 0.0);
  ArrayXd psid =
      ((x1 + x2) * psis + 0.5 * (rs2 * logr22 - rs1 * logr12 + x1 * x1 - x2 * x2)) * dxinv;

  ArrayXd psx1 = 0.5 * logr12;
  ArrayXd psx2 = -0.5 * logr22;
  ArrayXd psyy = t1 - t2;

  ArrayXd pdx1 = ((x1 + x2) * psx1 + psis - x1 * logr12 - psid) * dxinv;
  ArrayXd pdx2 = ((x1 + x2) * psx2 + psis + x2 * logr22 + psid) * dxinv;
  ArrayXd pdyy = ((x1 + x2) * psyy - yy * (logr12 - logr22)) * dxinv;

  Matrix2Xd gam_jo = gamu.middleCols(0, panels);
  Matrix2Xd gam_jp(2, panels);
  gam_jp.leftCols(panels - 1) = gamu.middleCols(1, panels - 1);
  gam_jp.col(panels - 1) = gamu.col(panels);

  Matrix2Xd gsum_vector = gam_jp + gam_jo;
  Matrix2Xd gdif_vector = gam_jp - gam_jo;

  ArrayXd sv_jo = surface_vortex.row(0).segment(0, panels).array();
  ArrayXd sv_jp(panels);
  sv_jp.head(panels - 1) = surface_vortex.row(0).segment(1, panels - 1).array();
  sv_jp(panels - 1) = surface_vortex(0, panels);
  ArrayXd gsum = sv_jp + sv_jo;
  ArrayXd gdif = sv_jp - sv_jo;

  ArrayXd xi = (s.transpose() * normal_vector).array();
  ArrayXd yyi = s.row(0).array() * normal_vector.y() - s.row(1).array() * normal_vector.x();

  ArrayXd psni = psx1 * xi + psx2 * xi + psyy * yyi;
  ArrayXd pdni = pdx1 * xi + pdx2 * xi + pdyy * yyi;

  psi_result.psi +=
      (1.0 / (4.0 * std::numbers::pi)) * (psis * gsum + psid * gdif).sum();
  psi_result.psi_ni +=
      (1.0 / (4.0 * std::numbers::pi)) * (gsum * psni + gdif * pdni).sum();

  Matrix2Xd q_contrib = gsum_vector.array().rowwise() * psni.transpose().array() +
                        gdif_vector.array().rowwise() * pdni.transpose().array();
  psi_result.qtan += (1.0 / (4.0 * std::numbers::pi)) * q_contrib.rowwise().sum();

  VectorXd dzdg_jo =
      (1.0 / (4.0 * std::numbers::pi)) * (psis - psid).matrix();
  VectorXd dzdg_jp =
      (1.0 / (4.0 * std::numbers::pi)) * (psis + psid).matrix();
  psi_result.dzdg.segment(0, panels) += dzdg_jo;
  psi_result.dzdg.segment(1, panels) += dzdg_jp;

  VectorXd dqdg_jo =
      (1.0 / (4.0 * std::numbers::pi)) * (psni - pdni).matrix();
  VectorXd dqdg_jp =
      (1.0 / (4.0 * std::numbers::pi)) * (psni + pdni).matrix();
  psi_result.dqdg.segment(0, panels) += dqdg_jo;
  psi_result.dqdg.segment(1, panels) += dqdg_jp;

  if (siglin) {
    for (int jo = 0; jo < panels; ++jo) {
      if (!valid_panel(jo)) continue;
      PsiResult sig_result = psisig(iNode, jo, point, normal_vector);
      psi_result = PsiResult::sum(psi_result, sig_result);
    }
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

PsiResult XFoil::psisig(int iNode, int jo, Vector2d point, Vector2d normal_vector) {
  PsiResult psi_result;
  psi_result.psi = 0;
  psi_result.psi_ni = 0;

  int io = iNode;

  int jp = (jo == n - 1) ? 0 : jo + 1;
  int jm = std::max(0, jo - 1);
  int jq = (jo >= n - 2) ? jp : jp + 1;

  double dso = (points.col(jo + INDEX_START_WITH) - points.col(jp + INDEX_START_WITH)).norm();

  double dsio = 1.0 / dso;

  double apan = apanel[jo];

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
  if (io != jo + INDEX_START_WITH && rs1 > 0.0) {
    logr12 = log(rs1);
    t1 = atan2(sgn * x1, sgn * yy) + (0.5 - 0.5 * sgn) * std::numbers::pi;
  } else {
    logr12 = 0.0;
    t1 = 0.0;
  }
  double logr22, t2;
  if (io != jp + INDEX_START_WITH && rs2 > 0.0) {
    logr22 = log(rs2);
    t2 = atan2(sgn * x2, sgn * yy) + (0.5 - 0.5 * sgn) * std::numbers::pi;
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
  double theta0 = atan2(sgn * x0, sgn * yy) + (0.5 - 0.5 * sgn) * std::numbers::pi;

  //------- calculate source contribution to psi        for  1-0  half-panel
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
  psi_result.dzdm[jm] += (1 / (4 * std::numbers::pi)) * (-psum * dsim + pdif * dsim);
  psi_result.dzdm[jo] += (1 / (4 * std::numbers::pi)) * (-psum / dso - pdif / dso);
  psi_result.dzdm[jp] += (1 / (4 * std::numbers::pi)) * (psum * (dsio + dsim) + pdif * (dsio - dsim));

  //------- dpsi/dni
  double psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
  double pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;

  psi_result.dqdm[jm] += (1 / (4 * std::numbers::pi)) * (-psni * dsim + pdni * dsim);
  psi_result.dqdm[jo] += (1 / (4 * std::numbers::pi)) * (-psni / dso - pdni / dso);
  psi_result.dqdm[jp] += (1 / (4 * std::numbers::pi)) * (psni * (dsio + dsim) + pdni * (dsio - dsim));

  //------- calculate source contribution to psi        for  0-2  half-panel
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

  double dsp = (points.col(jq + INDEX_START_WITH) - points.col(jo + INDEX_START_WITH)).norm();
  double dsip = 1.0 / dsp;

  //------- dpsi/dm
  psi_result.dzdm[jo] += (1 / (4 * std::numbers::pi)) * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
  psi_result.dzdm[jp] += (1 / (4 * std::numbers::pi)) * (psum / dso - pdif / dso);
  psi_result.dzdm[jq] += (1 / (4 * std::numbers::pi)) * (psum * dsip + pdif * dsip);

  //------- dpsi/dni
  psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
  pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;

  psi_result.dqdm[jo] += (1 / (4 * std::numbers::pi)) * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
  psi_result.dqdm[jp] += (1 / (4 * std::numbers::pi)) * (psni / dso - pdni / dso);
  psi_result.dqdm[jq] += (1 / (4 * std::numbers::pi)) * (psni * dsip + pdni * dsip);

  return psi_result;
}

PsiResult XFoil::psi_te(int iNode, Vector2d point, Vector2d normal_vector) {
  PsiResult psi_result = PsiResult();

  //------ skip null panel
  double apan = apanel[n - 1];

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
    blData1.param.tz = atan2(sgn * blData1.param.xz, sgn * yy) + (0.5 - 0.5 * sgn) * std::numbers::pi;
  } else {
    logr12 = 0.0;
    blData1.param.tz = 0.0;
  }

  if (iNode != 1 && rs2 > 0.0) {
    logr22 = log(rs2);
    blData2.param.tz = atan2(sgn * blData2.param.xz, sgn * yy) + (0.5 - 0.5 * sgn) * std::numbers::pi;
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

  sigte = 0.5 * scs * (surface_vortex(0, 0) - surface_vortex(0, n - 1));
  gamte = -0.5 * sds * (surface_vortex(0, 0) - surface_vortex(0, n - 1));

  //---- TE panel contribution to psi
  psi_result.psi += (1 / (2 * std::numbers::pi)) * (psig * sigte + pgam * gamte);

  //---- dpsi/dgam
  psi_result.dzdg[n - 1] += -(1 / (2 * std::numbers::pi)) * psig * scs * 0.5;
  psi_result.dzdg[0] += +(1 / (2 * std::numbers::pi)) * psig * scs * 0.5;

  psi_result.dzdg[n - 1] += +(1 / (2 * std::numbers::pi)) * pgam * sds * 0.5;
  psi_result.dzdg[0] += -(1 / (2 * std::numbers::pi)) * pgam * sds * 0.5;

  //---- dpsi/dni
  psi_result.psi_ni += (1 / (2 * std::numbers::pi)) * (psigni * sigte + pgamni * gamte);

  psi_result.qtan += (1 / (2 * std::numbers::pi)) * (psigni * sigte_vector + pgamni * gamte_vector);

  psi_result.dqdg[n - 1] += -(1 / (2 * std::numbers::pi)) * (psigni * 0.5 * scs - pgamni * 0.5 * sds);
  psi_result.dqdg[0] += +(1 / (2 * std::numbers::pi)) * (psigni * 0.5 * scs - pgamni * 0.5 * sds);

  return psi_result;
}

