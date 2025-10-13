#include <numbers>
#include <algorithm>
#include <cmath>
#include "XFoil.h"

static const FoilShape& resolveFoilShape(const Foil& foil, int node_index,
                                         int required_nodes) {
  const bool has_wake = foil.wake_shape.n > 0;
  const bool node_in_wake =
      has_wake && node_index >= 0 && node_index >= foil.foil_shape.n &&
      node_index < foil.wake_shape.n;
  const bool requires_wake_nodes =
      has_wake && required_nodes > foil.foil_shape.n;

  if (node_in_wake || requires_wake_nodes) {
    return foil.wake_shape;
  }

  return foil.foil_shape;
}

PsiResult psilin(const Foil& foil, int iNode, Vector2d point,
                 Vector2d normal_vector, bool siglin, int total_nodes,
                 const Matrix2Xd& gamu_values,
                 const Matrix2Xd& surface_vortex_values, double alfa_value,
                 double qinf_value, const VectorXd& apanel_values, bool sharp,
                 double ante, double dste, double aste) {
  const FoilShape& foil_shape = resolveFoilShape(foil, iNode, total_nodes);
  const Matrix2Xd& points = foil_shape.points;
  const VectorXd& spline_length_values = foil_shape.spline_length;

  PsiResult psi_result;

  //---- distance tolerance for determining if two points are the same
  const double seps = (spline_length_values[total_nodes - 1] -
                       spline_length_values[0]) *
                      0.00001;

  psi_result.psi = 0.0;
  psi_result.psi_ni = 0.0;
  psi_result.qtan = Vector2d::Zero();

  const int panels = total_nodes - 1;

  // panel end points
  Matrix2Xd p1 = points.leftCols(panels);
  Matrix2Xd p2(2, panels);
  p2.leftCols(panels - 1) = points.middleCols(1, panels - 1);
  p2.col(panels - 1) = points.col(panels);

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
  if (!(iNode >= 0 && iNode < total_nodes)) {
    sgn = (yy >= 0.0).select(ArrayXd::Ones(panels), ArrayXd::Constant(panels, -1.0));
  }

  ArrayXi idx = ArrayXi::LinSpaced(panels, 0, panels - 1);
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

  Matrix2Xd gam_jo = gamu_values.middleCols(0, panels);
  Matrix2Xd gam_jp(2, panels);
  gam_jp.leftCols(panels - 1) = gamu_values.middleCols(1, panels - 1);
  gam_jp.col(panels - 1) = gamu_values.col(panels);

  Matrix2Xd gsum_vector = gam_jp + gam_jo;
  Matrix2Xd gdif_vector = gam_jp - gam_jo;

  ArrayXd sv_jo =
      surface_vortex_values.row(0).segment(0, panels).array();
  ArrayXd sv_jp(panels);
  sv_jp.head(panels - 1) =
      surface_vortex_values.row(0).segment(1, panels - 1).array();
  sv_jp(panels - 1) = surface_vortex_values(0, panels);
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
      PsiResult sig_result = psisig(foil, iNode, jo, point, normal_vector,
                                    total_nodes, apanel_values);
      psi_result = PsiResult::sum(psi_result, sig_result);
    }
  }

  if ((points.col(total_nodes - 1) - points.col(0)).norm() > seps) {
    PsiResult te_result =
        psi_te(foil, iNode, normal_vector, total_nodes, apanel_values, sharp,
               ante, dste, aste, gamu_values, surface_vortex_values);
    psi_result = PsiResult::sum(psi_result, te_result);
  }

  //**** freestream terms
  Vector2d rotateVector = {cos(alfa_value), sin(alfa_value)};
  psi_result.psi += qinf_value * MathUtil::cross2(rotateVector, point);

  //---- dpsi/dn
  psi_result.psi_ni += qinf_value * MathUtil::cross2(rotateVector, normal_vector);

  psi_result.qtan.x() += qinf_value * normal_vector.y();
  psi_result.qtan.y() += -qinf_value * normal_vector.x();

  return psi_result;
}

PsiResult psisig(const Foil& foil, int iNode, int jo, Vector2d point,
                 Vector2d normal_vector, int total_nodes,
                 const VectorXd& apanel_values) {
  const FoilShape& foil_shape = resolveFoilShape(foil, iNode, total_nodes);
  const Matrix2Xd& points = foil_shape.points;

  PsiResult psi_result;
  psi_result.psi = 0;
  psi_result.psi_ni = 0;

  const int io = iNode;

  const int jp = (jo == total_nodes - 1) ? 0 : jo + 1;
  const int jm = std::max(0, jo - 1);
  const int jq = (jo >= total_nodes - 2) ? jp : jp + 1;

  const double dso = (points.col(jo) - points.col(jp)).norm();
  const double dsio = 1.0 / dso;

  const double apan = apanel_values[jo];

  const Vector2d r1 = point - points.col(jo);
  const Vector2d r2 = point - points.col(jp);
  const Vector2d s = (points.col(jp) - points.col(jo)).normalized();

  const double x1 = s.dot(r1);
  const double x2 = s.dot(r2);
  const double yy = MathUtil::cross2(s, r1);

  const double rs1 = r1.dot(r1);
  const double rs2 = r2.dot(r2);

  const double sgn = (io >= 0 && io < total_nodes) ? 1.0 : (yy >= 0.0 ? 1.0 : -1.0);

  double logr12 = 0.0;
  double t1 = 0.0;
  if (io != jo && rs1 > 0.0) {
    logr12 = log(rs1);
    t1 = atan2(sgn * x1, sgn * yy) + (0.5 - 0.5 * sgn) * std::numbers::pi;
  }

  double logr22 = 0.0;
  double t2 = 0.0;
  if (io != jp && rs2 > 0.0) {
    logr22 = log(rs2);
    t2 = atan2(sgn * x2, sgn * yy) + (0.5 - 0.5 * sgn) * std::numbers::pi;
  }

  const double x1i = s.dot(normal_vector);
  const double x2i = s.dot(normal_vector);
  const double yyi = MathUtil::cross2(s, normal_vector);

  const double x0 = 0.5 * (x1 + x2);
  const double rs0 = x0 * x0 + yy * yy;
  const double logr0 = log(rs0);
  const double theta0 =
      atan2(sgn * x0, sgn * yy) + (0.5 - 0.5 * sgn) * std::numbers::pi;

  double dxinv = 1.0 / (x1 - x0);
  double psum = x0 * (theta0 - apan) - x1 * (t1 - apan) +
                0.5 * yy * (logr12 - logr0);
  double pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) -
                 rs0 * (theta0 - apan) + (x0 - x1) * yy) *
                dxinv;

  double psx1 = -(t1 - apan);
  double psx0 = theta0 - apan;
  double psyy = 0.5 * (logr12 - logr0);

  double pdx1 =
      ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
  double pdx0 =
      ((x1 + x0) * psx0 + psum - 2.0 * x0 * (theta0 - apan) + pdif) * dxinv;
  double pdyy =
      ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - theta0))) * dxinv;

  const double dsm = (points.col(jp) - points.col(jm)).norm();
  const double dsim = 1.0 / dsm;

  const double cfac = 1.0 / (4 * std::numbers::pi);
  psi_result.dzdm[jm] += cfac * (-psum * dsim + pdif * dsim);
  psi_result.dzdm[jo] += cfac * (-psum / dso - pdif / dso);
  psi_result.dzdm[jp] += cfac * (psum * (dsio + dsim) + pdif * (dsio - dsim));

  double psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
  double pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;

  psi_result.dqdm[jm] += cfac * (-psni * dsim + pdni * dsim);
  psi_result.dqdm[jo] += cfac * (-psni / dso - pdni / dso);
  psi_result.dqdm[jp] += cfac * (psni * (dsio + dsim) + pdni * (dsio - dsim));

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

  const double dsp = (points.col(jq) - points.col(jo)).norm();
  const double dsip = 1.0 / dsp;

  psi_result.dzdm[jo] += cfac * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
  psi_result.dzdm[jp] += cfac * (psum / dso - pdif / dso);
  psi_result.dzdm[jq] += cfac * (psum * dsip + pdif * dsip);

  psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
  pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;

  psi_result.dqdm[jo] += cfac * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
  psi_result.dqdm[jp] += cfac * (psni / dso - pdni / dso);
  psi_result.dqdm[jq] += cfac * (psni * dsip + pdni * dsip);

  return psi_result;
}

PsiResult psi_te(const Foil& foil, int iNode, Vector2d normal_vector,
                 int total_nodes, const VectorXd& apanel_values, bool sharp,
                 double ante, double dste, double aste,
                 const Matrix2Xd& gamu_values,
                 const Matrix2Xd& surface_vortex_values) {
  const FoilShape& foil_shape = resolveFoilShape(foil, iNode, total_nodes);
  const Matrix2Xd& points = foil_shape.points;

  const Vector2d point = points.col(iNode);
  PsiResult psi_result;

  const double apan = apanel_values[total_nodes - 1];

  const Vector2d r1 = point - points.col(total_nodes - 1);
  const Vector2d r2 = point - points.col(0);
  const Vector2d s = (points.col(0) - points.col(total_nodes - 1)).normalized();

  const double xz1 = s.dot(r1);
  const double xz2 = s.dot(r2);
  const double yy = MathUtil::cross2(s, r1);

  const double rs1 = r1.dot(r1);
  const double rs2 = r2.dot(r2);

  const double sgn =
      (iNode >= 0 && iNode < total_nodes) ? 1.0 : (yy >= 0.0 ? 1.0 : -1.0);

  double logr12 = 0.0;
  double tz1 = 0.0;
  if (iNode != total_nodes - 1 && rs1 > 0.0) {
    logr12 = log(rs1);
    tz1 = atan2(sgn * xz1, sgn * yy) + (0.5 - 0.5 * sgn) * std::numbers::pi;
  }

  double logr22 = 0.0;
  double tz2 = 0.0;
  if (iNode != 0 && rs2 > 0.0) {
    logr22 = log(rs2);
    tz2 = atan2(sgn * xz2, sgn * yy) + (0.5 - 0.5 * sgn) * std::numbers::pi;
  }

  double scs = 1.0;
  double sds = 0.0;
  if (!sharp) {
    scs = ante / dste;
    sds = aste / dste;
  }

  const double x1i = s.dot(normal_vector);
  const double x2i = s.dot(normal_vector);
  const double yyi = MathUtil::cross2(s, normal_vector);

  const double psig = 0.5 * yy * (logr12 - logr22) +
                      xz2 * (tz2 - apan) - xz1 * (tz1 - apan);
  const double pgam = 0.5 * xz1 * logr12 - 0.5 * xz2 * logr22 + xz2 - xz1 +
                      yy * (tz1 - tz2);

  const double psigx1 = -(tz1 - apan);
  const double psigx2 = tz2 - apan;
  const double psigyy = 0.5 * (logr12 - logr22);
  const double pgamx1 = 0.5 * logr12;
  const double pgamx2 = -0.5 * logr22;
  const double pgamyy = tz1 - tz2;

  const double psigni = psigx1 * x1i + psigx2 * x2i + psigyy * yyi;
  const double pgamni = pgamx1 * x1i + pgamx2 * x2i + pgamyy * yyi;

  const Vector2d sigte_vector =
      0.5 * scs * (gamu_values.col(0) - gamu_values.col(total_nodes - 1));
  const Vector2d gamte_vector =
      -0.5 * sds * (gamu_values.col(0) - gamu_values.col(total_nodes - 1));

  const double sigte =
      0.5 * scs * (surface_vortex_values(0, 0) -
                   surface_vortex_values(0, total_nodes - 1));
  const double gamte =
      -0.5 * sds * (surface_vortex_values(0, 0) -
                    surface_vortex_values(0, total_nodes - 1));

  const double cfac = 1.0 / (2 * std::numbers::pi);
  psi_result.psi += cfac * (psig * sigte + pgam * gamte);

  psi_result.dzdg[total_nodes - 1] += -cfac * psig * scs * 0.5;
  psi_result.dzdg[0] += cfac * psig * scs * 0.5;

  psi_result.dzdg[total_nodes - 1] += cfac * pgam * sds * 0.5;
  psi_result.dzdg[0] += -cfac * pgam * sds * 0.5;

  psi_result.psi_ni += cfac * (psigni * sigte + pgamni * gamte);

  psi_result.qtan +=
      cfac * (psigni * sigte_vector + pgamni * gamte_vector);

  psi_result.dqdg[total_nodes - 1] +=
      -cfac * (psigni * 0.5 * scs - pgamni * 0.5 * sds);
  psi_result.dqdg[0] +=
      cfac * (psigni * 0.5 * scs - pgamni * 0.5 * sds);

  return psi_result;
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
PsiResult pswlin(const Foil& foil, int i, Vector2d point,
                 Vector2d normal_vector, int total_nodes, int wake_nodes,
                 const VectorXd& apanel_values) {
  const int required_nodes = total_nodes + wake_nodes;
  const FoilShape& foil_shape =
      resolveFoilShape(foil, i, required_nodes);
  const Matrix2Xd& points = foil_shape.points;

  PsiResult psi_result;
  psi_result.psi = 0.0;
  psi_result.psi_ni = 0.0;
  const int io = i;
  const int segs = wake_nodes - 1;
  if (segs <= 0)
    return psi_result;
  Matrix2Xd p0 = points.block(0, total_nodes, 2, segs);
  Matrix2Xd p1 = points.block(0, total_nodes + 1, 2, segs);
  Matrix2Xd svec = p1 - p0;
  VectorXd dso = svec.colwise().norm();
  ArrayXd dsio = dso.array().inverse();
  Matrix2Xd s = svec.array().rowwise() * dsio.transpose().array();
  Matrix2Xd r1 = point.replicate(1, segs) - p0;
  Matrix2Xd r2 = point.replicate(1, segs) - p1;
  ArrayXd x1 = (s.array() * r1.array()).colwise().sum();
  ArrayXd x2 = (s.array() * r2.array()).colwise().sum();
  ArrayXd yy = s.row(0).array() * r1.row(1).array() -
               s.row(1).array() * r1.row(0).array();
  ArrayXd rs1 = r1.colwise().squaredNorm().array();
  ArrayXd rs2 = r2.colwise().squaredNorm().array();
  ArrayXd sgn = ArrayXd::Ones(segs);
  if (!(io >= total_nodes + 1 && io <= total_nodes + wake_nodes)) {
    sgn = (yy >= 0).select(ArrayXd::Ones(segs), ArrayXd::Constant(segs, -1.0));
  }
  VectorXi jo = VectorXi::LinSpaced(segs, total_nodes, total_nodes + segs - 1);
  VectorXi jp = jo.array() + 1;
  VectorXi jm = jo.array() - 1;
  VectorXi jq = jp.array() + 1;
  jm(0) = jo(0);
  jq(segs - 1) = jp(segs - 1);
  Array<bool, Dynamic, 1> mask1 = (jo.array() != io) && (rs1 > 0.0);
  ArrayXd g1 = mask1.select(rs1.log(), 0.0);
  ArrayXd t1 = mask1.select(
      (sgn * x1).binaryExpr(
          sgn * yy, [](double a, double b) { return std::atan2(a, b); }) -
          (0.5 - 0.5 * sgn) * std::numbers::pi,
      0.0);
  Array<bool, Dynamic, 1> mask2 = (jp.array() != io) && (rs2 > 0.0);
  ArrayXd g2 = mask2.select(rs2.log(), 0.0);
  ArrayXd t2 = mask2.select(
      (sgn * x2).binaryExpr(
          sgn * yy, [](double a, double b) { return std::atan2(a, b); }) -
          (0.5 - 0.5 * sgn) * std::numbers::pi,
      0.0);
  VectorXd x1i = s.transpose() * normal_vector;
  VectorXd x2i = x1i;
  ArrayXd xsum = (x1i.array() + x2i.array()) * 0.5;
  ArrayXd yyi = s.row(0).array() * normal_vector.y() -
                s.row(1).array() * normal_vector.x();
  ArrayXd x0 = 0.5 * (x1 + x2);
  ArrayXd rs0 = x0.square() + yy.square();
  ArrayXd g0 = rs0.log();
  ArrayXd t0 = (sgn * x0).binaryExpr(sgn * yy, [](double a, double b) {
    return std::atan2(a, b);
  }) - (0.5 - 0.5 * sgn) * std::numbers::pi;
  ArrayXd dxinv = (x1 - x0).inverse();
  ArrayXd apan = apanel_values.segment(total_nodes, segs).array();
  ArrayXd psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
  ArrayXd pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan) +
                  (x0 - x1) * yy) *
                 dxinv;
  ArrayXd psx1 = -(t1 - apan);
  ArrayXd psx0 = t0 - apan;
  ArrayXd psyy = 0.5 * (g1 - g0);
  ArrayXd pdx1 =
      ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
  ArrayXd pdx0 =
      ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
  ArrayXd pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv;
  Matrix2Xd jm_p = points.block(0, total_nodes - 1, 2, segs);
  jm_p.col(0) = p0.col(0);
  Matrix2Xd jp_p = p1;
  VectorXd dsm = (jp_p - jm_p).colwise().norm();
  ArrayXd dsim = dsm.array().inverse();
  Matrix2Xd jq_p = points.block(0, total_nodes + 2, 2, segs);
  jq_p.col(segs - 1) = p1.col(segs - 1);
  VectorXd dsp = (jq_p - p0).colwise().norm();
  ArrayXd dsip = dsp.array().inverse();
  const double cfac = 1.0 / (4 * std::numbers::pi);
  ArrayXd dzdm_jm = cfac * (-psum * dsim + pdif * dsim);
  ArrayXd dzdm_jo1 = cfac * (-psum / dso.array() - pdif / dso.array());
  ArrayXd dzdm_jp1 = cfac * (psum * (dsio + dsim) + pdif * (dsio - dsim));
  ArrayXd psni = psx1 * x1i.array() + psx0 * xsum + psyy * yyi;
  ArrayXd pdni = pdx1 * x1i.array() + pdx0 * xsum + pdyy * yyi;
  ArrayXd dqdm_jm = cfac * (-psni * dsim + pdni * dsim);
  ArrayXd dqdm_jo1 = cfac * (-psni / dso.array() - pdni / dso.array());
  ArrayXd dqdm_jp1 = cfac * (psni * (dsio + dsim) + pdni * (dsio - dsim));
  ArrayXd dxinv2 = (x0 - x2).inverse();
  ArrayXd psum2 = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
  ArrayXd pdif2 = ((x0 + x2) * psum2 + rs0 * (t0 - apan) - rs2 * (t2 - apan) +
                   (x2 - x0) * yy) *
                  dxinv2;
  ArrayXd psx0_2 = -(t0 - apan);
  ArrayXd psx2 = t2 - apan;
  ArrayXd psyy2 = 0.5 * (g0 - g2);
  ArrayXd pdx0_2 =
      ((x0 + x2) * psx0_2 + psum2 + 2.0 * x0 * (t0 - apan) - pdif2) * dxinv2;
  ArrayXd pdx2 =
      ((x0 + x2) * psx2 + psum2 - 2.0 * x2 * (t2 - apan) + pdif2) * dxinv2;
  ArrayXd pdyy2 =
      ((x0 + x2) * psyy2 + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv2;
  ArrayXd dzdm_jo2 = cfac * (-psum2 * (dsip + dsio) - pdif2 * (dsip - dsio));
  ArrayXd dzdm_jp2 = cfac * (psum2 / dso.array() - pdif2 / dso.array());
  ArrayXd dzdm_jq = cfac * (psum2 * dsip + pdif2 * dsip);
  ArrayXd psni2 = psx0_2 * xsum + psx2 * x2i.array() + psyy2 * yyi;
  ArrayXd pdni2 = pdx0_2 * xsum + pdx2 * x2i.array() + pdyy2 * yyi;
  ArrayXd dqdm_jo2 = cfac * (-psni2 * (dsip + dsio) - pdni2 * (dsip - dsio));
  ArrayXd dqdm_jp2 = cfac * (psni2 / dso.array() - pdni2 / dso.array());
  ArrayXd dqdm_jq = cfac * (psni2 * dsip + pdni2 * dsip);
  VectorXd dzdm_acc = VectorXd::Zero(total_nodes + wake_nodes);
  VectorXd dqdm_acc = VectorXd::Zero(total_nodes + wake_nodes);
  dzdm_acc(jm) += dzdm_jm.matrix();
  MatrixXd dzdm_jo_mat(segs, 2);
  dzdm_jo_mat.col(0) = dzdm_jo1.matrix();
  dzdm_jo_mat.col(1) = dzdm_jo2.matrix();
  dzdm_acc(jo) += dzdm_jo_mat.rowwise().sum();
  MatrixXd dzdm_jp_mat(segs, 2);
  dzdm_jp_mat.col(0) = dzdm_jp1.matrix();
  dzdm_jp_mat.col(1) = dzdm_jp2.matrix();
  dzdm_acc(jp) += dzdm_jp_mat.rowwise().sum();
  dzdm_acc(jq) += dzdm_jq.matrix();
  dqdm_acc(jm) += dqdm_jm.matrix();
  MatrixXd dqdm_jo_mat(segs, 2);
  dqdm_jo_mat.col(0) = dqdm_jo1.matrix();
  dqdm_jo_mat.col(1) = dqdm_jo2.matrix();
  dqdm_acc(jo) += dqdm_jo_mat.rowwise().sum();
  MatrixXd dqdm_jp_mat(segs, 2);
  dqdm_jp_mat.col(0) = dqdm_jp1.matrix();
  dqdm_jp_mat.col(1) = dqdm_jp2.matrix();
  dqdm_acc(jp) += dqdm_jp_mat.rowwise().sum();
  dqdm_acc(jq) += dqdm_jq.matrix();
  psi_result.dzdm.segment(0, total_nodes + wake_nodes) = dzdm_acc;
  psi_result.dqdm.segment(0, total_nodes + wake_nodes) = dqdm_acc;
  return psi_result;
}
