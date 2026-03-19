#include "solver/boundary_layer/initialization/setbl.hpp"

#include <array>
#include <cmath>
#include <sstream>

#include "solver/xfoil/XFoil.h"
#include "numerics/math_util.hpp"
#include "infrastructure/logger.hpp"
#include "solver/boundary_layer/march_access.hpp"
#include "solver/boundary_layer/initialization/setbl_access.hpp"
#include "solver/march/march.hpp"
#include "solver/march/workflow_context.hpp"

namespace {

using Eigen::VectorXd;

struct LeTeSensitivities {
  SidePair<VectorXd> ule_m;
  SidePair<VectorXd> ute_m;
};

struct EdgeVelocitySensitivityResult {
  SidePair<VectorXd> edgeVelocity;
  SidePair<VectorXd> outputEdgeVelocity;
  SidePair<int> jvte{0, 0};
  SidePair<double> dule{0.0, 0.0};
  SidePair<VectorXd> ule_m;
  SidePair<VectorXd> ute_m;
  SidePair<double> ule_a{0.0, 0.0};
};

struct SetblStation {
  VectorXd u_m;
  VectorXd d_m;
  double u_a = 0.0;
  double d_a = 0.0;
  double due = 0.0;
  double dds = 0.0;
  double xi_ule = 0.0;

  void resizeSystem(int system_size) {
    u_m = VectorXd::Zero(system_size);
    d_m = VectorXd::Zero(system_size);
  }
};

struct SetblSideData {
  SidePair<int> jvte{0, 0};
  SidePair<VectorXd> usav;
  SidePair<VectorXd> ule_m;
  SidePair<VectorXd> ute_m;
  SidePair<double> ule_a{0.0, 0.0};
  SidePair<double> dule{0.0, 0.0};

  void resizeSystem(int system_size) {
    usav.top = VectorXd::Zero(system_size);
    usav.bottom = VectorXd::Zero(system_size);
    ule_m.top = VectorXd::Zero(system_size);
    ule_m.bottom = VectorXd::Zero(system_size);
    ute_m.top = VectorXd::Zero(system_size);
    ute_m.bottom = VectorXd::Zero(system_size);
  }
};

struct StationPrimaryVars {
  double xsi = 0.0;
  double uei = 0.0;
  double thi = 0.0;
  double mdi = 0.0;
  double dsi = 0.0;
  double dswaki = 0.0;
  double ami = 0.0;
  double cti = 0.0;
};

struct TeWakeCoefficients {
  double tte = 0.0;
  double cte = 0.0;
  double dte = 0.0;
  double tte_tte1 = 0.0;
  double tte_tte2 = 0.0;
  double dte_mte1 = 0.0;
  double dte_ute1 = 0.0;
  double dte_mte2 = 0.0;
  double dte_ute2 = 0.0;
  double cte_cte1 = 0.0;
  double cte_cte2 = 0.0;
  double cte_tte1 = 0.0;
  double cte_tte2 = 0.0;
};

struct SimilarityStationCoefficients {
  VectorXd u_m1;
  VectorXd d_m1;
};

struct SideSweepInitResult {
  double u_a1 = 0.0;
  double d_a1 = 0.0;
  double due1 = 0.0;
  double dds1 = 0.0;
  double xiforc = 0.0;
};

struct StationUpdateResult {
  BoundaryLayerState state;
  VectorXd u_m2;
  VectorXd d_m2;
  double u_a2 = 0.0;
  double d_a2 = 0.0;
  double due2 = 0.0;
  double dds2 = 0.0;
};

struct TeWakeUpdateResult {
  bool isStartOfWake = false;
  VectorXd d_m1;
  double due1 = 0.0;
  double dds1 = 0.0;
  TeWakeCoefficients coeffs{};
};

struct TeWakeJacobianAdjustments {
  double vz[3][2] = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
  Eigen::Matrix<double, 3, 2> vb = Eigen::Matrix<double, 3, 2>::Zero();
};

struct StationArraysAdvanceResult {
  VectorXd u_m1;
  VectorXd d_m1;
  double u_a1 = 0.0;
  double d_a1 = 0.0;
  double due1 = 0.0;
  double dds1 = 0.0;
};

struct SetblWorkingState {
  std::array<SetblStation, 2> stations{};
  SetblSideData sides;
  double cti = 0.0;
  double ami = 0.0;
  double re_clmr = 0.0;
  double msq_clmr = 0.0;
};

class BoundaryLayerSetblAssembler {
 public:
  BoundaryLayerSetblAssembler(BoundaryLayerSetblAccess access,
                              BoundaryLayerMarchAccess march_access)
      : access_(access), marchAccess_(march_access) {}

  SetblOutputView run(SidePairRef<const BoundaryLayerSideProfiles> profiles,
                      const FlowState &analysis_state,
                      const AeroCoefficients &aero_coeffs, double acrit,
                      const Foil &foil, const StagnationResult &stagnation,
                      const Eigen::MatrixXd &dij, bool bl_initialized) {
    SetblOutputView output{};
    MarcherUe marcher_ue;
    MarcherDu marcher_du;
    SetblWorkingState state;
    initializeSetblSystemStorage(output, state.stations, state.sides);

    double current_mach = analysis_state.currentMach;
    double current_re = analysis_state.currentRe;
    initializeSetblReferenceParams(analysis_state, aero_coeffs, acrit, output,
                                   state.re_clmr, state.msq_clmr, current_mach,
                                   current_re);

    WorkflowMarchContext marchContext{marchAccess_};
    if (!bl_initialized) {
      Logger::instance().write("   Initializing bl ...\n");
      marcher_ue.mrchue(marchContext, foil, stagnation);
    }
    marcher_du.mrchdu(marchContext, foil, stagnation);

    initializeSetblProfiles(output);
    initializeSetblEdgeVelocityState(profiles, dij, output, state.sides);

    for (int side = 1; side <= 2; ++side) {
      processSetblSide(side, foil, stagnation, analysis_state.controlByAlpha,
                       dij, output, state.stations, state.sides, state.cti,
                       state.ami, state.re_clmr, state.msq_clmr);
    }

    return output;
  }

 private:
  SimilarityStationCoefficients resetSimilarityStationCoefficients(
      const VectorXd &u_m1, const VectorXd &d_m1) const {
    SimilarityStationCoefficients result;
    result.u_m1 = u_m1;
    result.d_m1 = d_m1;
    for (int js = 1; js <= 2; ++js) {
      for (int jbl = 0; jbl < access_.lattice().get(js).stationCount - 1;
           ++jbl) {
        const int jv = access_.lattice().get(js).stationToSystem[jbl];
        result.u_m1[jv] = 0.0;
        result.d_m1[jv] = 0.0;
      }
    }
    return result;
  }

  SideSweepInitResult initializeSideSweepState(
      const Foil &foil, const StagnationResult &stagnation, int side) const {
    SideSweepInitResult result;
    result.xiforc = access_.xifset(foil, stagnation, side);
    return result;
  }

  StationPrimaryVars loadStationPrimaryVars(
      int side, int station, bool station_is_wake, const SetblOutputView &output,
      double ami, double cti) const {
    StationPrimaryVars vars;
    vars.xsi = access_.lattice().get(side).arcLengthCoordinates[station];

    if (station < output.profiles.get(side).transitionIndex) {
      vars.ami = output.profiles.get(side).skinFrictionCoeff[station];
      vars.cti = cti;
    } else {
      vars.cti = output.profiles.get(side).skinFrictionCoeff[station];
      vars.ami = ami;
    }

    vars.uei = output.profiles.get(side).edgeVelocity[station];
    vars.thi = output.profiles.get(side).momentumThickness[station];
    vars.mdi = output.profiles.get(side).massFlux[station];
    vars.dsi = vars.mdi / vars.uei;

    if (station_is_wake) {
      const int wake_index = station - access_.lattice().get(side).trailingEdgeIndex;
      vars.dswaki = access_.wgap()[wake_index - 1];
    }

    return vars;
  }

  StationUpdateResult updateStationMatricesAndState(
      int side, int station, int iv, const StationPrimaryVars &vars,
      const SidePair<VectorXd> &usav, const SetblOutputView &output,
      const BoundaryLayerState &base_state, int system_size,
      const Eigen::MatrixXd &dij) {
    StationUpdateResult result;
    const double d2_m2 = 1.0 / vars.uei;
    const double d2_u2 = -vars.dsi / vars.uei;

    result.u_m2 = VectorXd::Zero(system_size);
    result.d_m2 = VectorXd::Zero(system_size);
    for (int js = 1; js <= 2; ++js) {
      for (int jbl = 0; jbl < access_.lattice().get(js).stationCount - 1;
           ++jbl) {
        const int jv = access_.lattice().get(js).stationToSystem[jbl];
        result.u_m2[jv] =
            -access_.lattice().get(side).panelInfluenceFactor[station] *
            access_.lattice().get(js).panelInfluenceFactor[jbl] *
            dij(access_.lattice().get(side).stationToPanel[station],
                access_.lattice().get(js).stationToPanel[jbl]);
        result.d_m2[jv] = d2_u2 * result.u_m2[jv];
      }
    }
    result.d_m2[iv] += d2_m2;

    result.u_a2 =
        access_.lattice().get(side).inviscidEdgeVelocityMatrix(1, station);
    result.d_a2 = d2_u2 * result.u_a2;
    result.due2 = output.profiles.get(side).edgeVelocity[station] - usav.get(side)[station];
    result.dds2 = d2_u2 * result.due2;

    result.state = base_state;
    result.state.current() = access_.blprv(
        result.state.current(), vars.xsi, vars.ami, vars.cti, vars.thi,
        vars.dsi, vars.dswaki, vars.uei);
    access_.blkin(result.state);
    return result;
  }

  void buildTransitionLog(bool station_is_transition_candidate,
                          FlowRegimeEnum flow_regime) const {
    if (station_is_transition_candidate &&
        flow_regime != FlowRegimeEnum::Transition) {
      std::stringstream ss;
      ss << "setbl: xtr???  n1=" << access_.state().station1.param.amplz
         << " n2=" << access_.state().station2.param.amplz << ":\n";
      Logger::instance().write(ss.str());
    }
  }

  TeWakeUpdateResult computeTeWakeCoefficients(
      int side, int station, const SidePair<VectorXd> &usav,
      const SidePair<VectorXd> &ute_m, const SidePair<int> &jvte,
      const VectorXd &d_m1_template, const SetblOutputView &output,
      const Edge &edge) const {
    TeWakeUpdateResult result;
    if (station != access_.lattice().get(side).trailingEdgeIndex + 1) {
      return result;
    }
    result.isStartOfWake = true;

    result.coeffs.tte =
        output.profiles.get(1)
            .momentumThickness[access_.lattice().top.trailingEdgeIndex] +
        output.profiles.get(2)
            .momentumThickness[access_.lattice().bottom.trailingEdgeIndex];
    result.coeffs.dte =
        output.profiles.get(1)
            .displacementThickness[access_.lattice().top.trailingEdgeIndex] +
        output.profiles.get(2).displacementThickness
            [access_.lattice().bottom.trailingEdgeIndex] +
        edge.ante;
    result.coeffs.cte =
        (output.profiles.get(1)
             .skinFrictionCoeff[access_.lattice().top.trailingEdgeIndex] *
             output.profiles.get(1)
                 .momentumThickness[access_.lattice().top.trailingEdgeIndex] +
         output.profiles.get(2)
                 .skinFrictionCoeff[access_.lattice().bottom.trailingEdgeIndex] *
             output.profiles.get(2).momentumThickness
                 [access_.lattice().bottom.trailingEdgeIndex]) /
        result.coeffs.tte;

    result.coeffs.tte_tte1 = 1.0;
    result.coeffs.tte_tte2 = 1.0;
    result.coeffs.dte_mte1 =
        1.0 / output.profiles.top.edgeVelocity[access_.lattice().top.trailingEdgeIndex];
    result.coeffs.dte_ute1 =
        -output.profiles.get(1)
             .displacementThickness[access_.lattice().top.trailingEdgeIndex] /
        output.profiles.top.edgeVelocity[access_.lattice().top.trailingEdgeIndex];
    result.coeffs.dte_mte2 =
        1.0 / output.profiles.bottom.edgeVelocity[access_.lattice().bottom.trailingEdgeIndex];
    result.coeffs.dte_ute2 =
        -output.profiles.get(2)
             .displacementThickness[access_.lattice().bottom.trailingEdgeIndex] /
        output.profiles.bottom.edgeVelocity[access_.lattice().bottom.trailingEdgeIndex];
    result.coeffs.cte_cte1 =
        output.profiles.get(1)
            .momentumThickness[access_.lattice().top.trailingEdgeIndex] /
        result.coeffs.tte;
    result.coeffs.cte_cte2 =
        output.profiles.get(2)
            .momentumThickness[access_.lattice().bottom.trailingEdgeIndex] /
        result.coeffs.tte;
    result.coeffs.cte_tte1 =
        (output.profiles.get(1)
             .skinFrictionCoeff[access_.lattice().top.trailingEdgeIndex] -
         result.coeffs.cte) /
        result.coeffs.tte;
    result.coeffs.cte_tte2 =
        (output.profiles.get(2)
             .skinFrictionCoeff[access_.lattice().bottom.trailingEdgeIndex] -
         result.coeffs.cte) /
        result.coeffs.tte;

    result.d_m1 = d_m1_template;
    for (int js = 1; js <= 2; ++js) {
      for (int jbl = 0; jbl < access_.lattice().get(js).stationCount - 1;
           ++jbl) {
        const int jv = access_.lattice().get(js).stationToSystem[jbl];
        result.d_m1[jv] = result.coeffs.dte_ute1 * ute_m.get(1)[jv] +
                          result.coeffs.dte_ute2 * ute_m.get(2)[jv];
      }
    }
    result.d_m1[jvte.get(1)] += result.coeffs.dte_mte1;
    result.d_m1[jvte.get(2)] += result.coeffs.dte_mte2;

    result.dds1 =
        result.coeffs.dte_ute1 *
            (output.profiles.top.edgeVelocity[access_.lattice().top.trailingEdgeIndex] -
             usav.top[access_.lattice().top.trailingEdgeIndex]) +
        result.coeffs.dte_ute2 *
            (output.profiles.bottom.edgeVelocity[access_.lattice().bottom.trailingEdgeIndex] -
             usav.bottom[access_.lattice().bottom.trailingEdgeIndex]);

    return result;
  }

  TeWakeJacobianAdjustments computeTeWakeJacobianAdjustments(
      const TeWakeCoefficients &coeffs) const {
    TeWakeJacobianAdjustments result;
    result.vz[0][0] = access_.blc().a1(0, 0) * coeffs.cte_cte1;
    result.vz[0][1] = access_.blc().a1(0, 0) * coeffs.cte_tte1 +
                      access_.blc().a1(0, 1) * coeffs.tte_tte1;
    result.vb(0, 0) = access_.blc().a1(0, 0) * coeffs.cte_cte2;
    result.vb(0, 1) = access_.blc().a1(0, 0) * coeffs.cte_tte2 +
                      access_.blc().a1(0, 1) * coeffs.tte_tte2;

    result.vz[1][0] = access_.blc().a1(1, 0) * coeffs.cte_cte1;
    result.vz[1][1] = access_.blc().a1(1, 0) * coeffs.cte_tte1 +
                      access_.blc().a1(1, 1) * coeffs.tte_tte1;
    result.vb(1, 0) = access_.blc().a1(1, 0) * coeffs.cte_cte2;
    result.vb(1, 1) = access_.blc().a1(1, 0) * coeffs.cte_tte2 +
                      access_.blc().a1(1, 1) * coeffs.tte_tte2;

    result.vz[2][0] = access_.blc().a1(2, 0) * coeffs.cte_cte1;
    result.vz[2][1] = access_.blc().a1(2, 0) * coeffs.cte_tte1 +
                      access_.blc().a1(2, 1) * coeffs.tte_tte1;
    result.vb(2, 0) = access_.blc().a1(2, 0) * coeffs.cte_cte2;
    result.vb(2, 1) = access_.blc().a1(2, 0) * coeffs.cte_tte2 +
                      access_.blc().a1(2, 1) * coeffs.tte_tte2;
    return result;
  }

  EdgeVelocitySensitivityResult prepareEdgeVelocityAndSensitivities(
      SidePairRef<const BoundaryLayerSideProfiles> profiles,
      const Eigen::MatrixXd &dij, int nsys) const {
    EdgeVelocitySensitivityResult result;

    result.edgeVelocity = access_.ueset(dij);
    result.outputEdgeVelocity.top = profiles.top.edgeVelocity;
    result.outputEdgeVelocity.bottom = profiles.bottom.edgeVelocity;

    result.jvte.top = access_.lattice().top.stationToSystem
        [access_.lattice().top.trailingEdgeIndex];
    result.jvte.bottom = access_.lattice().bottom.stationToSystem
        [access_.lattice().bottom.trailingEdgeIndex];

    result.dule.top = result.outputEdgeVelocity.top[0] - result.edgeVelocity.top[0];
    result.dule.bottom =
        result.outputEdgeVelocity.bottom[0] - result.edgeVelocity.bottom[0];

    const auto le_te_sensitivities = computeLeTeSensitivities(
        access_.lattice().get(1).stationToPanel[0],
        access_.lattice().get(2).stationToPanel[0],
        access_.lattice().get(1)
            .stationToPanel[access_.lattice().top.trailingEdgeIndex],
        access_.lattice().get(2)
            .stationToPanel[access_.lattice().bottom.trailingEdgeIndex],
        nsys, dij);
    result.ule_m = le_te_sensitivities.ule_m;
    result.ute_m = le_te_sensitivities.ute_m;

    result.ule_a.top = access_.lattice().get(1).inviscidEdgeVelocityMatrix(1, 0);
    result.ule_a.bottom = access_.lattice().get(2).inviscidEdgeVelocityMatrix(1, 0);
    return result;
  }

  void assembleBlJacobianForStation(
      int iv, int nsys, const std::array<SetblStation, 2> &setblStations,
      const SetblSideData &setblSides, bool controlByAlpha, double re_clmr,
      double msq_clmr, SetblOutputView &output) {
    output.bl_newton_system.vb[iv] = access_.blc().a1.block(0, 0, 3, 2);
    output.bl_newton_system.va[iv] = access_.blc().a2.block(0, 0, 3, 2);

    Eigen::Matrix<double, 3, 4> A;
    A.col(0) = access_.blc().a1.col(3).head<3>();
    A.col(1) = access_.blc().a1.col(2).head<3>();
    A.col(2) = access_.blc().a2.col(3).head<3>();
    A.col(3) = access_.blc().a2.col(2).head<3>();

    Eigen::Matrix<double, 4, 2> B;
    B << setblStations[0].due, setblStations[0].u_a, setblStations[0].dds,
        setblStations[0].d_a, setblStations[1].due, setblStations[1].u_a,
        setblStations[1].dds, setblStations[1].d_a;
    const Eigen::Vector3d ax =
        (access_.blc().a1.col(4) + access_.blc().a2.col(4) + access_.blc().d_xi)
            .head<3>();
    const Eigen::RowVector2d xi =
        setblStations[0].xi_ule *
            Eigen::RowVector2d(setblSides.dule.get(1), setblSides.ule_a.get(1)) +
        setblStations[1].xi_ule *
            Eigen::RowVector2d(setblSides.dule.get(2), setblSides.ule_a.get(2));
    output.bl_newton_system.vdel[iv] = A * B + ax * xi;

    for (int jv = 1; jv < nsys; ++jv) {
      const Eigen::Vector4d m(setblStations[0].u_m(jv), setblStations[0].d_m(jv),
                              setblStations[1].u_m(jv), setblStations[1].d_m(jv));
      const double xi_m = setblStations[0].xi_ule * setblSides.ule_m.get(1)(jv) +
                          setblStations[1].xi_ule * setblSides.ule_m.get(2)(jv);
      const Eigen::Vector3d vm = A * m + ax * xi_m;
      output.bl_newton_system.vm.at(0, jv, iv) = vm[0];
      output.bl_newton_system.vm.at(1, jv, iv) = vm[1];
      output.bl_newton_system.vm.at(2, jv, iv) = vm[2];
    }

    if (controlByAlpha) {
      output.bl_newton_system.vdel[iv].col(1).head<3>() =
          access_.blc().d_re.head(3) * re_clmr +
          access_.blc().d_msq.head(3) * msq_clmr;
    }
  }

  StationArraysAdvanceResult advanceStationArrays(const VectorXd &u_m2,
                                                  const VectorXd &d_m2,
                                                  double u_a2, double d_a2,
                                                  double due2,
                                                  double dds2) const {
    StationArraysAdvanceResult result;
    result.u_m1 = u_m2;
    result.d_m1 = d_m2;
    result.u_a1 = u_a2;
    result.d_a1 = d_a2;
    result.due1 = due2;
    result.dds1 = dds2;
    return result;
  }

  BoundaryLayerReferenceParams computeBlReferenceParams(
      const FlowState &analysis_state, const AeroCoefficients &aero_coeffs,
      double acrit) const {
    return BoundaryLayerPhysics::buildReferenceParams(analysis_state,
                                                      aero_coeffs, acrit);
  }

  void initializeSetblReferenceParams(
      const FlowState &analysis_state, const AeroCoefficients &aero_coeffs,
      double acrit, SetblOutputView &output, double &re_clmr,
      double &msq_clmr, double &currentMach, double &currentRe) {
    const auto reference_params =
        computeBlReferenceParams(analysis_state, aero_coeffs, acrit);
    currentMach = reference_params.currentMach;
    currentRe = reference_params.currentRe;
    re_clmr = reference_params.re_clmr;
    msq_clmr = reference_params.msq_clmr;
    output.blCompressibility = reference_params.blCompressibility;
    output.blReynolds = reference_params.blReynolds;
    output.blTransition.amcrit = reference_params.amcrit;
    access_.blCompressibility() = output.blCompressibility;
    access_.blReynolds() = output.blReynolds;
    access_.blTransition().amcrit = output.blTransition.amcrit;
  }

  void initializeSetblSystemStorage(
      SetblOutputView &output, std::array<SetblStation, 2> &stations,
      SetblSideData &sideData) const {
    const auto zero = Eigen::Matrix<double, 3, 2>::Zero();
    output.bl_newton_system.vm.resize(access_.systemSize());
    output.bl_newton_system.va.resize(access_.systemSize(), zero);
    output.bl_newton_system.vb.resize(access_.systemSize(), zero);
    output.bl_newton_system.vdel.resize(access_.systemSize(), zero);
    stations[0].resizeSystem(access_.systemSize());
    stations[1].resizeSystem(access_.systemSize());
    sideData.resizeSystem(access_.systemSize());
  }

  void initializeSetblProfiles(SetblOutputView &output) const {
    output.profiles.top = access_.lattice().top.profiles;
    output.profiles.bottom = access_.lattice().bottom.profiles;
  }

  void initializeSetblEdgeVelocityState(
      SidePairRef<const BoundaryLayerSideProfiles> profiles,
      const Eigen::MatrixXd &dij, SetblOutputView &output,
      SetblSideData &sideData) const {
    const auto edge_result =
        prepareEdgeVelocityAndSensitivities(profiles, dij, access_.systemSize());
    sideData.usav = edge_result.edgeVelocity;
    sideData.jvte = edge_result.jvte;
    sideData.dule = edge_result.dule;
    sideData.ule_m = edge_result.ule_m;
    sideData.ute_m = edge_result.ute_m;
    sideData.ule_a = edge_result.ule_a;
    output.profiles.top.edgeVelocity = edge_result.outputEdgeVelocity.top;
    output.profiles.bottom.edgeVelocity = edge_result.outputEdgeVelocity.bottom;
  }

  void processSetblSide(
      int side, const Foil &foil, const StagnationResult &stagnation,
      bool controlByAlpha, const Eigen::MatrixXd &dij, SetblOutputView &output,
      std::array<SetblStation, 2> &stations, SetblSideData &sideData,
      double &cti, double &ami, double re_clmr, double msq_clmr) {
    const auto similarity_coeffs =
        resetSimilarityStationCoefficients(stations[0].u_m, stations[0].d_m);
    stations[0].u_m = similarity_coeffs.u_m1;
    stations[0].d_m = similarity_coeffs.d_m1;

    const auto sweep_init = initializeSideSweepState(foil, stagnation, side);
    stations[0].u_a = sweep_init.u_a1;
    stations[0].d_a = sweep_init.d_a1;
    stations[0].due = sweep_init.due1;
    stations[0].dds = sweep_init.dds1;
    output.blTransition.xiforc = sweep_init.xiforc;
    access_.blTransition().xiforc = output.blTransition.xiforc;

    for (int station = 0; station < access_.lattice().get(side).stationCount - 1;
         ++station) {
      const int iv = access_.lattice().get(side).stationToSystem[station];
      const bool station_is_wake =
          station > access_.lattice().get(side).trailingEdgeIndex;
      const bool station_is_transition_candidate =
          station == output.profiles.get(side).transitionIndex;

      output.flowRegime = access_.determineRegimeForStation(side, station);
      access_.flowRegime() = output.flowRegime;

      const auto vars = loadStationPrimaryVars(side, station, station_is_wake,
                                               output, ami, cti);
      ami = vars.ami;
      cti = vars.cti;

      auto station_update = updateStationMatricesAndState(
          side, station, iv, vars, sideData.usav, output, access_.state(),
          stations[1].u_m.size(), dij);
      stations[1].u_m = station_update.u_m2;
      stations[1].d_m = station_update.d_m2;
      stations[1].u_a = station_update.u_a2;
      stations[1].d_a = station_update.d_a2;
      stations[1].due = station_update.due2;
      stations[1].dds = station_update.dds2;
      access_.state() = station_update.state;

      if (station_is_transition_candidate) {
        access_.runTransitionCheck();
        ami = access_.state().station2.param.amplz;
      }
      buildTransitionLog(station_is_transition_candidate, output.flowRegime);

      const auto te_update = computeTeWakeCoefficients(
          side, station, sideData.usav, sideData.ute_m, sideData.jvte,
          stations[0].d_m, output, foil.edge);
      if (te_update.isStartOfWake) {
        access_.tesys(access_.lattice().top.profiles,
                      access_.lattice().bottom.profiles, foil.edge);
        stations[0].d_m = te_update.d_m1;
        stations[0].due = te_update.due1;
        stations[0].dds = te_update.dds1;
      } else {
        access_.blsys();
      }

      output.profiles.get(side).skinFrictionCoeffHistory[station] =
          access_.state().station2.cqz.scalar;

      if (side == 1) {
        stations[0].xi_ule = stagnation.sst_go;
        stations[1].xi_ule = -stagnation.sst_gp;
      } else {
        stations[0].xi_ule = -stagnation.sst_go;
        stations[1].xi_ule = stagnation.sst_gp;
      }

      assembleBlJacobianForStation(iv, access_.systemSize(), stations, sideData,
                                   controlByAlpha, re_clmr, msq_clmr, output);

      if (te_update.isStartOfWake) {
        const auto te_jacobian =
            computeTeWakeJacobianAdjustments(te_update.coeffs);
        for (int row = 0; row < 3; ++row) {
          output.bl_newton_system.vz[row][0] = te_jacobian.vz[row][0];
          output.bl_newton_system.vz[row][1] = te_jacobian.vz[row][1];
        }
        output.bl_newton_system.vb[iv] = te_jacobian.vb;
      }

      if (output.flowRegime == FlowRegimeEnum::Transition) {
        output.profiles.get(side).transitionIndex = station;
        access_.lattice().get(side).profiles.transitionIndex = station;
        output.flowRegime = FlowRegimeEnum::Turbulent;
        access_.flowRegime() = output.flowRegime;
      }

      if (station == access_.lattice().get(side).trailingEdgeIndex) {
        output.flowRegime = FlowRegimeEnum::Wake;
        access_.flowRegime() = output.flowRegime;
        access_.solveWakeState();
        access_.blmid(FlowRegimeEnum::Wake);
      }

      const auto advance = advanceStationArrays(
          stations[1].u_m, stations[1].d_m, stations[1].u_a, stations[1].d_a,
          stations[1].due, stations[1].dds);
      stations[0].u_m = advance.u_m1;
      stations[0].d_m = advance.d_m1;
      stations[0].u_a = advance.u_a1;
      stations[0].d_a = advance.d_a1;
      stations[0].due = advance.due1;
      stations[0].dds = advance.dds1;
      access_.state().stepbl();
    }
  }

  LeTeSensitivities computeLeTeSensitivities(int ile1, int ile2, int ite1,
                                             int ite2, int nsys,
                                             const Eigen::MatrixXd &dij) const {
    LeTeSensitivities sensitivities;
    sensitivities.ule_m.top = VectorXd::Zero(nsys);
    sensitivities.ule_m.bottom = VectorXd::Zero(nsys);
    sensitivities.ute_m.top = VectorXd::Zero(nsys);
    sensitivities.ute_m.bottom = VectorXd::Zero(nsys);
    for (int js = 1; js <= 2; ++js) {
      for (int jbl = 0; jbl < access_.lattice().get(js).stationCount - 1;
           ++jbl) {
        const int panel_index = access_.lattice().get(js).stationToPanel[jbl];
        const int system_index = access_.lattice().get(js).stationToSystem[jbl];
        const double panel_factor =
            access_.lattice().get(js).panelInfluenceFactor[jbl];
        sensitivities.ule_m.top[system_index] =
            -access_.lattice().top.panelInfluenceFactor[0] * panel_factor *
            dij(ile1, panel_index);
        sensitivities.ule_m.bottom[system_index] =
            -access_.lattice().bottom.panelInfluenceFactor[0] * panel_factor *
            dij(ile2, panel_index);
        sensitivities.ute_m.top[system_index] =
            -access_.lattice().top
                 .panelInfluenceFactor[access_.lattice().top.trailingEdgeIndex] *
            panel_factor * dij(ite1, panel_index);
        sensitivities.ute_m.bottom[system_index] =
            -access_.lattice().bottom.panelInfluenceFactor
                 [access_.lattice().bottom.trailingEdgeIndex] *
            panel_factor * dij(ite2, panel_index);
      }
    }
    return sensitivities;
  }

  BoundaryLayerSetblAccess access_;
  BoundaryLayerMarchAccess marchAccess_;
};

} // namespace

SetblOutputView runBoundaryLayerSetbl(
    BoundaryLayerWorkflow &workflow,
    SidePairRef<const BoundaryLayerSideProfiles> profiles,
    const FlowState &analysis_state, const AeroCoefficients &aero_coeffs,
    double acrit, const Foil &foil, const StagnationResult &stagnation,
    const Eigen::MatrixXd &dij, bool bl_initialized) {
  BoundaryLayerSetblAssembler assembler{
      makeBoundaryLayerSetblAccess(workflow),
      makeBoundaryLayerMarchAccess(workflow)};
  return assembler.run(profiles, analysis_state, aero_coeffs, acrit, foil,
                       stagnation, dij, bl_initialized);
}
