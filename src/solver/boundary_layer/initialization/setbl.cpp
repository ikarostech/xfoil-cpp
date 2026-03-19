#include "solver/boundary_layer/initialization/setbl.hpp"

#include <array>
#include <cmath>
#include <sstream>

#include "solver/xfoil/XFoil.h"
#include "numerics/math_util.hpp"
#include "infrastructure/logger.hpp"
#include "solver/boundary_layer/march_access.hpp"
#include "solver/boundary_layer/initialization/setbl_access.hpp"
#include "solver/boundary_layer/initialization/setbl_edge_velocity.hpp"
#include "solver/boundary_layer/initialization/setbl_jacobian.hpp"
#include "solver/boundary_layer/initialization/setbl_side_sweep.hpp"
#include "solver/boundary_layer/initialization/setbl_station_update.hpp"
#include "solver/boundary_layer/initialization/setbl_trailing_edge.hpp"
#include "solver/march/march.hpp"
#include "solver/march/workflow_context.hpp"

namespace {

using Eigen::VectorXd;

using SetblStation = setbl_jacobian::StationState;

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

using StationPrimaryVars = setbl_station_update::StationPrimaryVars;

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

using StationUpdateResult = setbl_station_update::StationUpdateResult;

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
  void resizeStation(SetblStation &station, int system_size) const {
    station.u_m = VectorXd::Zero(system_size);
    station.d_m = VectorXd::Zero(system_size);
  }

  SimilarityStationCoefficients resetSimilarityStationCoefficients(
      const VectorXd &u_m1, const VectorXd &d_m1) const {
    SimilarityStationCoefficients result;
    result.u_m1 = u_m1;
    result.d_m1 = d_m1;
    for (int js = 1; js <= 2; ++js) {
      for (int jbl = 0; jbl < access_.stationCount(js) - 1; ++jbl) {
        const int jv = access_.stationToSystem(js, jbl);
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
    vars.xsi = access_.arcLengthCoordinate(side, station);

    if (station < output.profiles.get(side).transitionIndex) {
      vars.ami = output.profiles.get(side).skinFrictionCoeff[station];
      vars.cti = cti;
    } else {
      vars.cti = output.profiles.get(side).skinFrictionCoeff[station];
      vars.ami = ami;
    }

    vars.uei = output.profiles.get(side).edgeVelocity[station];
    vars.thi = output.profiles.get(side).momentumThickness[station];
    vars.dsi = output.profiles.get(side).massFlux[station] / vars.uei;

    if (station_is_wake) {
      vars.dswaki = access_.wakeGapAt(side, station);
    }

    return vars;
  }

  StationUpdateResult updateStationMatricesAndState(
      int side, int station, int iv, const StationPrimaryVars &vars,
      const SidePair<VectorXd> &usav, const SetblOutputView &output,
      const BoundaryLayerState &base_state, int system_size,
      const Eigen::MatrixXd &dij) {
    return setbl_station_update::SetblStationUpdateOps::update(
        access_, side, station, iv, vars, usav, output.profiles, base_state,
        system_size, dij);
  }

  void buildTransitionLog(bool station_is_transition_candidate,
                          FlowRegimeEnum flow_regime) const {
    if (station_is_transition_candidate &&
        flow_regime != FlowRegimeEnum::Transition) {
      std::stringstream ss;
      ss << "setbl: xtr???  n1=" << access_.state().station1.param.amplz
         << " n2=" << access_.currentAmplification() << ":\n";
      Logger::instance().write(ss.str());
    }
  }

  setbl_te::TeWakeUpdateResult computeTeWakeCoefficients(
      int side, int station, const SidePair<VectorXd> &usav,
      const SidePair<VectorXd> &ute_m, const SidePair<int> &jvte,
      const VectorXd &d_m1_template, const SetblOutputView &output,
      const Edge &edge) const {
    return setbl_te::SetblTrailingEdgeOps::computeWakeUpdate(
        access_, side, station, usav, ute_m, jvte, d_m1_template,
        output.profiles, edge);
  }

  setbl_te::TeWakeJacobianAdjustments computeTeWakeJacobianAdjustments(
      const setbl_te::TeWakeCoefficients &coeffs) const {
    return setbl_te::SetblTrailingEdgeOps::computeJacobianAdjustments(access_,
                                                                      coeffs);
  }

  setbl_edge_velocity::EdgeVelocitySensitivityResult
  prepareEdgeVelocityAndSensitivities(
      SidePairRef<const BoundaryLayerSideProfiles> profiles,
      const Eigen::MatrixXd &dij, int nsys) const {
    return setbl_edge_velocity::SetblEdgeVelocityOps::prepare(access_, profiles,
                                                              dij, nsys);
  }

  void assembleBlJacobianForStation(
      int iv, int nsys, const std::array<SetblStation, 2> &setblStations,
      const SetblSideData &setblSides, bool controlByAlpha, double re_clmr,
      double msq_clmr, SetblOutputView &output) {
    std::array<setbl_jacobian::StationState, 2> stations;
    for (int i = 0; i < 2; ++i) {
      stations[i].u_m = setblStations[i].u_m;
      stations[i].d_m = setblStations[i].d_m;
      stations[i].u_a = setblStations[i].u_a;
      stations[i].d_a = setblStations[i].d_a;
      stations[i].due = setblStations[i].due;
      stations[i].dds = setblStations[i].dds;
      stations[i].xi_ule = setblStations[i].xi_ule;
    }

    setbl_jacobian::SideState sideState;
    sideState.ule_m = setblSides.ule_m;
    sideState.ule_a = setblSides.ule_a;
    sideState.dule = setblSides.dule;

    setbl_jacobian::SetblJacobianOps::assembleForStation(
        access_, iv, nsys, stations, sideState, controlByAlpha, re_clmr,
        msq_clmr, output.bl_newton_system);
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
    resizeStation(stations[0], access_.systemSize());
    resizeStation(stations[1], access_.systemSize());
    sideData.resizeSystem(access_.systemSize());
  }

  void initializeSetblProfiles(SetblOutputView &output) const {
    access_.copyProfilesTo(output.profiles);
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

    for (int station = 0; station < access_.stationCount(side) - 1;
         ++station) {
      const int iv = access_.stationToSystem(side, station);
      const bool station_is_wake = station > access_.trailingEdgeIndex(side);
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
      access_.replaceState(station_update.state);

      if (station_is_transition_candidate) {
        access_.runTransitionCheck();
        ami = access_.currentAmplification();
      }
      buildTransitionLog(station_is_transition_candidate, output.flowRegime);

      const auto te_update = computeTeWakeCoefficients(
          side, station, sideData.usav, sideData.ute_m, sideData.jvte,
          stations[0].d_m, output, foil.edge);
      if (te_update.isStartOfWake) {
        access_.solveTeSystemForCurrentProfiles(foil.edge);
        stations[0].d_m = te_update.d_m1;
        stations[0].due = te_update.due1;
        stations[0].dds = te_update.dds1;
      } else {
        access_.blsys();
      }

      output.profiles.get(side).skinFrictionCoeffHistory[station] =
          access_.currentSkinFrictionHistory();

      setbl_side_sweep::SetblSideSweepOps::assignXiUle(side, stagnation,
                                                       stations);

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

      setbl_side_sweep::SetblSideSweepOps::updateRegimeAfterStation(
          access_, side, station, output);
      setbl_side_sweep::SetblSideSweepOps::advanceStationWindow(stations[1],
                                                                stations[0],
                                                                access_);
    }
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
