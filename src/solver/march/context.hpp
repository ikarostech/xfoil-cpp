#pragma once

#include <string_view>

#include "model/foil/edge.hpp"
#include "model/foil/foil.hpp"
#include "solver/boundary_layer/boundary_layer_geometry.hpp"
#include "solver/boundary_layer/boundary_layer_state.hpp"
#include "solver/boundary_layer/viscous_types.hpp"

struct MarchContextTypes {
public:
  using MixedModeStationContext = BoundaryLayerMixedModeStationContext;
  using StationReadModel = BoundaryLayerStationReadModel;
  using TrailingEdgeReadModel = BoundaryLayerTrailingEdgeReadModel;
  using EdgeVelocityFallbackMode = BoundaryLayerEdgeVelocityFallbackMode;
};

class MarchCoreContext : public MarchContextTypes {
public:
  virtual ~MarchCoreContext() = default;
};

class MarchStateContext : public virtual MarchCoreContext {
public:
  virtual ~MarchStateContext() = default;

  virtual BoundaryLayerState &mutableState() = 0;

  virtual FlowRegimeEnum
  applyFlowRegimeCandidate(FlowRegimeEnum candidate) = 0;
  virtual FlowRegimeEnum currentFlowRegime() const = 0;
  virtual blData blprv(blData data, double xsi, double ami, double cti,
                       double thi, double dsi, double dswaki,
                       double uei) const = 0;
  virtual bool blkin(BoundaryLayerState &state) = 0;
};

class MarchStationDataAccess : public virtual MarchCoreContext {
public:
  virtual ~MarchStationDataAccess() = default;

  virtual int readSideStationCount(int side) const = 0;
  virtual StationReadModel readStationModel(int side, int stationIndex) const = 0;
};

class MarchTrailingEdgeAccess : public virtual MarchCoreContext {
public:
  virtual ~MarchTrailingEdgeAccess() = default;

  virtual TrailingEdgeReadModel readTrailingEdgeModel() const = 0;
};

class MarchEventSink : public virtual MarchCoreContext {
public:
  virtual ~MarchEventSink() = default;

  virtual void emitMarchInfoLog(std::string_view message) const = 0;
};

class MarchLifecycleContext : public virtual MarchCoreContext {
public:
  virtual ~MarchLifecycleContext() = default;

  virtual FlowRegimeEnum
  determineRegimeForStation(int side, int stationIndex) const = 0;
  virtual int resetSideState(int side, const Foil &foil,
                             const StagnationResult &stagnation) = 0;
  virtual void storeStationStateCommon(int side, int stationIndex,
                                       const MixedModeStationContext &ctx) = 0;
};

class MarchRecoveryContext : public virtual MarchCoreContext {
public:
  virtual ~MarchRecoveryContext() = default;

  virtual void recoverStationAfterFailure(
      int side, int stationIndex, MixedModeStationContext &ctx, double &ami,
      EdgeVelocityFallbackMode edgeMode, int laminarAdvance) = 0;
};

class MrchueContext : public MarchStateContext,
                      public MarchStationDataAccess,
                      public MarchTrailingEdgeAccess,
                      public MarchEventSink,
                      public MarchLifecycleContext,
                      public MarchRecoveryContext {
public:
  virtual ~MrchueContext() = default;

  virtual double readNewtonRhs(int row) const = 0;
  virtual void solveMrchueDirectNewtonSystem() = 0;
  virtual void solveMrchueInverseNewtonSystem(double htarg) = 0;
  virtual void runTransitionCheckForMrchue(int side, int stationIndex,
                                           double &ami, double &cti) = 0;
  virtual bool solveTeSystemForCurrentProfiles(const Edge &edge) = 0;
  virtual double readBlCompressibilityHstinv() const = 0;
  virtual double readBlCompressibilityGm1bl() const = 0;
  virtual double readBlReynoldsReybl() const = 0;
  virtual bool blsys() = 0;
  virtual double calcHtarg(int ibl, int is, bool wake) = 0;
};

class MrchduContext : public MarchStateContext,
                      public MarchStationDataAccess,
                      public MarchEventSink,
                      public MarchLifecycleContext,
                      public MarchRecoveryContext {
public:
  virtual ~MrchduContext() = default;

  virtual bool isStartOfWake(int side, int stationIndex) = 0;
  virtual void updateSystemMatricesForStation(const Edge &edge, int side,
                                              int stationIndex,
                                              MixedModeStationContext &ctx) = 0;
  virtual void initializeFirstIterationState(int side, int stationIndex,
                                             int previousTransition,
                                             MixedModeStationContext &ctx,
                                             double &ueref,
                                             double &hkref) = 0;
  virtual void configureSimilarityRow(double ueref) = 0;
  virtual void configureViscousRow(double hkref, double ueref, double senswt,
                                   bool resetSensitivity,
                                   bool averageSensitivity, double &sens,
                                   double &sennew) = 0;
  virtual bool applyMixedModeNewtonStep(int side, int stationIndex, double &ami,
                                        MixedModeStationContext &ctx) = 0;
  virtual void checkTransitionIfNeeded(int side, int stationIndex,
                                       bool skipCheck, int laminarAdvance,
                                       double &ami) = 0;
};

class MarchContext : public MrchueContext, public MrchduContext {
public:
  ~MarchContext() override = default;
};
