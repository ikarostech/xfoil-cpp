#pragma once

#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "Eigen/Core"
#include "solver/march/context.hpp"

class XFoil;

class Marcher
{
public:
  using Matrix3x2d = Eigen::Matrix<double, 3, 2>;
  using Matrix3x2dVector = std::vector<Matrix3x2d>;
  struct MarchEvent
  {
    enum class Kind
    {
      Info,
      Failure
    };

    Kind kind = Kind::Info;
    std::string message;
  };

  static MarchEvent makeFailureEvent(std::string_view phase, int side,
                                     int stationIndex, double residual)
  {
    std::ostringstream message;
    message << "     " << phase << ": convergence failed at " << stationIndex
            << ",  side " << side << ", res =" << std::fixed
            << std::setprecision(3) << residual << "\n";
    return {MarchEvent::Kind::Failure, message.str()};
  }

  template <typename StationContext>
  struct StationMarchResult
  {
    struct RecoveryPlan
    {
      bool required = false;
      std::string_view phase;
      MarchContextTypes::EdgeVelocityFallbackMode edgeMode =
          MarchContextTypes::EdgeVelocityFallbackMode::UsePreviousStation;
      int laminarAdvance = 0;
    };

    StationContext station;
    bool converged = false;
    RecoveryPlan recovery;
    std::vector<MarchEvent> events;
  };

  template <typename StationContext, typename FailureFn, typename EventFn,
            typename CommitFn>
  static StationContext finalizeStationResult(
      StationMarchResult<StationContext> result, FailureFn &&handleFailure,
      EventFn &&publishEvent, CommitFn &&commitStation)
  {
    for (const MarchEvent &event : result.events)
    {
      publishEvent(event);
    }
    if (!result.converged && result.recovery.required)
    {
      handleFailure(result);
    }
    commitStation(result.station);
    return result.station;
  }
};
