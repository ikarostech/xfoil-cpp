#pragma once

struct FlowState {
  double alpha = 0.0;
  double qinf = 1.0;
  double referenceRe = 0.0;
  double referenceMach = 0.0;
  double currentRe = 0.0;
  double currentMach = 0.0;
  double clspec = 0.0;
  bool controlByAlpha = true;
  bool viscous = false;
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
  ReynoldsType reynoldsType = ReynoldsType::CONSTANT;
  MachType machType = MachType::CONSTANT;
};
