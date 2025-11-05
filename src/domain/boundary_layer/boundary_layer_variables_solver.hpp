#pragma once

#include "domain/boundary_layer.hpp"
#include "domain/flow_regime.hpp"

class BoundaryLayerVariablesSolver {
 public:
  BoundaryLayerVariablesSolver() = default;

  blData solve(blData data, FlowRegimeEnum flowRegimeType) const;

 private:
  blData computeShapeParameters(const blData& ref,
                                FlowRegimeEnum flowRegimeType) const;
  blData computeShearCoefficients(const blData& ref,
                                  FlowRegimeEnum flowRegimeType) const;
  blData computeSkinFrictionCoefficients(const blData& ref,
                                         FlowRegimeEnum flowRegimeType) const;
  blData computeDissipation(const blData& ref,
                            FlowRegimeEnum flowRegimeType) const;
  blData computeThickness(const blData& ref,
                          FlowRegimeEnum flowRegimeType) const;

  static constexpr double kGbcon = 0.75;
  static constexpr double kGccon = 18.0;
  static constexpr double kGacon = 6.70;
  static constexpr double kCtcon = 0.5 / (kGacon * kGacon * kGbcon);
};
