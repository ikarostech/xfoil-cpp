#pragma once

#include "domain/boundary_layer.hpp"
#include "domain/flow_regime.hpp"

class BoundaryLayerVariablesSolver {
 public:
  BoundaryLayerVariablesSolver(double gbcon, double gccon, double ctcon);

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
  
  double gbcon_;
  double gccon_;
  double ctcon_;
};
