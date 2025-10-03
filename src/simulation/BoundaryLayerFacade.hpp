#pragma once

class XFoil;

/**
 * @brief Facade for boundary-layer related operations owned by XFoil.
 */
class BoundaryLayerFacade {
  public:
  explicit BoundaryLayerFacade(XFoil &owner);

  /**
   * @brief Reset boundary-layer state to an uninitialized state.
   */
  void reset();

  /**
   * @brief Ensure boundary-layer structures are initialized before marching.
   */
  bool ensureInitialized();

  /**
   * @brief Execute a single viscous iteration.
   */
  bool iterateOnce();

  /**
   * @brief Assemble the boundary-layer Newton system for the current state.
   */
  bool buildNewtonSystem();

  /**
   * @brief Initialize boundary layer marching using current edge velocities.
   */
  bool initializeWithCurrentUe();

  /**
   * @brief March boundary layers in displacement-thickness space.
   */
  bool marchDisplacement();

  /**
   * @brief Check if the latest march converged.
   */
  bool hasConverged() const;

  private:
  XFoil &xfoil_;
};
