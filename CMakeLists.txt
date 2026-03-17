cmake_minimum_required(VERSION 3.10)
project(XFoil)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

add_library(xfoil_numerics
    src/numerics/math_util.cpp
    src/numerics/spline.cpp
    src/numerics/boundary_layer_util.cpp
    src/numerics/coefficient/dissipation.cpp
    src/numerics/coefficient/skin_friction.cpp
)
target_include_directories(xfoil_numerics
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/numerics
)
target_link_libraries(xfoil_numerics
    PUBLIC
        Eigen3::Eigen
)

add_library(xfoil_model
    src/model/boundary_layer.cpp
    src/model/foil/edge.cpp
    src/model/foil/foil_shape.cpp
)
target_include_directories(xfoil_model
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/model
)
target_link_libraries(xfoil_model
    PUBLIC
        xfoil_numerics
        Eigen3::Eigen
)

add_library(xfoil_solver
    src/solver/xfoil/XFoil.cpp
    src/solver/xfoil/XFoil_aerodynamics.cpp
    src/solver/xfoil/XFoil_geometry.cpp
    src/solver/xfoil/XFoil_init.cpp
    src/solver/xfoil/XFoil_viscous.cpp
    src/solver/inviscid/InviscidSolver.cpp
    src/solver/inviscid/psi.cpp
    src/solver/boundary_layer/blsolve.cpp
    src/solver/boundary_layer/boundary_layer_workflow.cpp
    src/solver/boundary_layer/boundary_layer_store.cpp
    src/solver/boundary_layer/boundary_layer_geometry.cpp
    src/solver/boundary_layer/boundary_layer_transition.cpp
    src/solver/boundary_layer/boundary_layer_aerodynamics.cpp
    src/solver/boundary_layer/boundary_layer_mixed_mode.cpp
    src/solver/boundary_layer/boundary_layer_relaxation.cpp
    src/solver/boundary_layer/boundary_layer_runtime_state.cpp
    src/solver/boundary_layer/boundary_layer_setbl.cpp
    src/solver/boundary_layer/boundary_layer_solver_ops.cpp
    src/solver/march/march.cpp
    src/solver/march/du.cpp
    src/solver/march/ue.cpp
    src/solver/boundary_layer/boundary_layer_diff_solver.cpp
    src/solver/boundary_layer/boundary_layer_variables_solver.cpp
    src/model/foil/foil.cpp
)
target_include_directories(xfoil_solver
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/solver
        ${CMAKE_SOURCE_DIR}/src/infrastructure
)
target_link_libraries(xfoil_solver
    PUBLIC
        xfoil_model
        xfoil_numerics
        Eigen3::Eigen
)

add_executable(testExe
    apps/cli/main.cpp
)
target_include_directories(testExe
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/solver
)

target_link_libraries(testExe
    PRIVATE
        xfoil_solver
)

file(COPY sample/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/sample)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/out)
set_target_properties(testExe PROPERTIES OUTPUT_NAME "test")
set_target_properties(testExe PROPERTIES SUFFIX ".exe")
