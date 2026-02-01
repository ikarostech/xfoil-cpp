cmake_minimum_required(VERSION 3.10)
project(XFoil)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

add_library(xfoil_core
    src/core/math_util.cpp
    src/core/spline.cpp
    src/core/boundary_layer_util.cpp
)
target_include_directories(xfoil_core
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/core
)
target_link_libraries(xfoil_core
    PUBLIC
        Eigen3::Eigen
)

add_library(xfoil_domain
    src/domain/boundary_layer.cpp
    src/domain/boundary_layer/boundary_layer_diff_solver.cpp
    src/domain/boundary_layer/boundary_layer_variables_solver.cpp
    src/domain/coefficient/dissipation.cpp
    src/domain/coefficient/skin_friction.cpp
    src/domain/foil/edge.cpp
    src/domain/foil/foil_shape.cpp
)
target_include_directories(xfoil_domain
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/domain
)
target_link_libraries(xfoil_domain
    PUBLIC
        xfoil_core
        Eigen3::Eigen
)

add_library(xfoil_simulation
    src/simulation/XFoil.cpp
    src/simulation/XFoil_aerodynamics.cpp
    src/simulation/InviscidSolver.cpp
    src/simulation/Blsolve.cpp
    src/simulation/XFoil_geometry.cpp
    src/simulation/XFoil_init.cpp
    src/simulation/XFoil_viscous.cpp
    src/simulation/BoundaryLayer.cpp
    src/simulation/BoundaryLayer_geometry.cpp
    src/simulation/BoundaryLayerStore.cpp
    src/simulation/BoundaryLayer_march.cpp
    src/simulation/BoundaryLayer_transition.cpp
    src/simulation/psi.cpp
    src/domain/foil/foil.cpp
)
target_include_directories(xfoil_simulation
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/simulation
        ${CMAKE_SOURCE_DIR}/src/infrastructure
)
target_link_libraries(xfoil_simulation
    PUBLIC
        xfoil_domain
        xfoil_core
        Eigen3::Eigen
)

add_executable(testExe
    apps/cli/main.cpp
)
target_include_directories(testExe
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/simulation
)

target_link_libraries(testExe
    PRIVATE
        xfoil_simulation
)

file(COPY sample/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/sample)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/out)
set_target_properties(testExe PROPERTIES OUTPUT_NAME "test")
set_target_properties(testExe PROPERTIES SUFFIX ".exe")
