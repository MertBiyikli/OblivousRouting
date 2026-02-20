cmake_minimum_required(VERSION 3.20)
project(ObliviousRouting
        VERSION 0.0.1
        DESCRIPTION "A comparison of tree metrics and electrical flows as oblivious routing oracles"
        LANGUAGES CXX)

# --- Options ---
option(BUILD_TESTS "Build test targets" OFF)
option(OR_ENABLE_WARNINGS "Enable strict compiler warnings" OFF)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- Dependencies ---
find_package(Boost REQUIRED COMPONENTS program_options serialization)
find_package(ortools REQUIRED CONFIG)
find_package(Eigen3 REQUIRED)

# --- Target ---
add_executable(oblivious_routing)
target_sources(oblivious_routing PRIVATE
        main.cpp
        src/utils/my_math.cpp
        src/datastructures/GraphCSR.cpp
        experiments/performance/demands/GravityModel.cpp
        experiments/performance/demands/BimodalModel.cpp
        experiments/performance/demands/GaussianModel.cpp
        experiments/performance/demands/UniformModel.cpp
        src/tree_based/utils/ultrametric_tree.cpp
        src/electrical/electrical.cpp
        src/electrical/parallel/electrical_flow_par_batches.cpp
        src/electrical/laplacian_solvers/AMGSolver.cpp
        src/lp_solver/LPSolver.cpp
        src/lp_solver/MCCF_lp_solver.cpp
)

# Includes
target_include_directories(oblivious_routing
        PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/experiments
        ${CMAKE_SOURCE_DIR}/extern/amgcl
)

# Link libs
target_link_libraries(oblivious_routing PRIVATE
        ortools::ortools
        Eigen3::Eigen
        Boost::program_options
        Boost::serialization
)


message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ---- Tests ----
if(BUILD_TESTS)
    add_subdirectory(test)
endif()
