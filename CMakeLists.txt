cmake_minimum_required(VERSION 3.20)
project(ObliviousRouting CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Allow build type to be set via command line or CLion UI
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set sensible defaults per build type
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -fno-omit-frame-pointer" CACHE STRING "" FORCE) # add for macOs debugging
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os" CACHE STRING "" FORCE)

# Good stacks on Clang/GCC:
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-fno-omit-frame-pointer)
endif()


# only require what the main binary needs
find_package(Boost REQUIRED COMPONENTS program_options serialization)

# if you have tests that need unit_test_framework, guard it:
option(BUILD_TESTS "Build tests" ON)
if(BUILD_TESTS)
    find_package(Boost 1.70 REQUIRED COMPONENTS unit_test_framework)
endif()

# ---- amgcl (vendored) ----
add_subdirectory(extern/amgcl)

list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/local/ortools")
find_package(ortools REQUIRED CONFIG)

find_package(Eigen3 REQUIRED)




# ---- Your target ----
set(SRC_FILES
        main.cpp
        src/graph.cpp
        src/lp_solver/LPSolver.cpp
        src/lp_solver/MCCF_lp_solver.cpp
        src/tree_based/frt/raecke_frt.cpp
        src/tree_based/frt/mcct/Node.cpp
        src/tree_based/frt/mcct/mcct_derandomized_weighted_solver.cpp
        src/tree_based/frt/raecke_frt_transform.h
        src/electrical/electrical_flow_naive.cpp
        #src/electrical/AMGSolver.cpp
        src/electrical/laplacian_solvers/AMGSolver.h
        experiments/performance/demands/GravityModel.cpp
        experiments/performance/demands/GravityModel.h
        experiments/performance/Performance.cpp
        experiments/performance/Performance.h
        experiments/performance/LP_Oblivious_Ratio.cpp
        experiments/performance/LP_Oblivious_Ratio.h
        src/tree_based/frt/raecke_frt_solver.h
        src/parse_parameter.h
        experiments/performance/demands/DemandModel.h
        experiments/performance/demands/BimodalModel.cpp
        experiments/performance/demands/BimodalModel.h
        experiments/performance/demands/GaussianModel.cpp
        experiments/performance/demands/GaussianModel.h
        experiments/performance/demands/UniformModel.cpp
        experiments/performance/demands/UniformModel.h
        src/tree_based/random_mst/mst.h
        src/tree_based/random_mst/raecke_random_mst.cpp
        src/tree_based/random_mst/raecke_random_mst.h
        src/tree_based/random_mst/raecke_mst_solver.cpp
        src/tree_based/random_mst/raecke_mst_solver.h
        src/tree_based/random_mst/raecke_mst_transform.h
        src/tree_based/raecke_framework.h
        src/tree_based/raecke_base.h
        src/tree_based/raecke_transform_base.h
        experiments/performance/linear_oblivious_routing_ratio.h
        experiments/performance/linear_oblivious_routing_ratio.cpp
        src/electrical/laplacian_solvers/LaplacianSolver.h
        src/electrical/laplacian_solvers/AMGBiCGStab.cpp
        src/electrical/laplacian_solvers/AMGBiCGStab.h
        src/electrical/laplacian_solvers/LaplacianSolverFactory.h
        src/electrical/electrical_flow_optimized.cpp
        src/electrical/electrical_flow_optimized.h
        src/utils/my_math.h
        src/utils/my_math.cpp
)


option(USE_OPENMP "Enable OpenMP support (used in AMGSolver)" OFF)

# Only link OpenMP to AMGSolver if explicitly enabled
if(USE_OPENMP)
    add_executable(oblivious_routing ${SRC_FILES})

    find_package(OpenMP REQUIRED)
    message(STATUS "✅ OpenMP will be enabled for AMG solver")

    # Only build AMGSolver as separate object library with OpenMP
    add_library(amg_solver_obj OBJECT src/electrical/laplacian_solvers/AMGSolver.cpp
            src/electrical/laplacian_solvers/AMGSolver.h
            src/electrical/laplacian_solvers/LaplacianSolver.h
            src/electrical/laplacian_solvers/AMGBiCGStab.cpp
            src/electrical/laplacian_solvers/AMGBiCGStab.h
            src/electrical/laplacian_solvers/LaplacianSolverFactory.h)
    target_link_libraries(amg_solver_obj PRIVATE OpenMP::OpenMP_CXX Eigen3::Eigen)
    target_include_directories(amg_solver_obj PRIVATE
            ${CMAKE_SOURCE_DIR}/extern/amgcl
    )

    # macOS: ensure correct OpenMP runtime is linked
    if(APPLE)
        find_library(OMP_LIBRARY omp HINTS /opt/homebrew/lib /usr/local/lib)
        if(OMP_LIBRARY)
            message(STATUS "🔗 Found OpenMP runtime on macOS: ${OMP_LIBRARY}")
            target_link_libraries(oblivious_routing PRIVATE ${OMP_LIBRARY})
        else()
            message(FATAL_ERROR "❌ OpenMP library not found on macOS")
        endif()
    endif()

    # Inject compiled object into binary
    target_sources(oblivious_routing PRIVATE $<TARGET_OBJECTS:amg_solver_obj>)

else()
    message(STATUS "🚫 OpenMP disabled — compiling AMGSolver in serial mode")
    list(APPEND SRC_FILES src/electrical/laplacian_solvers/AMGSolver.cpp
            src/electrical/laplacian_solvers/AMGSolver.h
            src/electrical/laplacian_solvers/LaplacianSolver.h
            src/electrical/laplacian_solvers/AMGBiCGStab.cpp
            src/electrical/laplacian_solvers/AMGBiCGStab.h
            src/electrical/laplacian_solvers/LaplacianSolverFactory.h)
    add_executable(oblivious_routing ${SRC_FILES})
endif()


# target_link_libraries(oblivious_routing ortools::ortools)
target_include_directories(oblivious_routing PRIVATE ${CMAKE_SOURCE_DIR}/extern/amgcl)
target_link_libraries(oblivious_routing PRIVATE ortools::ortools)
target_link_libraries(oblivious_routing
        PRIVATE
        Boost::boost
        Boost::program_options
        Boost::serialization
)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
