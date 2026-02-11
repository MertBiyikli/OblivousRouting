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
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O2 -fno-omit-frame-pointer -fno-inline-functions")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os" CACHE STRING "" FORCE)

# Good stacks on Clang/GCC:
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-fno-omit-frame-pointer)
endif()


# only require what the main binary needs
find_package(Boost REQUIRED COMPONENTS program_options serialization)

# if you have tests that need unit_test_framework, guard it:
option(BUILD_TESTS "Build test targets" OFF)

if(BUILD_TESTS)
    find_package(Boost 1.70 REQUIRED COMPONENTS unit_test_framework)
endif()

# ---- amgcl (vendored) ----
#add_subdirectory(extern/amgcl)
# add_subdirectory(extern/xcut)


# list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/local/ortools")
# include_directories(BEFORE ${ORTOOLS_ROOT}/include)
find_package(ortools REQUIRED CONFIG)
# include_directories(BEFORE /Users/halilibrahim/local/ortools/include)

find_package(Eigen3 REQUIRED)


set(SRC_FILES
        main.cpp
        src/utils/my_math.cpp

        src/datastructures/GraphADJ.cpp
        src/datastructures/GraphADJ.h
        src/datastructures/GraphCSR.h
        src/datastructures/GraphCSR.cpp

        experiments/performance/demands/GravityModel.cpp
        experiments/performance/demands/GravityModel.h
        experiments/performance/demands/BimodalModel.cpp
        experiments/performance/demands/BimodalModel.h
        experiments/performance/demands/GaussianModel.cpp
        experiments/performance/demands/GaussianModel.h
        experiments/performance/demands/UniformModel.cpp
        experiments/performance/demands/UniformModel.h

        src/tree_based/ckr/efficient_oracle_ckr.h
        src/tree_based/ckr/efficient_oracle_ckr.cpp

        src/tree_based/frt/efficient_frt.h
        src/tree_based/frt/efficient_frt.cpp

        src/tree_based/raecke_transform.h
        src/tree_based/raecke_transform.cpp

        src/tree_based/ckr/ckr_partition.cpp
        src/tree_based/ckr/ckr_partition.h

        src/tree_based/ckr/ckr_tree_decomposer.cpp
        src/tree_based/ckr/ckr_tree_decomposer.h

        src/tree_based/ckr/utils/ultrametric_tree.h
        src/tree_based/ckr/utils/ultrametric_tree.cpp

        src/electrical/electrical_flow_optimized.cpp
        src/electrical/electrical_flow_optimized.h

        src/electrical/electrical_flow_par_batches.cpp
        src/electrical/electrical_flow_par_batches.h

        src/electrical/laplacian_solvers/AMGSolver.cpp
        src/electrical/laplacian_solvers/AMGSolver.h
        src/electrical/laplacian_solvers/AMGSolverParallel.cpp
        src/electrical/laplacian_solvers/AMGSolverParallel.h
        src/electrical/laplacian_solvers/LaplacianSolver.h
        src/electrical/laplacian_solvers/LaplacianSolverFactory.h

        #src/expanders/xcut_expander_hierarchy.cpp
        #src/expanders/xcut_expander_hierarchy.h

        src/lp_solver/LPSolver.cpp
        src/lp_solver/MCCF_lp_solver.cpp
        src/tree_based/tree_mwu.h
        src/tree_based/frt.h
        src/tree_based/fast_ckr.h
        src/tree_based/hst.h
        src/tree_based/tree_oracle.h
        src/tree_based/tree_transform.h
        src/tree_based/random_mst.h
)


# ---------------------------------------------------------------------------
# OpenMP configuration (portable: Linux + macOS)
# ---------------------------------------------------------------------------
option(USE_OPENMP "Enable OpenMP support (used in AMGSolver)" ON)

if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
endif()

add_executable(oblivious_routing ${SRC_FILES})
# Ensure OR-Tools headers are searched BEFORE anything OpenMP/Homebrew injects
get_target_property(ORT_INC ortools::ortools INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(oblivious_routing SYSTEM BEFORE PRIVATE ${ORT_INC})

target_include_directories(oblivious_routing PRIVATE ${CMAKE_SOURCE_DIR}/extern/amgcl)
# add_subdirectory(extern/xcut)
# target_link_libraries(oblivious_routing PRIVATE xcut_lib)  # or whatever the target is called


# Link OR-Tools first
target_link_libraries(oblivious_routing PRIVATE ortools::ortools)

# Then OpenMP
if(USE_OPENMP)
    target_link_libraries(oblivious_routing PRIVATE OpenMP::OpenMP_CXX)
endif()

if(APPLE)
    # Force Homebrew libomp explicitly (Intel: /usr/local, Apple Silicon: /opt/homebrew)
    if(EXISTS "/usr/local/opt/libomp/lib/libomp.dylib")
        target_link_libraries(oblivious_routing PRIVATE /usr/local/opt/libomp/lib/libomp.dylib)
        target_include_directories(oblivious_routing SYSTEM PRIVATE /usr/local/opt/libomp/include)
    elseif(EXISTS "/opt/homebrew/opt/libomp/lib/libomp.dylib")
        target_link_libraries(oblivious_routing PRIVATE /opt/homebrew/opt/libomp/lib/libomp.dylib)
        target_include_directories(oblivious_routing SYSTEM PRIVATE /opt/homebrew/opt/libomp/include)
    endif()
endif()


target_link_libraries(oblivious_routing
        PRIVATE
        Boost::boost
        Boost::program_options
        Boost::serialization
)

target_link_libraries(oblivious_routing PRIVATE Eigen3::Eigen)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ---- Tests ----
if(BUILD_TESTS)
    add_subdirectory(test)
endif()
