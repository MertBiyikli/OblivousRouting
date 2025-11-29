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

# list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/local/ortools")
# include_directories(BEFORE ${ORTOOLS_ROOT}/include)
find_package(ortools REQUIRED CONFIG)
# include_directories(BEFORE /Users/halilibrahim/local/ortools/include)

find_package(Eigen3 REQUIRED)




# ---- Your target ----
set(SRC_FILES
        main.cpp
        src/datastructures/graph.cpp
        src/lp_solver/LPSolver.cpp
        src/lp_solver/MCCF_lp_solver.cpp
        src/tree_based/random_mst/raecke_mst_solver.cpp
        src/tree_based/random_mst/raecke_random_mst.cpp
        src/tree_based/frt/raecke_frt.cpp
        src/tree_based/frt/mcct/Node.cpp
        src/tree_based/frt/mcct/mcct_derandomized_weighted_solver.cpp
        src/electrical/electrical_flow_naive.cpp
        src/electrical/electrical_flow_optimized.cpp
        src/electrical/electrical_flow_optimized.h
        src/electrical/electrical_flow_parallel.cpp
        src/electrical/electrical_flow_parallel.h
        src/electrical/laplacian_solvers/AMGSolver.cpp
        src/electrical/laplacian_solvers/AMGSolver.h
        src/electrical/laplacian_solvers/LaplacianSolver.h
        src/electrical/laplacian_solvers/LaplacianSolverFactory.h
        experiments/performance/Performance.cpp
        experiments/performance/Performance.h
        experiments/performance/LP_Oblivious_Ratio.cpp
        experiments/performance/LP_Oblivious_Ratio.h
        experiments/performance/linear_oblivious_routing_ratio.cpp
        experiments/performance/linear_oblivious_routing_ratio.h
        experiments/performance/demands/GravityModel.cpp
        experiments/performance/demands/GravityModel.h
        experiments/performance/demands/BimodalModel.cpp
        experiments/performance/demands/BimodalModel.h
        experiments/performance/demands/GaussianModel.cpp
        experiments/performance/demands/GaussianModel.h
        experiments/performance/demands/UniformModel.cpp
        experiments/performance/demands/UniformModel.h
        src/utils/my_math.cpp
        src/utils/my_math.h
        src/electrical/laplacian_solvers/AMGSolverParallel.cpp
        src/electrical/laplacian_solvers/AMGSolverParallel.h
        src/electrical/electrical_flow_parallel_on_the_fly.cpp
        src/electrical/electrical_flow_parallel_on_the_fly.h
        src/tree_based/ckr/ckr_partition.cpp
        src/tree_based/ckr/ckr_partition.h
        src/utils/priority_queue.h
        src/utils/lca_datstructure.h
        src/tree_based/ckr/ckr_tree_decomposer.cpp
        src/tree_based/ckr/ckr_tree_decomposer.h
        src/tree_based/ckr/raecke_ckr_transform.h
        src/tree_based/ckr/raecke_ckr_solver.h
        # src/tree_based/ckr/raecke_ckr_optimized.cpp
        # src/tree_based/ckr/raecke_ckr_optimized.h
        src/tree_based/optimized_versions/raecke_optimized_base.h
        src/tree_based/optimized_versions/frt/raecke_frt_opt.h
        src/tree_based/optimized_versions/frt/raecke_frt_transform_opt.h
        # src/tree_based/optimized_versions/frt/raecke_frt_opt.cpp
        src/datastructures/graph_csr.h
        src/datastructures/graph_csr.cpp
        # src/tree_based/optimized_versions/frt/mcct_optimized.cpp
        # src/tree_based/optimized_versions/frt/mcct_optimized.h
        src/tree_based/ckr/utils/ultrametric_tree.h
        src/tree_based/ckr/utils/ultrametric_tree.cpp
        src/tree_based/ckr/optimized/efficient_oracle_ckr.h
        src/tree_based/ckr/optimized/efficient_oracle_ckr.cpp
        src/tree_based/ckr/optimized/efficient_raecke_ckr.h
        src/tree_based/ckr/optimized/efficient_raecke_ckr_transform.h
        src/tree_based/ckr/optimized/efficient_raecke_ckr_transform.cpp

)


# ---------------------------------------------------------------------------
# OpenMP configuration (portable: Linux + macOS)
# ---------------------------------------------------------------------------
option(USE_OPENMP "Enable OpenMP support (used in AMGSolver)" ON)

if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
endif()

add_executable(oblivious_routing ${SRC_FILES})

if(USE_OPENMP)
    target_link_libraries(oblivious_routing PRIVATE OpenMP::OpenMP_CXX)
endif()

# Optional: keep a small macOS tweak if you still need it locally
if(APPLE AND USE_OPENMP)
    target_compile_options(oblivious_routing PRIVATE -Xclang -fopenmp)
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

target_link_libraries(oblivious_routing PRIVATE Eigen3::Eigen)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ---- Tests ----
if(BUILD_TESTS)
    add_subdirectory(test)
endif()
