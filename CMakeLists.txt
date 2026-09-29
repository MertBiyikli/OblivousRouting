cmake_minimum_required(VERSION 3.25)

project(E-Routing VERSION 0.0.1 DESCRIPTION "E-Routing: Network Resilience and Capacity Planning Engine." LANGUAGES CXX)

# Set CMake policies early to suppress warnings
if(POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW)
endif()
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

# ============================================================================
# Options
# ============================================================================
option(BUILD_TESTING "Build unit and integration tests" OFF)
option(OR_ENABLE_WARNINGS "Enable strict compiler warnings" OFF)
option(OR_ENABLE_OPENMP "Enable OpenMP support for parallel algorithms" ON)
option(OR_ENABLE_ASAN "Enable AddressSanitizer" OFF)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# ============================================================================
# Build helpers and dependencies
# ============================================================================
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

include(ProjectOptions)
include(Dependencies)

# ============================================================================
# Core library: all production code except app/main.cpp
# ============================================================================
add_library(oblivious_routing_core STATIC
        include/routing/routing_table.h
        include/routing/routing_engine.h

        include/io/graph_io.h
        include/io/solver_io.h
        include/io/parse_argument_io.h
        include/io/demand_io.h

        include/data_structures/graph/Igraph.h
        include/data_structures/graph/graph_adj.h
        include/data_structures/graph/graph_csr.h
        source/data_structures/graph/graph_adj.cpp
        source/data_structures/graph/graph_csr.cpp

        include/data_structures/mendel_scaling/quotient_graph.h
        include/data_structures/mendel_scaling/ultrametric_tree.h
        source/data_structures/mendel_scaling/quotient_graph.cpp
        source/data_structures/mendel_scaling/ultrametric_tree.cpp

        include/data_structures/hst/pointer_hst.h
        include/data_structures/hst/flat_hst.h
        source/data_structures/hst/pointer_hst.cpp
        source/data_structures/hst/flat_hst.cpp

        include/data_structures/union_find/union_find.h
        source/data_structures/union_find/union_find.cpp


        include/algorithms/oblivious/mwu/flow_sparsifier_mwu.h

        include/algorithms/oblivious/mwu/electrical_mwu.h
        source/algorithms/oblivious/mwu/electrical_mwu.cpp
        include/algorithms/oblivious/mwu/tree_mwu.h
        source/algorithms/oblivious/mwu/tree_mwu.cpp

        include/algorithms/oblivious/mwu/oracle/electrical/laplacian_solver.h
        source/algorithms/oblivious/mwu/oracle/electrical/laplacian_solver.cpp
        include/algorithms/oblivious/mwu/oracle/electrical/graph_to_laplacian.h
        source/algorithms/oblivious/mwu/oracle/electrical/graph_to_laplacian.cpp

        include/algorithms/oblivious/mwu/oracle/tree/tree_transform.h
        include/algorithms/oblivious/mwu/oracle/tree/tree_oracle.h
        include/algorithms/oblivious/mwu/oracle/tree/mst/mst_oracle.h
        include/algorithms/oblivious/mwu/oracle/tree/mst/mst_algo.h
        source/algorithms/oblivious/mwu/oracle/tree/mst/mst_algo.cpp
        include/algorithms/oblivious/mwu/oracle/tree/frt/frt.h
        include/algorithms/oblivious/mwu/oracle/tree/fast_ckr/fast_ckr.h

        include/algorithms/lp/lp_base.h
        include/algorithms/lp/lp_ac.h
        include/algorithms/lp/lp_mcf.h
        source/algorithms/lp/lp_ac.cpp
        source/algorithms/lp/lp_mcf.cpp

        include/utils/time_tracking.h
        include/utils/demands.h
        source/utils/demands.cpp
        include/utils/hash.h
        include/algorithms/semi_oblivious/semi_routing_engine.h
        include/algorithms/semi_oblivious/preprocessing/candidate_routing_scheme.h
        include/algorithms/semi_oblivious/postprocessing/load_optimizer.h
        include/algorithms/semi_oblivious/semi_oblivious_solver.h
        include/algorithms/semi_oblivious/postprocessing/or_tools_optimizer.h
        source/algorithms/semi_oblivious/postprocessing/or_tools_optimizer.cpp
        include/core/config.h
        include/routing/routing_result.h
        include/core/solver.h
        include/algorithms/lp/offline_solver.h
        include/core/types.h
        source/routing/routing_engine.cpp
        include/routing/storage/linear_routing_table.h
        include/routing/storage/allpair_routing_table.h
        source/routing/storage/linear_routing_table.cpp
        source/routing/storage/allpair_routing_table.cpp
        include/routing/routing_runner.h
        include/routing/routing_validation.h
        source/routing/routing_validation.cpp
        include/algorithms/semi_oblivious/semi_routing_runner.h
        include/algorithms/oblivious/oblivious_routing_runner.h
        source/algorithms/semi_oblivious/semi_routing_engine.cpp
        source/algorithms/semi_oblivious/semi_routing_runner.cpp
        source/algorithms/semi_oblivious/semi_oblivious_solver.cpp
        include/algorithms/semi_oblivious/expander_hierarchy/hierarchy_results.h
        include/algorithms/semi_oblivious/expander_hierarchy/preprocessing/hierarchy_preprocessor.h
        source/algorithms/oblivious/oblivious_routing_runner.cpp
        include/core/errors.h
        include/io/result_io.h
        include/algorithms/semi_oblivious/expander_hierarchy/routing_oracle/local_electrical_solver.h
        source/algorithms/semi_oblivious/expander_hierarchy/preprocessing/hierarchy_preprocessor.cpp
        source/algorithms/semi_oblivious/expander_hierarchy/routing_oracle/local_electrical_solver.cpp
        include/algorithms/semi_oblivious/expander_hierarchy/tree_sparsifier.h
        source/algorithms/semi_oblivious/expander_hierarchy/tree_sparsifier.cpp
        include/algorithms/semi_oblivious/expander_hierarchy/tree_flow_router.h
        source/algorithms/semi_oblivious/expander_hierarchy/tree_flow_router.cpp
        source/algorithms/semi_oblivious/expander_hierarchy/tree_flow_electrical_embedder.cpp
        include/algorithms/semi_oblivious/expander_hierarchy/tree_flow_electrical_embedder.h
        include/algorithms/semi_oblivious/expander_hierarchy/tree_sparsifier_solver.h
        # include/algorithms/oblivious/flow_sparsifier/boundary_linked_exp_decomp/expander_decomp.h
        # source/algorithms/oblivious/flow_sparsifier/boundary_linked_exp_decomp/expander_decomp.cpp
        # include/algorithms/oblivious/flow_sparsifier/boundary_linked_exp_decomp/augmented_cluster_view.h
        # source/algorithms/oblivious/flow_sparsifier/boundary_linked_exp_decomp/augmented_cluster_view.cpp
        # include/algorithms/oblivious/flow_sparsifier/boundary_linked_exp_decomp/practical_adapter.h
        # source/algorithms/oblivious/flow_sparsifier/boundary_linked_exp_decomp/practical_adapter.cpp
        # source/algorithms/oblivious/flow_sparsifier/boundary_linked_exp_decomp/practical_logger.cpp
        # include/algorithms/oblivious/flow_sparsifier/boundary_linked_exp_decomp/validation.h
        # source/algorithms/oblivious/flow_sparsifier/boundary_linked_exp_decomp/validation.cpp
        source/algorithms/oblivious/mwu/flow_sparsifier_mwu.cpp
        include/visualization/visualization_result.h
        source/visualization/visualization_result.cpp
        include/visualization/json_exporter.h
        source/visualization/json_exporter.cpp
        include/visualization/failure_analysis.h
        source/visualization/failure_analysis.cpp
        include/data_structures/graph/graph.h
)

target_include_directories(oblivious_routing_core
        PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/source
        ${PROJECT_SOURCE_DIR}/experiments
        ${PROJECT_SOURCE_DIR}/extern/amgcl
        ${PROJECT_SOURCE_DIR}/extern/xcut/include
        ${ORTOOLS_INCLUDE_DIR}
)

target_link_libraries(oblivious_routing_core
        PUBLIC
        xcut_lib
        ortools::ortools
        Eigen3::Eigen
        Boost::program_options
        Boost::serialization
)

if(OR_ENABLE_OPENMP)
    target_link_libraries(oblivious_routing_core PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(oblivious_routing_core PUBLIC OR_ENABLE_OPENMP)
endif()

target_compile_definitions(oblivious_routing_core PUBLIC EIGEN_DONT_PARALLELIZE)
# ============================================================================
# Main executable
# ============================================================================
add_executable(oblivious_routing
        app/main.cpp
)

target_link_libraries(oblivious_routing
        PRIVATE
        oblivious_routing_core
)

or_apply_project_options(oblivious_routing_core)
or_apply_project_options(oblivious_routing)

# See cmake/ProjectOptions.cmake: these translation units directly touch
# OR-Tools' MPSolver/MPObjective internals and must be compiled without ASan
# to remain ABI-compatible with the prebuilt (non-ASan) OR-Tools library.
or_exempt_sources_from_asan(
        source/algorithms/lp/lp_ac.cpp
        source/algorithms/lp/lp_mcf.cpp
        source/algorithms/semi_oblivious/postprocessing/or_tools_optimizer.cpp
)

# ============================================================================
# Testing
# ============================================================================
include(CTest)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

# ============================================================================
# Configuration summary
# ============================================================================
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "OpenMP enabled: ${OR_ENABLE_OPENMP}")
message(STATUS "Tests enabled: ${BUILD_TESTING}")
message(STATUS "Warnings enabled: ${OR_ENABLE_WARNINGS}")
message(STATUS "ASan enabled: ${OR_ENABLE_ASAN}")
