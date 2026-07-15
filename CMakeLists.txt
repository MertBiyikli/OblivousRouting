cmake_minimum_required(VERSION 3.20)

project(ObliviousRouting
        VERSION 0.0.1
        DESCRIPTION "A comparison of tree metrics and electrical flows as oblivious routing oracles"
        LANGUAGES CXX)

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
option(BUILD_TESTS "Build unit and integration tests" OFF)
option(OR_ENABLE_WARNINGS "Enable strict compiler warnings" OFF)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ============================================================================
# Dependencies — Homebrew setup
# ============================================================================
set(CMAKE_FIND_PACKAGE_PREFER_CONFIG ON)

# Homebrew is at /usr/local on this Intel Mac.
set(PKG_CONFIG_EXECUTABLE "/usr/local/bin/pkg-config"
        CACHE FILEPATH "Homebrew pkg-config executable" FORCE)

set(ENV{PKG_CONFIG_PATH}
        "/usr/local/lib/pkgconfig:/usr/local/share/pkgconfig:$ENV{PKG_CONFIG_PATH}")

find_package(PkgConfig REQUIRED)
pkg_check_modules(Cbc       REQUIRED IMPORTED_TARGET cbc)
pkg_check_modules(OsiCbc    REQUIRED IMPORTED_TARGET osi-cbc)
pkg_check_modules(Clp       REQUIRED IMPORTED_TARGET clp)
pkg_check_modules(OsiClp    REQUIRED IMPORTED_TARGET osi-clp)
pkg_check_modules(Cgl       REQUIRED IMPORTED_TARGET cgl)
pkg_check_modules(Osi       REQUIRED IMPORTED_TARGET osi)
pkg_check_modules(CoinUtils REQUIRED IMPORTED_TARGET coinutils)
# This is the imported target Homebrew's OR-Tools export requires.
pkg_check_modules(Cbc REQUIRED IMPORTED_TARGET cbc)

find_package(ortools CONFIG REQUIRED)
get_filename_component(ORTOOLS_PREFIX "${ortools_DIR}/../../.." REALPATH)
set(ORTOOLS_INCLUDE_DIR "${ORTOOLS_PREFIX}/include")

if(NOT EXISTS "${ORTOOLS_INCLUDE_DIR}/ortools/linear_solver/linear_solver.h")
    message(FATAL_ERROR "OR-Tools headers were not found at: ${ORTOOLS_INCLUDE_DIR}")
endif()
find_package(Boost REQUIRED COMPONENTS program_options serialization)
find_package(Eigen3 REQUIRED)


# Never let CMake select Anaconda's fmt installation.
find_package(fmt CONFIG REQUIRED
        PATHS "/usr/local/lib/cmake/fmt"
        NO_DEFAULT_PATH)

# xcut fetches spdlog. Make that spdlog use the same external fmt.
set(SPDLOG_FMT_EXTERNAL ON CACHE BOOL
        "Use external fmt for spdlog" FORCE)

add_subdirectory(${CMAKE_SOURCE_DIR}/extern/xcut EXCLUDE_FROM_ALL)

target_link_libraries(xcut_lib PUBLIC fmt::fmt)

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
)

target_include_directories(oblivious_routing_core
        PUBLIC
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/source
        ${CMAKE_SOURCE_DIR}/experiments
        ${CMAKE_SOURCE_DIR}/extern/amgcl
        ${CMAKE_SOURCE_DIR}/extern/xcut/include
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

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ============================================================================
# Testing
# ============================================================================
if(BUILD_TESTS)
    enable_testing()

    find_package(Catch2 3 QUIET)

    if(NOT Catch2_FOUND)
        message(STATUS "Catch2 not found, fetching from GitHub...")
        include(FetchContent)
        FetchContent_Declare(
                Catch2
                GIT_REPOSITORY https://github.com/catchorg/Catch2.git
                GIT_TAG v3.5.0
        )
        FetchContent_MakeAvailable(Catch2)
    endif()

    include(CTest)
    include(Catch)

    # ------------------------------------------------------------------------
    # Unit tests
    # ------------------------------------------------------------------------
    add_executable(unit_tests
            tests/common/utils.h
            # Graph data structures
            tests/unit/graph/test_graph_common.cpp
            tests/unit/graph/test_graph_csr.cpp
            tests/unit/graph/test_graph_adj.cpp
            # Tree oracles
            tests/unit/tree/test_tree_oracle_common.cpp
            tests/unit/tree/test_fast_ckr.cpp
            tests/unit/tree/test_frt.cpp
            tests/unit/tree/test_mst.cpp
            tests/unit/test_electrical_flow.cpp
            tests/unit/test_routing_table.cpp
    )

    target_include_directories(
            unit_tests PRIVATE tests/unit/graph/ tests/unit/tree/ tests/common/ include/
    )

    target_link_libraries(unit_tests
            PRIVATE
            oblivious_routing_core
            Catch2::Catch2WithMain
    )

    target_compile_definitions(unit_tests PRIVATE
            PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
            PROJECT_BINARY_DIR="${PROJECT_BINARY_DIR}"
    )

    catch_discover_tests(unit_tests
            TEST_PREFIX "unit."
            REPORTER console
    )


    # ------------------------------------------------------------------------
    # Integration tests
    # ------------------------------------------------------------------------
    add_executable(integration_tests
            tests/common/utils.h
            tests/integration/electrical/test_electrical_mwu.cpp
            tests/integration/electrical/test_electrical_mwu_io.cpp

            tests/integration/tree/test_tree_oracle.cpp
            tests/integration/tree/test_tree_transform.cpp
            tests/integration/tree/test_mst.cpp
            tests/integration/tree/test_tree_mwu.cpp

            tests/integration/lp/test_lp.cpp
            tests/integration/test_routing_engine.cpp
    )

    target_link_libraries(integration_tests
            PRIVATE
            oblivious_routing_core
            Catch2::Catch2WithMain
    )

    target_include_directories(
            integration_tests PRIVATE tests/common/ include/
    )

    target_compile_definitions(integration_tests
            PRIVATE
            PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
            PROJECT_BINARY_DIR="${CMAKE_BINARY_DIR}"
    )

    catch_discover_tests(integration_tests
            TEST_PREFIX "integration."
            PROPERTIES LABELS "integration"
            REPORTER console
    )

    # ------------------------------------------------------------------------
    # Convenience target
    # ------------------------------------------------------------------------
    add_custom_target(run_tests
            COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
            DEPENDS unit_tests integration_tests
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )

    message(STATUS "Tests enabled.")
    message(STATUS "Build with: cmake --build .")
    message(STATUS "Run with:   cmake --build . --target run_tests")
else()
    message(STATUS "Tests disabled. Enable with: -DBUILD_TESTS=ON")
endif()