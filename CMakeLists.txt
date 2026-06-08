cmake_minimum_required(VERSION 3.20)

project(ObliviousRouting
        VERSION 0.0.1
        DESCRIPTION "A comparison of tree metrics and electrical flows as oblivious routing oracles"
        LANGUAGES CXX)

# ============================================================================
# Options
# ============================================================================
option(BUILD_TESTS "Build unit and integration tests" OFF)
option(OR_ENABLE_WARNINGS "Enable strict compiler warnings" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_FIND_PACKAGE_PREFER_CONFIG ON)

if(DEFINED ENV{ORTOOLS_ROOT})
    list(PREPEND CMAKE_PREFIX_PATH "$ENV{ORTOOLS_ROOT}")
endif()

# ============================================================================
# Dependencies
# ============================================================================
find_package(Protobuf CONFIG REQUIRED)
find_package(ortools REQUIRED CONFIG)
find_package(Boost REQUIRED COMPONENTS program_options serialization)
find_package(Eigen3 REQUIRED)

# ============================================================================
# Core library: all production code except app/main.cpp
# ============================================================================
add_library(oblivious_routing_core STATIC
        include/core/routing_table.h
        source/core/routing_table.cpp
        include/core/routing_scheme.h
        source/core/routing_scheme.cpp

        include/io/graph_io.h
        source/io/graph_io.cpp
        include/io/solver_io.h
        include/io/parse_argurment_io.h
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

        include/algorithms/mwu/electrical_mwu.h
        source/algorithms/mwu/electrical_mwu.cpp
        include/algorithms/mwu/tree_mwu.h
        source/algorithms/mwu/tree_mwu.cpp

        include/algorithms/mwu/oracle/electrical/laplacian_solver.h
        source/algorithms/mwu/oracle/electrical/laplacian_solver.cpp
        include/algorithms/mwu/oracle/electrical/graph_to_laplacian.h
        source/algorithms/mwu/oracle/electrical/graph_to_laplacian.cpp

        include/algorithms/mwu/oracle/tree/tree_transform.h
        include/algorithms/mwu/oracle/tree/tree_oracle.h
        include/algorithms/mwu/oracle/tree/mst/mst_oracle.h
        include/algorithms/mwu/oracle/tree/mst/mst_algo.h
        source/algorithms/mwu/oracle/tree/mst/mst_algo.cpp
        include/algorithms/mwu/oracle/tree/frt/frt.h
        include/algorithms/mwu/oracle/tree/fast_ckr/fast_ckr.h

        include/algorithms/lp/lp_base.h
        include/algorithms/lp/lp_ac.h
        include/algorithms/lp/lp_mcf.h
        source/algorithms/lp/lp_base.cpp
        source/algorithms/lp/lp_ac.cpp
        source/algorithms/lp/lp_mcf.cpp

        include/utils/time_tracking.h
        include/utils/demands.h
        source/utils/demands.cpp
        include/utils/hash.h
        include/utils/config.h
)

target_include_directories(oblivious_routing_core
        PUBLIC
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/source
        ${CMAKE_SOURCE_DIR}/experiments
        ${CMAKE_SOURCE_DIR}/extern/amgcl
)

target_link_libraries(oblivious_routing_core
        PUBLIC
        ortools::ortools
        protobuf::libprotobuf
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
            unit_tests PRIVATE tests/unit/graph/ tests/unit/tree/
    )

    target_link_libraries(unit_tests
            PRIVATE
            oblivious_routing_core
            Catch2::Catch2WithMain
    )

    catch_discover_tests(unit_tests
            TEST_PREFIX "unit."
            REPORTER console
    )

    # ------------------------------------------------------------------------
    # Integration tests
    # ------------------------------------------------------------------------
    add_executable(integration_tests
            tests/integration/electrical/test_electrical_mwu.cpp
            tests/integration/electrical/test_electrical_mwu_io.cpp
    )

    target_link_libraries(integration_tests
            PRIVATE
            oblivious_routing_core
            Catch2::Catch2WithMain
    )

    target_include_directories(
            integration_tests PRIVATE tests/integration/electrical/
    )

    target_compile_definitions(integration_tests
            PRIVATE
            PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
            PROJECT_BINARY_DIR="${CMAKE_BINARY_DIR}"
    )

    catch_discover_tests(integration_tests
            TEST_PREFIX "integration."
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