# Test targets are intentionally kept out of the top-level CMakeLists.txt.

include(Catch)

add_executable(unit_tests
    common/utils.h
    unit/graph/test_graph_common.cpp
    unit/graph/test_graph_csr.cpp
    unit/graph/test_graph_adj.cpp
    unit/tree/test_tree_oracle_common.cpp
    unit/tree/test_fast_ckr.cpp
    unit/tree/test_frt.cpp
    unit/tree/test_mst.cpp
    unit/test_electrical_flow.cpp
    unit/test_routing_table.cpp
)

target_include_directories(unit_tests PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/unit/graph"
    "${CMAKE_CURRENT_SOURCE_DIR}/unit/tree"
    "${CMAKE_CURRENT_SOURCE_DIR}/common"
    "${PROJECT_SOURCE_DIR}/include"
)

target_link_libraries(unit_tests PRIVATE
    oblivious_routing_core
    Catch2::Catch2WithMain
)

target_compile_definitions(unit_tests PRIVATE
    PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
    PROJECT_BINARY_DIR="${PROJECT_BINARY_DIR}"
)

or_apply_project_options(unit_tests)

catch_discover_tests(unit_tests
    TEST_PREFIX "unit."
    REPORTER console
)

add_executable(integration_tests
    integration/electrical/test_electrical_mwu.cpp
    integration/electrical/test_electrical_mwu_io.cpp
    integration/tree/test_tree_oracle.cpp
    integration/tree/test_tree_transform.cpp
    integration/tree/test_mst.cpp
    integration/tree/test_tree_mwu.cpp
    integration/lp/test_lp.cpp
    integration/test_routing_engine.cpp
)

target_include_directories(integration_tests PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/common"
    "${PROJECT_SOURCE_DIR}/include"
)

target_link_libraries(integration_tests PRIVATE
    oblivious_routing_core
    Catch2::Catch2WithMain
)

target_compile_definitions(integration_tests PRIVATE
    PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
    PROJECT_BINARY_DIR="${PROJECT_BINARY_DIR}"
)

or_apply_project_options(integration_tests)

catch_discover_tests(integration_tests
    TEST_PREFIX "integration."
    PROPERTIES LABELS "integration"
    REPORTER console
)

add_custom_target(run_tests
    COMMAND "${CMAKE_CTEST_COMMAND}" --output-on-failure
    DEPENDS unit_tests integration_tests
    WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
)
