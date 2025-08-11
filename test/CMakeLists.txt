add_executable(smoke_test smoke/smoke_test.cpp
        ../src/utils/hash.h)
target_link_libraries(smoke_test PRIVATE core gtest_main)

# Add correct include paths
target_include_directories(smoke_test PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${gtest_SOURCE_DIR}/include
)

add_test(NAME SmokeTest COMMAND smoke_test)

add_executable(integration_test
        integration/lp_solver_test.cpp
        integration/mcct_node.cpp
        integration/raecke_test.cpp
        ../src/utils/hash.h
        ../src/tree_based/raecke_transform.cpp
        ../src/tree_based/raecke_transform.h
        ../src/lp_solver/LPSolver.cpp
        ../src/lp_solver/MCCF_lp_solver.h
        ../src/lp_solver/MCCF_lp_solver.cpp
)
target_link_libraries(integration_test PRIVATE core gtest gtest_main)


target_include_directories(integration_test PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/integration
        ${gtest_SOURCE_DIR}/src/lp_solver

)

# MCCTs test
add_test(NAME mcctTest COMMAND mcct_test)

add_executable(mcct_test integration/mcct_node.cpp ../src/tree_based/mcct/Node.h ../src/tree_based/mcct/Tree.h ../src/tree_based/mcct/mcct_derandomized_weighted_solver.h
        ../src/tree_based/raecke_tree_decomp.cpp
        ../src/tree_based/raecke_tree_decomp.h
        ../src/tree_based/mcct/Node.cpp
        ../src/tree_based/mcct/mcct_derandomized_weighted_solver.cpp
        ../src/utils/hash.h)


target_link_libraries(mcct_test PRIVATE core gtest gtest_main)


target_include_directories(mcct_test PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/integration
        ${gtest_SOURCE_DIR}/src/lp_solver
        ${gtest_SOURCE_DIR}/src/tree_based

)


add_test(NAME IntegrationSolverTest COMMAND integration_test)

include(FetchContent)

# Try to find system-installed OR-Tools first
# find_package(ortools QUIET)

if(NOT ortools_FOUND)
    message(STATUS "OR-Tools not found. Fetching from GitHub...")

    FetchContent_Declare(
            ortools
            GIT_REPOSITORY https://github.com/google/or-tools.git
            GIT_TAG        v9.9 # you can update this to the latest stable release
    )

    set(BUILD_DEPS ON)
    set(BUILD_SAMPLES OFF)
    set(BUILD_TESTING OFF)
    set(BUILD_PYTHON OFF)
    set(BUILD_JAVA OFF)
    set(BUILD_CSHARP OFF)


    FetchContent_MakeAvailable(ortools)
endif()

target_link_libraries(integration_test PRIVATE ortools::ortools)
