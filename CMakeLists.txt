cmake_minimum_required(VERSION 3.20)
project(ObliviousRouting CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(Boost_NO_BOOST_CMAKE ON)  # prefer system/Homebrew Boost over any packaged config

# only require what the main binary needs
find_package(Boost REQUIRED COMPONENTS program_options serialization)

# if you have tests that need unit_test_framework, guard it:
option(BUILD_TESTS "Build tests" ON)
if(BUILD_TESTS)
    find_package(Boost 1.88.0 REQUIRED COMPONENTS unit_test_framework)
endif()



# ---- amgcl (vendored) ----
add_subdirectory(extern/amgcl)



find_package(ortools CONFIG QUIET)


include(FetchContent)
if(NOT ortools_FOUND)
    message(STATUS "OR-Tools not found. Fetching from GitHub...")

    FetchContent_Declare(
            ortools
            GIT_REPOSITORY https://github.com/google/or-tools.git
            GIT_TAG        v9.9 # you can update this to the latest stable release
    )

    set(BUILD_DEPS ON)            # OR-Tools builds protobuf/abseil/etc. itself
    set(BUILD_SAMPLES OFF)
    set(BUILD_EXAMPLES OFF)
    set(BUILD_TESTING OFF)
    set(BUILD_PYTHON OFF)
    set(BUILD_JAVA OFF)
    set(BUILD_CSHARP OFF)
    # Disable heavy/paid solvers you don't need
    set(BUILD_SCIP OFF)
    set(USE_GUROBI OFF)
    set(USE_XPRESS OFF)
    set(USE_HIGHS OFF)

    FetchContent_MakeAvailable(ortools)
endif()




# ---- Your target ----
add_executable(oblivious_routing
        main.cpp
        src/graph.cpp
        src/lp_solver/LPSolver.cpp
        src/lp_solver/MCCF_lp_solver.cpp
        src/tree_based/raecke_tree_decomp.cpp
        src/tree_based/mcct/Node.cpp
        src/tree_based/mcct/mcct_derandomized_weighted_solver.cpp
        src/tree_based/raecke_transform.cpp
        src/electrical/electrical_flow_naive.cpp
        src/electrical/AMGSolver.cpp
        # (add any other .cpp files you have)
)


# target_link_libraries(oblivious_routing ortools::ortools)
target_include_directories(oblivious_routing PRIVATE ${CMAKE_SOURCE_DIR}/extern/amgcl)
target_link_libraries(oblivious_routing
        PRIVATE
        ortools::ortools
        Boost::boost
        Boost::program_options
        Boost::serialization
)
