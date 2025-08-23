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
    find_package(Boost 1.70 REQUIRED COMPONENTS unit_test_framework)
endif()


# ---- amgcl (vendored) ----
add_subdirectory(extern/amgcl)

list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/local/ortools")
find_package(ortools REQUIRED CONFIG)






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
        src/tree_based/raecke_transform.h
        src/electrical/electrical_flow_naive.cpp
        src/electrical/AMGSolver.cpp
        experiments/performance/demands/GravityModel.cpp
        experiments/performance/demands/GravityModel.h
        experiments/performance/Performance.cpp
        experiments/performance/Performance.h
        experiments/performance/LP_Oblivious_Ratio.cpp
        experiments/performance/LP_Oblivious_Ratio.h
        src/tree_based/raecke_solver.cpp
        src/tree_based/raecke_solver.h
        src/parse_parameter.h

        # (add any other .cpp files you have)
)


# target_link_libraries(oblivious_routing ortools::ortools)
target_include_directories(oblivious_routing PRIVATE ${CMAKE_SOURCE_DIR}/extern/amgcl)
target_link_libraries(oblivious_routing PRIVATE ortools::ortools)
target_link_libraries(oblivious_routing
        PRIVATE
        Boost::boost
        Boost::program_options
        Boost::serialization
)
