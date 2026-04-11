cmake_minimum_required(VERSION 3.20)
project(ObliviousRouting
        VERSION 0.0.1
        DESCRIPTION "A comparison of tree metrics and electrical flows as oblivious routing oracles"
        LANGUAGES CXX)

# --- Options ---
option(BUILD_TESTS "Build test targets" OFF)
option(OR_ENABLE_WARNINGS "Enable strict compiler warnings" OFF)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- Dependencies ---
set(CMAKE_FIND_PACKAGE_PREFER_CONFIG ON)

# Point CMake to your protobuf 6.31.1 package config
set(Protobuf_DIR "$ENV{HOME}/local/protobuf-6.31.1/lib/cmake/protobuf" CACHE PATH "")

find_package(Protobuf CONFIG REQUIRED)
find_package(ortools REQUIRED CONFIG)
find_package(Boost REQUIRED COMPONENTS program_options serialization)
find_package(ortools REQUIRED CONFIG)
find_package(Eigen3 REQUIRED)


# --- Target ---
add_executable(oblivious_routing
        app/main.cpp

        include/core/routing_table.h
        source/core/routing_table.cpp
        include/core/routing_scheme.h
        source/core/routing_scheme.cpp
        include/io/graph_io.h
        source/io/graph_io.cpp

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
        include/data_structures/union_find/union_find.h
        source/data_structures/union_find/union_find.cpp
        include/io/solver_io.h
        include/io/parse_argurment_io.h
        include/io/demand_io.h
        source/data_structures/hst/flat_hst.cpp
        source/data_structures/hst/pointer_hst.cpp
        include/algorithms/mwu/oracle/tree/fast_ckr/engineering_fast_ckr.h
        include/core/routing_table.h
)


target_include_directories(oblivious_routing BEFORE PRIVATE
        "$ENV{HOME}/local/protobuf-6.31.1/include"
)

target_include_directories(oblivious_routing
        PRIVATE
        ${CMAKE_SOURCE_DIR}/source
        ${CMAKE_SOURCE_DIR}/experiments
        ${CMAKE_SOURCE_DIR}/extern/amgcl
)

# Link libs
target_link_libraries(oblivious_routing PRIVATE
        ortools::ortools
        protobuf::libprotobuf
        Eigen3::Eigen
        Boost::program_options
        Boost::serialization
)


message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

