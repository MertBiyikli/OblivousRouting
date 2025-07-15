file(GLOB SRC_FILES *.cpp)
add_library(core ${SRC_FILES}
        solver/solver.h
        utils/hash.h
        tree_based/mcct/Node.h
        tree_based/mcct/Tree.h
        tree_based/mcct/mcct_derandomized_weighted_solver.h
        tree_based/raecke_tree_decomp.cpp
        tree_based/raecke_tree_decomp.h
        tree_based/mcct/Node.cpp
        tree_based/mcct/mcct_derandomized_weighted_solver.h
        tree_based/mcct/mcct_derandomized_weighted_solver.cpp
        tree_based/raecke_tree_decomp.h
        tree_based/raecke_tree_decomp.cpp
        tree_based/raecke_transform.cpp
        tree_based/raecke_transform.h
        lp_solver/MCCF_lp_solver.cpp
        lp_solver/MCCF_lp_solver.h
        electrical/electrical_flow_naive.cpp
        electrical/electrical_flow_naive.h
        electrical/laplcian_solver.cpp
        electrical/laplcian_solver.h
)
target_include_directories(core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
