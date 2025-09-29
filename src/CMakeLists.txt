file(GLOB SRC_FILES *.cpp)
add_library(core ${SRC_FILES}
        solver/solver.h
        utils/hash.h
        tree_based/frt/mcct/Node.h
        tree_based/frt/mcct/Tree.h
        tree_based/frt/mcct/mcct_derandomized_weighted_solver.h
        tree_based/frt/raecke_frt.cpp
        tree_based/frt/raecke_frt.h
        tree_based/frt/mcct/Node.cpp
        tree_based/frt/mcct/mcct_derandomized_weighted_solver.h
        tree_based/frt/mcct/mcct_derandomized_weighted_solver.cpp
        tree_based/frt/raecke_frt.h
        tree_based/frt/raecke_frt.cpp
        tree_based/frt/raecke_frt_transform.h
        lp_solver/MCCF_lp_solver.cpp
        lp_solver/MCCF_lp_solver.h
        electrical/electrical_flow_naive.cpp
        electrical/electrical_flow_naive.h
        electrical/laplacian_solvers/AMGSolver.cpp
        electrical/laplacian_solvers/AMGSolver.h
)
target_include_directories(core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
