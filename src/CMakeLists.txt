file(GLOB SRC_FILES *.cpp)
add_library(core ${SRC_FILES}
        solver/solver.h
        utils/hash.h
        tree_based/frt/mcct/Node.h
        tree_based/frt/mcct/mcct_derandomized_weighted_solver.h
        tree_based/frt/raecke_frt.cpp
        tree_based/frt/mcct/mcct_derandomized_weighted_solver.h
        tree_based/frt/mcct/mcct_derandomized_weighted_solver.cpp
        lp_solver/MCCF_lp_solver.cpp
        lp_solver/MCCF_lp_solver.h
        electrical/laplacian_solvers/AMGSolver.cpp
        electrical/laplacian_solvers/AMGSolver.h
)
target_include_directories(core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
