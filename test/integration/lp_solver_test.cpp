//
// Created by Mert Biyikli on 13.05.25.
//



#include "gtest/gtest.h"
#include "../../src/lp_solver/LPSolver.h"
#include <stdexcept>

TEST(LPSolver_TriangleGraph, NodeAndEdgeExistence) {
    LPSolver solver;

    Graph g(3);
    g.addEdge(0, 1, 1.0);
    g.addEdge(1, 2, 1.0);
    g.addEdge(0, 2, 1.0);

    solver.solve(g);
    //solver.PrintObliviousRoutingSolution(g.GetDiGraph());

}
