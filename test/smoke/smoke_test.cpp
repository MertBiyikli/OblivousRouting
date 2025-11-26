//
// Created by Mert Biyikli on 11.05.25.
//

#include "gtest/gtest.h"
#include "../../src/datastructures/graph.h"
#include "../../src/solver/solver.h"
#include "../../src/electrical/electrical_flow_optimized.h"
#include <stdexcept>

#include "../../../../../opt/anaconda3/include/gtest/gtest.h"
#include "../../src/parse_parameter.h"


// * --- Testing solvers --- *
TEST(SmokeTest, ElectricalFlowOptimizedBasicTest) {
    Graph g(3);
    g.addEdge(0, 1, 1.0);
    g.addEdge(1, 2, 1.0);
    g.addEdge(0, 2, 1.0);

    std::unique_ptr<ObliviousRoutingSolver> &solver;
    solver->init(g);

    std::vector<std::pair<int, int>> demands = { {0, 2} };
    std::vector<double> demand_values = { 1.0 };

    solver->runSolve(g);

    int argc = 3;
    char* argv[] = { (char*)"electrical_optimized", (char*)"input_graph.txt", (char*)"gravity" };
    std::string err;
    auto cfg = parse_parameter(argc, argv, &err);
    if (!cfg) { std::cerr << err; return -1; } // EX_USAGE


    auto [offline, solver_cong] = HandleDemandModel(argc, argv, cfg, g, solver);

    EXPECT_TRUE((solver_cong >= offline));
}