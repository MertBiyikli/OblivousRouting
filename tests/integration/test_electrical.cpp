//
// Created by Mert Biyikli on 05.06.26.
//
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <vector>
#include <Eigen/Dense>
#include "data_structures/graph/graph_csr.h"
#include "utils/config.h"
#include "algorithms/mwu/oracle/electrical/laplacian_solver.h"

using Catch::Approx;
using Eigen::VectorXd;

TEST_CASE("Electrical Flow - Simple Oblivious Routing solve", "[ElectricalFlowMWU]")
{
    auto g = std::make_unique<GraphCSR>(3);
    g->addEdge(0, 1, 10.0);
    g->addEdge(1, 2, 10.0);
    g->addEdge(0, 2, 5.0);
    g->finalize();

    Config cfg;
    cfg.solvers = {SolverType::ELECTRICAL_SKETCHING};
    cfg.graph_format = GraphFormat::CSR;

    RoutingEngine engine;
    auto result = engine.solve(*g, cfg, cfg.solvers[0]);

    REQUIRE(result);
    double obl_ratio = result->oblivious_ratio;
    REQUIRE(obl_ratio >= 1.0); // Oblivious ratio should be at least 1
}