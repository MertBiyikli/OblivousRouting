//
// Created by Mert Biyikli on 08.06.26.
//

#ifndef OBLIVIOUSROUTING_TEST_ELECTRICAL_MWU_HELPER_H
#define OBLIVIOUSROUTING_TEST_ELECTRICAL_MWU_HELPER_H

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <vector>
#include <Eigen/Dense>
#include "data_structures/graph/graph_csr.h"
#include "utils/config.h"
#include "algorithms/mwu/oracle/electrical/laplacian_solver.h"
#include "io/parse_argurment_io.h"

using Catch::Approx;
using Eigen::VectorXd;

inline std::unique_ptr<GraphCSR> makeTriangleGraph()
{
    auto g = std::make_unique<GraphCSR>(3);

    g->addEdge(0, 1, 10.0);
    g->addEdge(1, 2, 10.0);
    g->addEdge(0, 2, 5.0);

    g->finalize();
    return g;
}

inline std::unique_ptr<GraphCSR> makePathGraph()
{
    auto g = std::make_unique<GraphCSR>(4);

    g->addEdge(0, 1, 10.0);
    g->addEdge(1, 2, 10.0);
    g->addEdge(2, 3, 10.0);

    g->finalize();
    return g;
}

inline Config makeElectricalConfig()
{
    Config cfg;
    cfg.solvers = {SolverType::ELECTRICAL_SKETCHING};
    cfg.graph_format = GraphFormat::CSR;

    return cfg;
}

template <typename ResultPtr>
void requireValidElectricalResult(const ResultPtr& result)
{
    REQUIRE(result);

    REQUIRE(std::isfinite(result->oblivious_ratio));
    REQUIRE(result->oblivious_ratio >= 1.0);

    // Integration sanity bound.
    // This should be loose enough to avoid numerical brittleness,
    // but tight enough to catch totally broken behavior.
    REQUIRE(result->oblivious_ratio < 100.0);
}

#endif //OBLIVIOUSROUTING_TEST_ELECTRICAL_MWU_HELPER_H