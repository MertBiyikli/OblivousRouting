//
// Created by Mert Biyikli on 09.06.26.
//
#include <catch2/catch_test_macros.hpp>

#include <memory>

#include "../../common/utils.h"

#include "algorithms/mwu/tree_mwu.h"
#include "algorithms/mwu/oracle/tree/frt/frt.h"
#include "algorithms/mwu/oracle/tree/fast_ckr/fast_ckr.h"
#include "algorithms/mwu/oracle/tree/mst/mst_oracle.h"

#include "data_structures/hst/flat_hst.h"
#include "core/routing_table.h"

using namespace integration;

TEST_CASE("TreeMWU runs end-to-end with FRT flat HST oracle",
          "[integration][tree][mwu][frt][flat-hst]")
{
    auto graph = makeCycleGraph();

    auto oracle = std::make_unique<FRT<FlatHST>>(*graph);
    TreeMWU<FlatHST> solver(*graph, 0, std::move(oracle));

    LinearRoutingTable table;
    solver.computeBasisFlows(table);

    requireNonEmptyLinearTable(table, *graph);
    REQUIRE(solver.getIterationCount() > 0);
}

TEST_CASE("TreeMWU runs end-to-end with FastCKR flat HST oracle",
          "[integration][tree][mwu][ckr][flat-hst]")
{
    auto graph = makeCycleGraph();

    auto oracle = std::make_unique<FastCKR<FlatHST>>(*graph);
    TreeMWU<FlatHST> solver(*graph, 0, std::move(oracle));

    LinearRoutingTable table;
    solver.computeBasisFlows(table);

    requireNonEmptyLinearTable(table, *graph);
    REQUIRE(solver.getIterationCount() > 0);
}

TEST_CASE("TreeMWU runs end-to-end with MST flat HST oracle",
          "[integration][tree][mwu][mst][flat-hst]")
{
    auto graph = makeCycleGraph();

    auto oracle = std::make_unique<TreeMST<FlatHST>>(*graph);
    TreeMWU<FlatHST> solver(*graph, 0, std::move(oracle));

    LinearRoutingTable table;
    solver.computeBasisFlows(table);

    requireNonEmptyLinearTable(table, *graph);
    REQUIRE(solver.getIterationCount() > 0);
}