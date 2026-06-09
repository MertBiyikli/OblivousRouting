//
// Created by Mert Biyikli on 09.06.26.
//

#include <catch2/catch_test_macros.hpp>

#include <set>
#include <utility>
#include <vector>

#include "../common/utils.h"

#include "algorithms/mwu/oracle/tree/mst/mst_algo.h"

using namespace integration;

TEST_CASE("MST algorithm returns n-1 edges on a connected graph",
          "[integration][tree][mst]")
{
    auto graph = makeCycleGraph();

    MST mst(*graph);
    auto edges = mst.computeMST();

    REQUIRE(static_cast<int>(edges.size()) == graph->getNumNodes() - 1);

    std::set<std::pair<int, int>> unique_edges;
    for (auto [u, v] : edges) {
        REQUIRE(u >= 0);
        REQUIRE(v >= 0);
        REQUIRE(u < graph->getNumNodes());
        REQUIRE(v < graph->getNumNodes());
        REQUIRE(u != v);

        REQUIRE(graph->getEdgeId(u, v) != INVALID_EDGE_ID);

        unique_edges.insert(normalizedEdge(u, v));
    }

    REQUIRE(unique_edges.size() == edges.size());
}