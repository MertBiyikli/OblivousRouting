//
// Created by Mert Biyikli on 08.06.26.
//
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <vector>

#include "tree_oracle_test_helpers.h"

using Catch::Approx;

// -----------------------------------------------------------------------------
// TreeMST-specific tests
//
// TreeMST is different from FRT and FastCKR:
//   - it does not use computeLevelPartition(...)
//   - it overrides getTree(...)
//   - it builds an MST-based decomposition tree directly
// -----------------------------------------------------------------------------

TEST_CASE("TreeMST can be constructed from a graph", "[tree-oracle][mst][construction]") {
    auto graph = makeSimplePathGraph5();

    TreeMST<FlatHST> oracle(graph);

    REQUIRE(oracle.graph.getNumNodes() == 5);
}

TEST_CASE("TreeMST builds a tree on a simple path graph", "[tree-oracle][mst][get-tree]") {
    auto graph = makeSimplePathGraph5();

    TreeMST<FlatHST> oracle(graph);

    auto distances = currentGraphDistances(graph);

    REQUIRE_NOTHROW(oracle.getTree(distances));
}

TEST_CASE("TreeMST builds a tree on a grid graph", "[tree-oracle][mst][get-tree][grid]") {
    auto graph = makeGridGraph(3);

    TreeMST<FlatHST> oracle(graph);

    auto distances = currentGraphDistances(graph);

    REQUIRE_NOTHROW(oracle.getTree(distances));
}

TEST_CASE("TreeMST builds a tree on a weighted cycle graph", "[tree-oracle][mst][get-tree][cycle]") {
    auto graph = makeWeightedCycleGraph6();

    TreeMST<FlatHST> oracle(graph);

    auto distances = currentGraphDistances(graph);

    REQUIRE_NOTHROW(oracle.getTree(distances));
}

TEST_CASE("TreeMST updates graph distances before computing the MST tree", "[tree-oracle][mst][distance-update]") {
    auto graph = makeSimplePathGraph5();

    TreeMST<FlatHST> oracle(graph);

    auto distances = increasingDistances(graph);

    REQUIRE_NOTHROW(oracle.getTree(distances));

    for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
        REQUIRE(graph.getEdgeDistance(e) == Approx(distances[e]));
    }
}

TEST_CASE("TreeMST records at least one scale after building the tree", "[tree-oracle][mst][scales]") {
    auto graph = makeSimplePathGraph5();

    TreeMST<FlatHST> oracle(graph);

    auto distances = currentGraphDistances(graph);

    REQUIRE(oracle.scales.empty());

    REQUIRE_NOTHROW(oracle.getTree(distances));

    REQUIRE(!oracle.scales.empty());
}

TEST_CASE("TreeMST can rebuild after distances change", "[tree-oracle][mst][rebuild]") {
    auto graph = makeWeightedCycleGraph6();

    TreeMST<FlatHST> oracle(graph);

    auto first_distances = currentGraphDistances(graph);
    auto second_distances = increasingDistances(graph);

    REQUIRE_NOTHROW(oracle.getTree(first_distances));
    REQUIRE_NOTHROW(oracle.getTree(second_distances));

    for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
        REQUIRE(graph.getEdgeDistance(e) == Approx(second_distances[e]));
    }

    REQUIRE(!oracle.scales.empty());
}

TEST_CASE("TreeMST computeLevelPartition is intentionally unused", "[tree-oracle][mst][level-partition]") {
    auto graph = makeSimplePathGraph5();

    TreeMST<FlatHST> oracle(graph);

    const auto permutation = identityPermutation(graph.getNumNodes());
    HSTLevel level;

    oracle.computeLevelPartition(graph, level, permutation, 2.0);

    REQUIRE(level.centers.empty());
    REQUIRE(level.owner.empty());
}