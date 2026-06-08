#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "data_structures/graph/graph_csr.h"
#include "algorithms/mwu/oracle/tree/frt/frt.h"
#include <vector>

#include "tree_oracle_test_helpers.h"
#include "catch2/catch_template_test_macros.hpp"

using Catch::Approx;

// -----------------------------------------------------------------------------
// FRT-specific tests
//
// These tests should only describe behavior that is special to FRT.
// Generic TreeOracle behavior belongs in test_tree_oracle_common.cpp.
// -----------------------------------------------------------------------------

TEST_CASE("FRT stores the exact delta value as level radius", "[tree-oracle][frt][radius]") {
    auto graph = makeSimplePathGraph5();
    FRT<FlatHST> oracle(graph);

    const std::vector<int> permutation = identityPermutation(graph.getNumNodes());
    HSTLevel level;

    const double delta = 2.5;

    oracle.computeLevelPartition(graph, level, permutation, delta);

    REQUIRE(level.R == Approx(delta));
}

TEST_CASE("FRT with a large radius and one center assigns all nodes to that center", "[tree-oracle][frt][single-center]") {
    auto graph = makeSimplePathGraph5();
    FRT<FlatHST> oracle(graph);

    const std::vector<int> permutation = {0};
    HSTLevel level;

    oracle.computeLevelPartition(graph, level, permutation, 1000.0);

    REQUIRE(level.centers.size() == 1);
    REQUIRE(level.centers[0] == 0);

    for (int owner : level.owner) {
        REQUIRE(owner == 0);
    }
}

TEST_CASE("FRT creates multiple centers for a small radius on a path graph", "[tree-oracle][frt][multi-center]") {
    auto graph = makePathGraph(6);
    FRT<FlatHST> oracle(graph);

    const std::vector<int> permutation = identityPermutation(graph.getNumNodes());
    HSTLevel level;

    oracle.computeLevelPartition(graph, level, permutation, 0.5);

    REQUIRE(level.centers.size() > 1);
    requirePartitionShapeIsValid(level, graph.getNumNodes(), permutation);
    requireAssignedOwnersAreCenters(level);
}

TEST_CASE("FRT assigns all nodes when the permutation contains all nodes", "[tree-oracle][frt][owners]") {
    auto graph = makeSimplePathGraph5();
    FRT<FlatHST> oracle(graph);

    const std::vector<int> permutation = identityPermutation(graph.getNumNodes());
    HSTLevel level;

    oracle.computeLevelPartition(graph, level, permutation, 5.0);

    requireAllOwnersAreValidNodes(level, graph.getNumNodes());
    requireAssignedOwnersAreCenters(level);
}

TEST_CASE("FRT radius is stable across different delta values", "[tree-oracle][frt][delta]") {
    auto graph = makeSimplePathGraph5();
    FRT<FlatHST> oracle(graph);

    const std::vector<int> permutation = identityPermutation(graph.getNumNodes());

    for (double delta : {0.5, 1.0, 2.0, 5.0}) {
        HSTLevel level;

        oracle.computeLevelPartition(graph, level, permutation, delta);

        REQUIRE(level.R == Approx(delta));
    }
}