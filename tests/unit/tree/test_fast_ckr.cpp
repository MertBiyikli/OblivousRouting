#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <vector>

#include "tree_oracle_test_helpers.h"

using Catch::Approx;

// -----------------------------------------------------------------------------
// FastCKR-specific tests
//
// These tests should only describe behavior that is special to FastCKR.
// Generic TreeOracle behavior belongs in test_tree_oracle_common.cpp.
// -----------------------------------------------------------------------------

TEST_CASE("FastCKR samples radius in the interval delta over four to delta over two", "[tree-oracle][fast-ckr][radius]") {
    auto graph = makeSimplePathGraph5();
    FastCKR<FlatHST> oracle(graph);

    const std::vector<int> permutation = identityPermutation(graph.getNumNodes());
    HSTLevel level;

    const double delta = 2.0;

    oracle.computeLevelPartition(graph, level, permutation, delta);

    REQUIRE(level.R >= delta / 4.0);
    REQUIRE(level.R <= delta / 2.0);
}

TEST_CASE("FastCKR radius remains in the valid sampling interval for different deltas", "[tree-oracle][fast-ckr][delta]") {
    auto graph = makeSimplePathGraph5();
    FastCKR<FlatHST> oracle(graph);

    const std::vector<int> permutation = identityPermutation(graph.getNumNodes());

    for (double delta : {0.5, 1.0, 2.0, 5.0}) {
        HSTLevel level;

        oracle.computeLevelPartition(graph, level, permutation, delta);

        REQUIRE(level.R >= delta / 4.0);
        REQUIRE(level.R <= delta / 2.0);
    }
}

TEST_CASE("FastCKR with a full permutation eventually assigns every node", "[tree-oracle][fast-ckr][owners]") {
    auto graph = makeSimplePathGraph5();
    FastCKR<FlatHST> oracle(graph);

    const std::vector<int> permutation = identityPermutation(graph.getNumNodes());
    HSTLevel level;

    oracle.computeLevelPartition(graph, level, permutation, 1.0);

    requireAllOwnersAreValidNodes(level, graph.getNumNodes());
    requireAssignedOwnersAreCenters(level);
}

TEST_CASE("FastCKR may leave nodes unassigned when the permutation does not cover all nodes", "[tree-oracle][fast-ckr][partial-permutation]") {
    auto graph = makeSimplePathGraph5();
    FastCKR<FlatHST> oracle(graph);

    const std::vector<int> permutation = {0};
    HSTLevel level;

    oracle.computeLevelPartition(graph, level, permutation, 1.0);

    REQUIRE(level.centers.size() == 1);
    REQUIRE(level.centers[0] == 0);
    REQUIRE(level.owner.size() == static_cast<size_t>(graph.getNumNodes()));

    requireOwnersAreValidOrUnassigned(level, graph.getNumNodes());
}

TEST_CASE("FastCKR produces valid partitions despite randomized radius sampling", "[tree-oracle][fast-ckr][randomized]") {
    auto graph = makeSimplePathGraph5();
    FastCKR<FlatHST> oracle(graph);

    const std::vector<int> permutation = identityPermutation(graph.getNumNodes());

    for (int trial = 0; trial < 10; ++trial) {
        HSTLevel level;

        oracle.computeLevelPartition(graph, level, permutation, 2.0);

        REQUIRE(!level.centers.empty());
        REQUIRE(level.owner.size() == static_cast<size_t>(graph.getNumNodes()));

        requireCentersAreValidNodes(level, graph.getNumNodes());
        requireCentersAreContainedInPermutation(level, permutation);
        requireAllOwnersAreValidNodes(level, graph.getNumNodes());
        requireAssignedOwnersAreCenters(level);
    }
}

TEST_CASE("FastCKR handles larger grid graphs", "[tree-oracle][fast-ckr][large-graph]") {
    auto graph = makeGridGraph(4);
    FastCKR<FlatHST> oracle(graph);

    const std::vector<int> permutation = identityPermutation(graph.getNumNodes());
    HSTLevel level;

    oracle.computeLevelPartition(graph, level, permutation, 3.0);

    REQUIRE(level.owner.size() == 16);
    REQUIRE(!level.centers.empty());

    requireCentersAreValidNodes(level, graph.getNumNodes());
    requireCentersAreContainedInPermutation(level, permutation);
    requireAllOwnersAreValidNodes(level, graph.getNumNodes());
    requireAssignedOwnersAreCenters(level);
}
