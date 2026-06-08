#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <vector>

#include "tree_oracle_test_helpers.h"
#include "catch2/catch_template_test_macros.hpp"

using Catch::Approx;

// -----------------------------------------------------------------------------
// Common tests for all tree oracles.
//
// This file now includes:
//   - FRT
//   - FastCKR
//   - TreeMST
//
// Important design distinction:
//   - FRT and FastCKR are partition-based.
//   - TreeMST is not partition-based; it overrides getTree(...) directly.
// -----------------------------------------------------------------------------

TEMPLATE_TEST_CASE(
    "TreeOracle implementations can be constructed from a graph",
    "[tree-oracle][common][construction]",
    FRTOracleCase,
    FastCKROracleCase,
    MSTOracleCase
) {
    auto graph = makeSimplePathGraph5();

    typename TestType::Oracle oracle(graph);

    REQUIRE(oracle.graph.getNumNodes() == 5);
}

TEMPLATE_TEST_CASE(
    "TreeOracle implementations keep a reference to the original graph",
    "[tree-oracle][common][graph-reference]",
    FRTOracleCase,
    FastCKROracleCase,
    MSTOracleCase
) {
    auto graph = makeGridGraph(3);

    typename TestType::Oracle oracle(graph);

    REQUIRE(oracle.graph.getNumNodes() == 9);
    REQUIRE(oracle.graph.getNumUndirectedEdges() == graph.getNumUndirectedEdges());
    REQUIRE(oracle.graph.getNumDirectedEdges() == graph.getNumDirectedEdges());
}

TEMPLATE_TEST_CASE(
    "TreeOracle implementations can build a flat tree using current graph distances",
    "[tree-oracle][common][get-tree]",
    FRTOracleCase,
    FastCKROracleCase,
    MSTOracleCase
) {
    auto graph = makeSimplePathGraph5();

    typename TestType::Oracle oracle(graph);

    auto distances = currentGraphDistances(graph);

    REQUIRE_NOTHROW(oracle.getTree(distances));
}

TEMPLATE_TEST_CASE(
    "TreeOracle implementations can build a flat tree using unit distances",
    "[tree-oracle][common][get-tree]",
    FRTOracleCase,
    FastCKROracleCase,
    MSTOracleCase
) {
    auto graph = makePathGraph(6);

    typename TestType::Oracle oracle(graph);

    auto distances = unitDistances(graph);

    REQUIRE_NOTHROW(oracle.getTree(distances));
}

TEMPLATE_TEST_CASE(
    "TreeOracle implementations can build a flat tree on a grid graph",
    "[tree-oracle][common][get-tree][grid]",
    FRTOracleCase,
    FastCKROracleCase,
    MSTOracleCase
) {
    auto graph = makeGridGraph(3);

    typename TestType::Oracle oracle(graph);

    auto distances = currentGraphDistances(graph);

    REQUIRE_NOTHROW(oracle.getTree(distances));
}

TEMPLATE_TEST_CASE(
    "TreeOracle implementations update graph distances through getTree",
    "[tree-oracle][common][distance-update]",
    FRTOracleCase,
    FastCKROracleCase,
    MSTOracleCase
) {
    auto graph = makeSimplePathGraph5();

    typename TestType::Oracle oracle(graph);

    auto distances = increasingDistances(graph);

    REQUIRE_NOTHROW(oracle.getTree(distances));

    for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
        REQUIRE(graph.getEdgeDistance(e) == Approx(distances[e]));
    }
}

TEMPLATE_TEST_CASE(
    "TreeOracle implementations populate scales after getTree",
    "[tree-oracle][common][scales]",
    FRTOracleCase,
    FastCKROracleCase,
    MSTOracleCase
) {
    auto graph = makeSimplePathGraph5();

    typename TestType::Oracle oracle(graph);

    auto distances = currentGraphDistances(graph);

    REQUIRE(oracle.scales.empty());

    REQUIRE_NOTHROW(oracle.getTree(distances));

    REQUIRE(!oracle.scales.empty());
}

TEMPLATE_TEST_CASE(
    "Partition-based TreeOracle implementations compute a basic level partition",
    "[tree-oracle][common][level-partition]",
    FRTOracleCase,
    FastCKROracleCase
) {
    auto graph = makeSimplePathGraph5();

    typename TestType::Oracle oracle(graph);

    const std::vector<int> permutation = identityPermutation(graph.getNumNodes());
    HSTLevel level;

    oracle.computeLevelPartition(graph, level, permutation, 2.0);

    REQUIRE(!level.centers.empty());
    requirePartitionShapeIsValid(level, graph.getNumNodes(), permutation);
}

TEMPLATE_TEST_CASE(
    "Partition-based TreeOracle implementations produce unique centers",
    "[tree-oracle][common][centers]",
    FRTOracleCase,
    FastCKROracleCase
) {
    auto graph = makeSimplePathGraph5();

    typename TestType::Oracle oracle(graph);

    const std::vector<int> permutation = identityPermutation(graph.getNumNodes());
    HSTLevel level;

    oracle.computeLevelPartition(graph, level, permutation, 2.0);

    REQUIRE(hasUniqueValues(level.centers));
}

TEMPLATE_TEST_CASE(
    "Partition-based TreeOracle implementations choose centers only from the permutation",
    "[tree-oracle][common][centers]",
    FRTOracleCase,
    FastCKROracleCase
) {
    auto graph = makeSimplePathGraph5();

    typename TestType::Oracle oracle(graph);

    const std::vector<int> permutation = {0, 2, 4};
    HSTLevel level;

    oracle.computeLevelPartition(graph, level, permutation, 2.0);

    requireCentersAreContainedInPermutation(level, permutation);
}

TEMPLATE_TEST_CASE(
    "Partition-based TreeOracle implementations assign every node when the full permutation is used",
    "[tree-oracle][common][owners]",
    FRTOracleCase,
    FastCKROracleCase
) {
    auto graph = makeSimplePathGraph5();

    typename TestType::Oracle oracle(graph);

    const std::vector<int> permutation = identityPermutation(graph.getNumNodes());
    HSTLevel level;

    oracle.computeLevelPartition(graph, level, permutation, 2.0);

    requireAllOwnersAreValidNodes(level, graph.getNumNodes());
}

TEMPLATE_TEST_CASE(
    "Partition-based TreeOracle implementations assign owners that correspond to chosen centers",
    "[tree-oracle][common][owners]",
    FRTOracleCase,
    FastCKROracleCase
) {
    auto graph = makeSimplePathGraph5();

    typename TestType::Oracle oracle(graph);

    const std::vector<int> permutation = identityPermutation(graph.getNumNodes());
    HSTLevel level;

    oracle.computeLevelPartition(graph, level, permutation, 2.0);

    requireAssignedOwnersAreCenters(level);
}

TEMPLATE_TEST_CASE(
    "Partition-based TreeOracle implementations handle an empty permutation without creating centers",
    "[tree-oracle][common][empty-permutation]",
    FRTOracleCase,
    FastCKROracleCase
) {
    auto graph = makeSimplePathGraph5();

    typename TestType::Oracle oracle(graph);

    const std::vector<int> permutation = {};
    HSTLevel level;

    oracle.computeLevelPartition(graph, level, permutation, 2.0);

    REQUIRE(level.centers.empty());
    REQUIRE(level.owner.size() == static_cast<size_t>(graph.getNumNodes()));
}

TEMPLATE_TEST_CASE(
    "Partition-based TreeOracle implementations handle duplicate permutation entries gracefully",
    "[tree-oracle][common][permutation]",
    FRTOracleCase,
    FastCKROracleCase
) {
    auto graph = makeSimplePathGraph5();

    typename TestType::Oracle oracle(graph);

    const std::vector<int> permutation = {0, 1, 0, 2, 2, 4};
    HSTLevel level;

    oracle.computeLevelPartition(graph, level, permutation, 2.0);

    REQUIRE(!level.centers.empty());
    requirePartitionShapeIsValid(level, graph.getNumNodes(), permutation);
}

TEMPLATE_TEST_CASE(
    "Partition-based TreeOracle implementations work for different delta values",
    "[tree-oracle][common][delta]",
    FRTOracleCase,
    FastCKROracleCase
) {
    auto graph = makeSimplePathGraph5();

    typename TestType::Oracle oracle(graph);

    const std::vector<int> permutation = identityPermutation(graph.getNumNodes());

    for (double delta : {0.5, 1.0, 2.0, 5.0}) {
        HSTLevel level;

        oracle.computeLevelPartition(graph, level, permutation, delta);

        REQUIRE(!level.centers.empty());
        requirePartitionShapeIsValid(level, graph.getNumNodes(), permutation);
    }
}

TEMPLATE_TEST_CASE(
    "Partition-based TreeOracle implementations work on a grid graph",
    "[tree-oracle][common][grid]",
    FRTOracleCase,
    FastCKROracleCase
) {
    auto graph = makeGridGraph(4);

    typename TestType::Oracle oracle(graph);

    const std::vector<int> permutation = identityPermutation(graph.getNumNodes());
    HSTLevel level;

    oracle.computeLevelPartition(graph, level, permutation, 3.0);

    REQUIRE(!level.centers.empty());
    requirePartitionShapeIsValid(level, graph.getNumNodes(), permutation);
}

TEMPLATE_TEST_CASE(
    "Partition-based TreeOracle implementations work on a fully connected graph",
    "[tree-oracle][common][complete-graph]",
    FRTOracleCase,
    FastCKROracleCase
) {
    auto graph = makeFullyConnectedSmallGraph();

    typename TestType::Oracle oracle(graph);

    const std::vector<int> permutation = identityPermutation(graph.getNumNodes());
    HSTLevel level;

    oracle.computeLevelPartition(graph, level, permutation, 1.0);

    REQUIRE(!level.centers.empty());
    requirePartitionShapeIsValid(level, graph.getNumNodes(), permutation);
}

TEMPLATE_TEST_CASE(
    "Partition-based TreeOracle implementations support Mendel-scaling construction flag",
    "[tree-oracle][common][mendel-scaling]",
    FRTOracleCase,
    FastCKROracleCase
) {
    auto graph = makeSimplePathGraph5();

    typename TestType::Oracle oracle(graph, true);

    REQUIRE(oracle.applyMendelScaling);
    REQUIRE(oracle.graph.getNumNodes() == 5);
}

TEMPLATE_TEST_CASE(
    "Partition-based TreeOracle implementations can compute a level partition when Mendel scaling is enabled",
    "[tree-oracle][common][mendel-scaling]",
    FRTOracleCase,
    FastCKROracleCase
) {
    auto graph = makeSimplePathGraph5();

    typename TestType::Oracle oracle(graph, true);

    const std::vector<int> permutation = identityPermutation(graph.getNumNodes());
    HSTLevel level;

    oracle.computeLevelPartition(graph, level, permutation, 2.0);

    REQUIRE(!level.centers.empty());
    requirePartitionShapeIsValid(level, graph.getNumNodes(), permutation);
}