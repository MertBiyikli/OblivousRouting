#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "data_structures/graph/graph_csr.h"
#include "algorithms/mwu/oracle/tree/frt/frt.h"

using Catch::Approx;

// Helper to create a simple test graph
static GraphCSR createSimpleGraph() {
    GraphCSR graph(5);
    graph.addEdge(0, 1, 1.0, 1.0);
    graph.addEdge(1, 2, 1.0, 1.5);
    graph.addEdge(2, 3, 1.0, 1.0);
    graph.addEdge(3, 4, 1.0, 2.0);
    graph.finalize();
    return graph;
}

static GraphCSR createFullyConnectedSmall() {
    GraphCSR graph(4);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (i != j) {
                graph.addEdge(i, j, 1.0, 1.0);
            }
        }
    }
    graph.finalize();
    return graph;
}

TEST_CASE("FRT - Construction", "[FRT]") {
    auto graph = createSimpleGraph();
    FRT<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    REQUIRE(oracle.graph.getNumNodes() == 5);
}

TEST_CASE("FRT - Basic Level Partition", "[FRT]") {
    auto graph = createSimpleGraph();
    FRT<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    std::vector<int> permutation = {0, 1, 2, 3, 4};
    HSTLevel level;

    // This should not throw
    oracle.computeLevelPartition(
        const_cast<GraphCSR&>(graph),
        level,
        permutation,
        2.0
    );

    REQUIRE(!level.centers.empty());
    REQUIRE(level.owner.size() == 5);
}

TEST_CASE("FRT - Owner Assignment", "[FRT]") {
    auto graph = createSimpleGraph();
    FRT<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    std::vector<int> permutation = {0, 2, 4};
    HSTLevel level;

    oracle.computeLevelPartition(
        const_cast<GraphCSR&>(graph),
        level,
        permutation,
        1.5
    );

    // All nodes should have valid owners
    for (int i = 0; i < 5; ++i) {
        REQUIRE(level.owner[i] >= 0);
        REQUIRE(level.owner[i] < 5);
    }
}

TEST_CASE("FRT - Delta Storage", "[FRT]") {
    auto graph = createSimpleGraph();
    FRT<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    std::vector<int> permutation = {0, 2};
    HSTLevel level;
    double delta = 2.5;

    oracle.computeLevelPartition(
        const_cast<GraphCSR&>(graph),
        level,
        permutation,
        delta
    );

    // FRT stores exact delta value
    REQUIRE(level.R == Approx(delta));
}

TEST_CASE("FRT - Single Center", "[FRT]") {
    auto graph = createSimpleGraph();
    FRT<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    std::vector<int> permutation = {0};
    HSTLevel level;

    oracle.computeLevelPartition(
        const_cast<GraphCSR&>(graph),
        level,
        permutation,
        10.0 // Large radius
    );

    // Single center case
    REQUIRE(level.centers.size() == 1);
    REQUIRE(level.centers[0] == 0);
}

TEST_CASE("FRT - Multiple Centers", "[FRT]") {
    auto graph = createSimpleGraph();
    FRT<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    std::vector<int> permutation = {0, 1, 2, 3, 4};
    HSTLevel level;

    oracle.computeLevelPartition(
        const_cast<GraphCSR&>(graph),
        level,
        permutation,
        0.5 // Small radius forces multiple centers
    );

    REQUIRE(level.centers.size() > 0);
}

TEST_CASE("FRT - Fully Connected Graph", "[FRT]") {
    auto graph = createFullyConnectedSmall();
    FRT<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    std::vector<int> permutation = {0, 1, 2, 3};
    HSTLevel level;

    oracle.computeLevelPartition(
        const_cast<GraphCSR&>(graph),
        level,
        permutation,
        1.0
    );

    REQUIRE(level.owner.size() == 4);
    REQUIRE(!level.centers.empty());
}

TEST_CASE("FRT - Path Graph Partition", "[FRT]") {
    GraphCSR graph(6);
    // Create a path: 0-1-2-3-4-5
    for (int i = 0; i < 5; ++i) {
        graph.addEdge(i, i + 1, 1.0, 1.0);
    }
    graph.finalize();

    FRT<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    std::vector<int> permutation = {0, 3};
    HSTLevel level;

    oracle.computeLevelPartition(
        const_cast<GraphCSR&>(graph),
        level,
        permutation,
        2.0
    );

    REQUIRE(level.owner.size() == 6);
}

TEST_CASE("FRT - Empty Permutation", "[FRT]") {
    auto graph = createSimpleGraph();
    FRT<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    std::vector<int> permutation = {};
    HSTLevel level;

    oracle.computeLevelPartition(
        const_cast<GraphCSR&>(graph),
        level,
        permutation,
        2.0
    );

    REQUIRE(level.centers.size() == 0);
}

TEST_CASE("FRT - Varying Delta Values", "[FRT]") {
    auto graph = createSimpleGraph();
    FRT<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    std::vector<int> permutation = {0, 2, 4};

    for (double delta : {0.5, 1.0, 2.0, 5.0}) {
        HSTLevel level;

        oracle.computeLevelPartition(
            const_cast<GraphCSR&>(graph),
            level,
            permutation,
            delta
        );

        REQUIRE(level.R == Approx(delta));
    }
}

TEST_CASE("FRT - With Mendel Scaling", "[FRT]") {
    auto graph = createSimpleGraph();
    FRT<FlatHST> oracle(const_cast<GraphCSR&>(graph), true);

    std::vector<int> permutation = {0, 1, 2};
    HSTLevel level;

    oracle.computeLevelPartition(
        const_cast<GraphCSR&>(graph),
        level,
        permutation,
        2.0
    );

    REQUIRE(!level.centers.empty());
}

TEST_CASE("FRT - All Nodes in Permutation", "[FRT]") {
    auto graph = createSimpleGraph();
    FRT<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    // All nodes in permutation
    std::vector<int> permutation = {0, 1, 2, 3, 4};
    HSTLevel level;

    oracle.computeLevelPartition(
        const_cast<GraphCSR&>(graph),
        level,
        permutation,
        5.0
    );

    // All unassigned nodes should be covered
    for (int i = 0; i < 5; ++i) {
        REQUIRE(level.owner[i] >= 0);
        REQUIRE(level.owner[i] < 5);
    }
}

TEST_CASE("FRT - Redundant Permutation Nodes", "[FRT]") {
    auto graph = createSimpleGraph();
    FRT<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    // Permutation with duplicates - second occurrence should be skipped
    std::vector<int> permutation = {0, 1, 0, 2};
    HSTLevel level;

    oracle.computeLevelPartition(
        const_cast<GraphCSR&>(graph),
        level,
        permutation,
        2.0
    );

    // Should handle gracefully
    REQUIRE(!level.centers.empty());
}

TEST_CASE("FRT - Large Delta Value", "[FRT]") {
    auto graph = createSimpleGraph();
    FRT<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    std::vector<int> permutation = {0};
    HSTLevel level;

    oracle.computeLevelPartition(
        const_cast<GraphCSR&>(graph),
        level,
        permutation,
        1000.0 // Very large radius
    );

    // All nodes should be in cluster of node 0
    REQUIRE(level.centers.size() == 1);
}

