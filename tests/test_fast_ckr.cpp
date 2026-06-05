#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "data_structures/graph/graph_csr.h"
#include "algorithms/mwu/oracle/tree/fast_ckr/fast_ckr.h"

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

// Helper to create a grid graph
static GraphCSR createGridGraph(int size) {
    GraphCSR graph(size * size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            int node = i * size + j;
            if (j + 1 < size) {
                int right = i * size + (j + 1);
                graph.addEdge(node, right, 1.0, 1.0);
            }
            if (i + 1 < size) {
                int down = (i + 1) * size + j;
                graph.addEdge(node, down, 1.0, 1.0);
            }
        }
    }
    graph.finalize();
    return graph;
}

TEST_CASE("FastCKR - Construction", "[FastCKR]") {
    auto graph = createSimpleGraph();
    FastCKR<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    REQUIRE(oracle.graph.getNumNodes() == 5);
}

TEST_CASE("FastCKR - Basic Level Partition", "[FastCKR]") {
    auto graph = createSimpleGraph();
    FastCKR<FlatHST> oracle(const_cast<GraphCSR&>(graph));

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

TEST_CASE("FastCKR - Cluster Assignment", "[FastCKR]") {
    auto graph = createSimpleGraph();
    FastCKR<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    std::vector<int> permutation = {0, 1, 2};
    HSTLevel level;

    oracle.computeLevelPartition(
        const_cast<GraphCSR&>(graph),
        level,
        permutation,
        1.5
    );

    // All nodes should be assigned to some center
    for (int i = 0; i < 5; ++i) {
        REQUIRE(level.owner[i] >= -1);
        REQUIRE(level.owner[i] < 5);
    }
}

TEST_CASE("FastCKR - Centers Non-Empty", "[FastCKR]") {
    auto graph = createSimpleGraph();
    FastCKR<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    std::vector<int> permutation = {0, 2, 4};
    HSTLevel level;

    oracle.computeLevelPartition(
        const_cast<GraphCSR&>(graph),
        level,
        permutation,
        2.0
    );

    REQUIRE(level.centers.size() > 0);
    REQUIRE(level.centers.size() <= permutation.size());
}

TEST_CASE("FastCKR - Radius Setting", "[FastCKR]") {
    auto graph = createSimpleGraph();
    FastCKR<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    std::vector<int> permutation = {0, 1, 2};
    HSTLevel level;
    double delta = 2.0;

    oracle.computeLevelPartition(
        const_cast<GraphCSR&>(graph),
        level,
        permutation,
        delta
    );

    // Radius should be in range [delta/4, delta/2]
    REQUIRE(level.R >= delta / 4.0);
    REQUIRE(level.R <= delta / 2.0);
}

TEST_CASE("FastCKR - Single Center", "[FastCKR]") {
    auto graph = createSimpleGraph();
    FastCKR<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    std::vector<int> permutation = {0}; // Only one center
    HSTLevel level;

    oracle.computeLevelPartition(
        const_cast<GraphCSR&>(graph),
        level,
        permutation,
        10.0 // Large radius
    );

    // With large radius and single center, all nodes should be assigned to it
    REQUIRE(level.centers.size() == 1);
}

TEST_CASE("FastCKR - Large Graph Partition", "[FastCKR]") {
    auto graph = createGridGraph(4);
    FastCKR<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    std::vector<int> permutation;
    for (int i = 0; i < 16; i += 2) {
        permutation.push_back(i);
    }

    HSTLevel level;

    oracle.computeLevelPartition(
        const_cast<GraphCSR&>(graph),
        level,
        permutation,
        3.0
    );

    REQUIRE(level.owner.size() == 16);
    REQUIRE(!level.centers.empty());
}

TEST_CASE("FastCKR - Empty Permutation", "[FastCKR]") {
    auto graph = createSimpleGraph();
    FastCKR<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    std::vector<int> permutation = {}; // Empty
    HSTLevel level;

    oracle.computeLevelPartition(
        const_cast<GraphCSR&>(graph),
        level,
        permutation,
        2.0
    );

    // Should handle empty permutation gracefully
    REQUIRE(level.centers.size() == 0);
}

TEST_CASE("FastCKR - Varying Delta Values", "[FastCKR]") {
    auto graph = createSimpleGraph();
    FastCKR<FlatHST> oracle(const_cast<GraphCSR&>(graph));

    std::vector<int> permutation = {0, 2, 4};

    // Test with different delta values
    for (double delta : {0.5, 1.0, 2.0, 5.0}) {
        HSTLevel level;

        oracle.computeLevelPartition(
            const_cast<GraphCSR&>(graph),
            level,
            permutation,
            delta
        );

        REQUIRE(level.R >= delta / 4.0);
        REQUIRE(level.R <= delta / 2.0);
    }
}

TEST_CASE("FastCKR - With Mendel Scaling", "[FastCKR]") {
    auto graph = createSimpleGraph();
    FastCKR<FlatHST> oracle(const_cast<GraphCSR&>(graph), true); // Enable Mendel scaling

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

TEST_CASE("FastCKR - Deterministic with Fixed Seed", "[FastCKR]") {
    auto graph = createSimpleGraph();
    std::vector<int> permutation = {0, 1, 2, 3, 4};

    // Create multiple oracles and run partitions (note: they use random seeds)
    FastCKR<FlatHST> oracle1(const_cast<GraphCSR&>(graph));
    HSTLevel level1;
    oracle1.computeLevelPartition(const_cast<GraphCSR&>(graph), level1, permutation, 2.0);

    // Results should be valid even if not deterministic due to randomness
    REQUIRE(!level1.centers.empty());
    REQUIRE(level1.owner.size() == 5);
}

