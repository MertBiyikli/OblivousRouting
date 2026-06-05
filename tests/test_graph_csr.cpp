#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <vector>
#include "data_structures/graph/graph_csr.h"

using Catch::Approx;

TEST_CASE("GraphCSR - Basic Construction", "[GraphCSR]") {
    GraphCSR graph(5);
    REQUIRE(graph.getNumNodes() == 5);
    REQUIRE(graph.getNumDirectedEdges() == 0);
}

TEST_CASE("GraphCSR - Edge Addition", "[GraphCSR]") {
    GraphCSR graph(4);
    graph.addEdge(0, 1, 2.0, 1.5);
    graph.addEdge(1, 2, 3.0, 2.5);
    graph.addEdge(2, 3, 4.0, 3.5);

    REQUIRE(graph.getNumUndirectedEdges() == 3);
}

TEST_CASE("GraphCSR - Finalize and Edge Access", "[GraphCSR]") {
    GraphCSR graph(3);
    graph.addEdge(0, 1, 5.0, 1.0);
    graph.addEdge(1, 2, 6.0, 2.0);
    graph.finalize();

    REQUIRE(graph.getEdgeDistance(0, 1) == Approx(1.0));
    REQUIRE(graph.getEdgeCapacity(0, 1) == Approx(5.0));
    REQUIRE(graph.getEdgeDistance(1, 2) == Approx(2.0));
    REQUIRE(graph.getEdgeCapacity(1, 2) == Approx(6.0));
}

TEST_CASE("GraphCSR - Neighbors", "[GraphCSR]") {
    GraphCSR graph(4);
    graph.addEdge(0, 1, 1.0);
    graph.addEdge(0, 2, 1.0);
    graph.addEdge(1, 3, 1.0);
    graph.finalize();

    auto neighbors_0 = graph.neighbors(0);
    REQUIRE(neighbors_0.size() == 2);

    auto neighbors_1 = graph.neighbors(1);
    REQUIRE(neighbors_1.size() == 2);

    auto neighbors_3 = graph.neighbors(3);
    REQUIRE(neighbors_3.size() == 1);
}

TEST_CASE("GraphCSR - Shortest Path Simple", "[GraphCSR]") {
    GraphCSR graph(4);
    graph.addEdge(0, 1, 1.0, 1.0);
    graph.addEdge(1, 2, 1.0, 1.0);
    graph.addEdge(2, 3, 1.0, 1.0);
    graph.finalize();

    auto path = graph.getShortestPath(0, 3);
    REQUIRE(path.size() == 4); // 0 -> 1 -> 2 -> 3
    REQUIRE(path[0] == 0);
    REQUIRE(path[3] == 3);
}

TEST_CASE("GraphCSR - Shortest Path with Multiple Routes", "[GraphCSR]") {
    GraphCSR graph(5);
    // Create diamond graph
    graph.addEdge(0, 1, 1.0, 2.0);
    graph.addEdge(0, 2, 1.0, 2.0);
    graph.addEdge(1, 4, 1.0, 1.0);
    graph.addEdge(2, 4, 1.0, 1.0);
    graph.addEdge(0, 4, 1.0, 5.0); // longer direct route
    graph.finalize();

    auto path = graph.getShortestPath(0, 4);
    REQUIRE(!path.empty());
    REQUIRE(path.front() == 0);
    REQUIRE(path.back() == 4);
}

TEST_CASE("GraphCSR - Edge Distance Update", "[GraphCSR]") {
    GraphCSR graph(3);
    graph.addEdge(0, 1, 5.0, 1.0);
    graph.finalize();

    REQUIRE(graph.getEdgeDistance(0, 1) == Approx(1.0));
    graph.updateEdgeDistance(0, 1, 3.0);
    REQUIRE(graph.getEdgeDistance(0, 1) == Approx(3.0));
}

TEST_CASE("GraphCSR - Bidirectional Dijkstra", "[GraphCSR]") {
    GraphCSR graph(6);
    graph.addEdge(0, 1, 1.0, 1.0);
    graph.addEdge(1, 2, 1.0, 1.0);
    graph.addEdge(2, 3, 1.0, 1.0);
    graph.addEdge(3, 4, 1.0, 1.0);
    graph.addEdge(4, 5, 1.0, 1.0);
    graph.finalize();

    auto path = graph.getShortestPathBidirectionalSearch(0, 5);
    REQUIRE(!path.empty());
    REQUIRE(path.front() == 0);
    REQUIRE(path.back() == 5);
    REQUIRE(path.size() == 6); // Should be shortest path
}

TEST_CASE("GraphCSR - Custom Distance Vector", "[GraphCSR]") {
    GraphCSR graph(4);
    graph.addEdge(0, 1, 1.0, 1.0);
    graph.addEdge(1, 2, 1.0, 2.0);
    graph.addEdge(0, 2, 1.0, 5.0);
    graph.finalize();

    // Create custom distance vector that overrides edge distances
    std::vector<double> custom_dist(graph.getNumDirectedEdges(), 1.0);
    auto path = graph.getShortestPath(0, 2, custom_dist);
    REQUIRE(!path.empty());
}


TEST_CASE("GraphCSR - Edge Endpoints", "[GraphCSR]") {
    GraphCSR graph(3);
    graph.addEdge(0, 1, 1.0);
    graph.addEdge(1, 2, 1.0);
    graph.finalize();

    auto endpoints = graph.getEdgeEndpoints(0);
    REQUIRE(endpoints.first == 0);
    REQUIRE(endpoints.second == 1);
}

TEST_CASE("GraphCSR - Large Graph Performance", "[GraphCSR]") {
    GraphCSR graph(100);
    // Create a grid-like structure
    for (int i = 0; i < 100; ++i) {
        if (i + 1 < 100) {
            graph.addEdge(i, i + 1, 1.0, 1.0);
        }
        if (i + 10 < 100) {
            graph.addEdge(i, i + 10, 1.0, 1.0);
        }
    }
    graph.finalize();

    REQUIRE(graph.getNumNodes() == 100);
    auto path = graph.getShortestPath(0, 99);
    REQUIRE(!path.empty());
}

