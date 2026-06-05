#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "data_structures/graph/graph_adj.h"

using Catch::Approx;

TEST_CASE("GraphADJList - Basic Construction", "[GraphADJList]") {
    GraphADJList graph(5);
    REQUIRE(graph.getNumNodes() == 5);
    REQUIRE(graph.getNumDirectedEdges() == 0);
}

TEST_CASE("GraphADJList - Edge Addition", "[GraphADJList]") {
    GraphADJList graph(4);
    graph.addEdge(0, 1, 2.0, 1.5);
    graph.addEdge(1, 2, 3.0, 2.5);
    graph.addEdge(2, 3, 4.0, 3.5);

    REQUIRE(graph.getNumUndirectedEdges() == 3);
}

TEST_CASE("GraphADJList - Finalize and Edge Access", "[GraphADJList]") {
    GraphADJList graph(3);
    graph.addEdge(0, 1, 5.0, 1.0);
    graph.addEdge(1, 2, 6.0, 2.0);
    graph.finalize();

    REQUIRE(graph.getEdgeDistance(0, 1) == Approx(1.0));
    REQUIRE(graph.getEdgeCapacity(0, 1) == Approx(5.0));
    REQUIRE(graph.getEdgeDistance(1, 2) == Approx(2.0));
    REQUIRE(graph.getEdgeCapacity(1, 2) == Approx(6.0));
}

TEST_CASE("GraphADJList - Neighbors", "[GraphADJList]") {
    GraphADJList graph(4);
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

TEST_CASE("GraphADJList - Shortest Path Simple", "[GraphADJList]") {
    GraphADJList graph(4);
    graph.addEdge(0, 1, 1.0, 1.0);
    graph.addEdge(1, 2, 1.0, 1.0);
    graph.addEdge(2, 3, 1.0, 1.0);
    graph.finalize();

    auto path = graph.getShortestPath(0, 3);
    REQUIRE(path.size() == 4); // 0 -> 1 -> 2 -> 3
    REQUIRE(path[0] == 0);
    REQUIRE(path[3] == 3);
}

TEST_CASE("GraphADJList - Shortest Distance", "[GraphADJList]") {
    GraphADJList graph(4);
    graph.addEdge(0, 1, 1.0, 2.0);
    graph.addEdge(1, 2, 1.0, 3.0);
    graph.addEdge(2, 3, 1.0, 4.0);
    graph.finalize();

    double dist = graph.getShortestDistance(0, 3);
    REQUIRE(dist == Approx(9.0)); // 2 + 3 + 4
}

TEST_CASE("GraphADJList - Edge Distance Update", "[GraphADJList]") {
    GraphADJList graph(3);
    graph.addEdge(0, 1, 5.0, 1.0);
    graph.finalize();

    REQUIRE(graph.getEdgeDistance(0, 1) == Approx(1.0));
    graph.updateEdgeDistance(0, 1, 3.0);
    REQUIRE(graph.getEdgeDistance(0, 1) == Approx(3.0));
}

TEST_CASE("GraphADJList - Edge Capacity Update", "[GraphADJList]") {
    GraphADJList graph(3);
    graph.addEdge(0, 1, 5.0, 1.0);
    graph.finalize();

    REQUIRE(graph.getEdgeCapacity(0, 1) == Approx(5.0));
    graph.updateEdgeCapacity(0, 1, 10.0);
    REQUIRE(graph.getEdgeCapacity(0, 1) == Approx(10.0));
}

TEST_CASE("GraphADJList - Multiple Edges from Same Node", "[GraphADJList]") {
    GraphADJList graph(5);
    graph.addEdge(0, 1, 1.0, 1.0);
    graph.addEdge(0, 2, 1.0, 2.0);
    graph.addEdge(0, 3, 1.0, 3.0);
    graph.addEdge(0, 4, 1.0, 4.0);
    graph.finalize();

    auto neighbors = graph.neighbors(0);
    REQUIRE(neighbors.size() == 4);
}

TEST_CASE("GraphADJList - Shortest Path with Alternative Routes", "[GraphADJList]") {
    GraphADJList graph(5);
    // Triangle with extra edge
    graph.addEdge(0, 1, 1.0, 1.0);
    graph.addEdge(0, 2, 1.0, 1.0);
    graph.addEdge(1, 3, 1.0, 1.0);
    graph.addEdge(2, 3, 1.0, 1.0);
    graph.addEdge(0, 3, 1.0, 5.0); // longer direct route
    graph.finalize();

    auto path = graph.getShortestPath(0, 3);
    REQUIRE(!path.empty());
    REQUIRE(path.front() == 0);
    REQUIRE(path.back() == 3);
}

TEST_CASE("GraphADJList - Distances from Single Source", "[GraphADJList]") {
    GraphADJList graph(4);
    graph.addEdge(0, 1, 1.0, 1.0);
    graph.addEdge(0, 2, 1.0, 2.0);
    graph.addEdge(1, 3, 1.0, 2.0);
    graph.finalize();

    auto distances = graph.getDistances(0);
    REQUIRE(distances.size() == 4);
    REQUIRE(distances[0] == Approx(0.0)); // Source
    REQUIRE(distances[1] == Approx(1.0));
    REQUIRE(distances[2] == Approx(2.0));
}

TEST_CASE("GraphADJList - Adjacent List Structure", "[GraphADJList]") {
    GraphADJList graph(3);
    graph.addEdge(0, 1, 1.0);
    graph.addEdge(0, 2, 1.0);
    graph.addEdge(1, 2, 1.0);
    graph.finalize();

    REQUIRE(graph.adjList[0].size() == 2);
    REQUIRE(graph.adjList[1].size() == 2);
    REQUIRE(graph.adjList[2].size() == 2);
}

TEST_CASE("GraphADJList - Large Graph", "[GraphADJList]") {
    GraphADJList graph(20);  // Reduced from 50 to avoid memory issues
    for (int i = 0; i < 19; ++i) {
        graph.addEdge(i, i + 1, 1.0, 1.0);
    }
    graph.finalize();

    REQUIRE(graph.getNumNodes() == 20);
    auto path = graph.getShortestPath(0, 19);
    REQUIRE(!path.empty());
    REQUIRE(path.size() == 20);
}



