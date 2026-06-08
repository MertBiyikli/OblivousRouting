#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "data_structures/graph/graph_adj.h"

#include "graph_test_helpers.h"
using Catch::Approx;

TEST_CASE("GraphADJList stores the expected adjacency-list structure", "[graph][adj][structure]") {
    GraphADJList graph(3);
    graph.addEdge(0, 1, 1.0);
    graph.addEdge(0, 2, 1.0);
    graph.addEdge(1, 2, 1.0);
    graph.finalize();

    REQUIRE(graph.adjList.size() == 3);
    REQUIRE(graph.adjList[0].size() == 2);
    REQUIRE(graph.adjList[1].size() == 2);
    REQUIRE(graph.adjList[2].size() == 2);
}

TEST_CASE("GraphADJList updates edge capacity by node pair", "[graph][adj][capacity-update]") {
    GraphADJList graph(3);
    graph.addEdge(0, 1, 5.0, 1.0);
    graph.finalize();

    REQUIRE(graph.getEdgeCapacity(0, 1) == Approx(5.0));
    REQUIRE(graph.updateEdgeCapacity(0, 1, 10.0));
    REQUIRE(graph.getEdgeCapacity(0, 1) == Approx(10.0));
}

TEST_CASE("GraphADJList computes all distances from a single source", "[graph][adj][single-source-distances]") {
    GraphADJList graph(4);
    graph.addEdge(0, 1, 1.0, 1.0);
    graph.addEdge(0, 2, 1.0, 2.0);
    graph.addEdge(1, 3, 1.0, 2.0);
    graph.finalize();

    auto distances = graph.getDistances(0);

    REQUIRE(distances.size() == 4);
    REQUIRE(distances[0] == Approx(0.0));
    REQUIRE(distances[1] == Approx(1.0));
    REQUIRE(distances[2] == Approx(2.0));
    REQUIRE(distances[3] == Approx(3.0));
}

