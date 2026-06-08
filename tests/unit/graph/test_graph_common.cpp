//
// Created by Mert Biyikli on 08.06.26.
//
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <vector>
#include <cmath>

#include "graph_test_helpers.h"
#include "catch2/catch_template_test_macros.hpp"
#include "data_structures/graph/graph_adj.h"
#include "data_structures/graph/graph_csr.h"
#include "data_structures/graph/Igraph.h"

using Catch::Approx;

// These tests define the common behavioral contract of every IGraph implementation.
// If a new graph implementation is added later, add it to this list and it must pass
// the same interface-level tests.

TEMPLATE_TEST_CASE(
    "IGraph implementations start with the requested number of nodes and no edges",
    "[graph][common][construction]",
    GraphADJList,
    GraphCSR
) {
    TestType graph(5);

    REQUIRE(graph.getNumNodes() == 5);
    REQUIRE(graph.getVertices().size() == 5);
    REQUIRE(graph.getVertices()[0] == 0);
    REQUIRE(graph.getVertices()[4] == 4);
    REQUIRE(graph.getNumDirectedEdges() == 0);
    REQUIRE(graph.getNumUndirectedEdges() == 0);
}

TEMPLATE_TEST_CASE(
    "IGraph implementations count undirected and directed edges consistently",
    "[graph][common][edge-count]",
    GraphADJList,
    GraphCSR
) {
    TestType graph(4);

    graph.addEdge(0, 1, 2.0, 1.5);
    graph.addEdge(1, 2, 3.0, 2.5);
    graph.addEdge(2, 3, 4.0, 3.5);
    graph.finalize();

    REQUIRE(graph.getNumUndirectedEdges() == 3);
    REQUIRE(graph.getNumDirectedEdges() == 6);
}

TEMPLATE_TEST_CASE(
    "IGraph implementations expose capacities and distances by node pair",
    "[graph][common][node-edge-access]",
    GraphADJList,
    GraphCSR
) {
    TestType graph(3);

    graph.addEdge(0, 1, 5.0, 1.0);
    graph.addEdge(1, 2, 6.0, 2.0);
    graph.finalize();

    REQUIRE(graph.getEdgeCapacity(0, 1) == Approx(5.0));
    REQUIRE(graph.getEdgeDistance(0, 1) == Approx(1.0));
    REQUIRE(graph.getEdgeCapacity(1, 2) == Approx(6.0));
    REQUIRE(graph.getEdgeDistance(1, 2) == Approx(2.0));

    REQUIRE(graph.getEdgeCapacity(1, 0) == Approx(5.0));
    REQUIRE(graph.getEdgeDistance(1, 0) == Approx(1.0));
}

TEMPLATE_TEST_CASE(
    "IGraph implementations expose valid neighbor ranges",
    "[graph][common][neighbors]",
    GraphADJList,
    GraphCSR
) {
    TestType graph(4);

    graph.addEdge(0, 1, 1.0);
    graph.addEdge(0, 2, 1.0);
    graph.addEdge(1, 3, 1.0);
    graph.finalize();

    auto neighbors_0 = graph.neighbors(0);
    auto neighbors_1 = graph.neighbors(1);
    auto neighbors_2 = graph.neighbors(2);
    auto neighbors_3 = graph.neighbors(3);

    REQUIRE(neighbors_0.size() == 2);
    REQUIRE(containsNeighbor(neighbors_0, 1));
    REQUIRE(containsNeighbor(neighbors_0, 2));

    REQUIRE(neighbors_1.size() == 2);
    REQUIRE(containsNeighbor(neighbors_1, 0));
    REQUIRE(containsNeighbor(neighbors_1, 3));

    REQUIRE(neighbors_2.size() == 1);
    REQUIRE(containsNeighbor(neighbors_2, 0));

    REQUIRE(neighbors_3.size() == 1);
    REQUIRE(containsNeighbor(neighbors_3, 1));
}

TEMPLATE_TEST_CASE(
    "IGraph implementations compute shortest paths on a simple path graph",
    "[graph][common][shortest-path]",
    GraphADJList,
    GraphCSR
) {
    auto graph = makePathGraph4<TestType>();

    auto path = graph.getShortestPath(0, 3);

    REQUIRE(path == std::vector<int>{0, 1, 2, 3});
}

TEMPLATE_TEST_CASE(
    "IGraph implementations choose a cheaper indirect path over an expensive direct edge",
    "[graph][common][shortest-path]",
    GraphADJList,
    GraphCSR
) {
    auto graph = makeDiamondGraph<TestType>();

    auto path = graph.getShortestPath(0, 4);

    REQUIRE(isValidPathFromTo(path, 0, 4));
    REQUIRE(path.size() == 3);
    REQUIRE(graph.getShortestPathDistance(0, 4) == Approx(2.0));
}

TEMPLATE_TEST_CASE(
    "IGraph implementations return the shortest path distance derived from the path",
    "[graph][common][shortest-distance]",
    GraphADJList,
    GraphCSR
) {
    auto graph = makeWeightedPathGraph4<TestType>();

    REQUIRE(graph.getShortestPathDistance(0, 3) == Approx(9.0));
}

TEMPLATE_TEST_CASE(
    "IGraph implementations support custom edge-distance vectors for shortest paths",
    "[graph][common][custom-distances]",
    GraphADJList,
    GraphCSR
) {
    TestType graph(3);
    graph.addEdge(0, 1, 10.0, 10.0);
    graph.addEdge(1, 2, 10.0, 10.0);
    graph.addEdge(0, 2, 10.0, 1.0);
    graph.finalize();

    std::vector<double> custom_distances(graph.getNumDirectedEdges(), 100.0);

    const int e01 = graph.getEdgeId(0, 1);
    const int e12 = graph.getEdgeId(1, 2);
    const int e02 = graph.getEdgeId(0, 2);

    REQUIRE(e01 != INVALID_EDGE_ID);
    REQUIRE(e12 != INVALID_EDGE_ID);
    REQUIRE(e02 != INVALID_EDGE_ID);

    custom_distances[e01] = 1.0;
    custom_distances[e12] = 1.0;
    custom_distances[e02] = 10.0;

    auto path = graph.getShortestPath(0, 2, custom_distances);

    REQUIRE(path == std::vector<int>{0, 1, 2});
}

TEMPLATE_TEST_CASE(
    "IGraph implementations expose edge ids, endpoints and anti-edges consistently",
    "[graph][common][edge-id]",
    GraphADJList,
    GraphCSR
) {
    auto graph = makePathGraph4<TestType>();

    const int e01 = graph.getEdgeId(0, 1);
    const int e10 = graph.getEdgeId(1, 0);

    REQUIRE(e01 != INVALID_EDGE_ID);
    REQUIRE(e10 != INVALID_EDGE_ID);
    REQUIRE(e01 != e10);

    auto [u, v] = graph.getEdgeEndpoints(e01);
    REQUIRE(u == 0);
    REQUIRE(v == 1);

    const int anti = graph.getAntiEdge(e01);
    REQUIRE(anti == e10);

    auto [anti_u, anti_v] = graph.getEdgeEndpoints(anti);
    REQUIRE(anti_u == 1);
    REQUIRE(anti_v == 0);
}

TEMPLATE_TEST_CASE(
    "IGraph implementations expose capacities and distances by edge id",
    "[graph][common][edge-access]",
    GraphADJList,
    GraphCSR
) {
    auto graph = makePathGraph4<TestType>();

    const int e01 = graph.getEdgeId(0, 1);
    REQUIRE(e01 != INVALID_EDGE_ID);

    REQUIRE(graph.getEdgeCapacity(e01) == Approx(10.0));
    REQUIRE(graph.getEdgeDistance(e01) == Approx(1.0));
}

TEMPLATE_TEST_CASE(
    "IGraph implementations update edge distances by node pair and by edge id",
    "[graph][common][update-distance]",
    GraphADJList,
    GraphCSR
) {
    auto graph = makePathGraph4<TestType>();

    REQUIRE(graph.updateEdgeDistance(0, 1, 3.0));
    REQUIRE(graph.getEdgeDistance(0, 1) == Approx(3.0));

    const int e12 = graph.getEdgeId(1, 2);
    REQUIRE(e12 != INVALID_EDGE_ID);

    REQUIRE(graph.updateEdgeDistance(e12, 4.0));
    REQUIRE(graph.getEdgeDistance(e12) == Approx(4.0));
    REQUIRE(graph.getEdgeDistance(1, 2) == Approx(4.0));
}

TEMPLATE_TEST_CASE(
    "IGraph implementations reset all edge distances to one",
    "[graph][common][reset-distance]",
    GraphADJList,
    GraphCSR
) {
    auto graph = makeWeightedPathGraph4<TestType>();

    REQUIRE(graph.updateEdgeDistance(0, 1, 9.0));
    REQUIRE(graph.updateEdgeDistance(1, 2, 8.0));

    graph.resetEdgeDistance();

    for (int u = 0; u < graph.getNumNodes(); ++u) {
        for (int v : graph.neighbors(u)) {
            REQUIRE(graph.getEdgeDistance(u, v) == Approx(1.0));
        }
    }
}

TEMPLATE_TEST_CASE(
    "IGraph implementations compute exact diameter from shortest-path distances",
    "[graph][common][diameter]",
    GraphADJList,
    GraphCSR
) {
    auto graph = makeWeightedTriangleGraph<TestType>();

    REQUIRE(graph.getDiameter() == Approx(3.0));
}

TEMPLATE_TEST_CASE(
    "IGraph implementations return a usable approximate diameter",
    "[graph][common][diameter-approx]",
    GraphADJList,
    GraphCSR
) {
    auto graph = makeWeightedTriangleGraph<TestType>();

    const double approx_diameter = graph.getDiameterApprox();

    REQUIRE(std::isfinite(approx_diameter));
    REQUIRE(approx_diameter >= 0.0);
}