#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <vector>

#include "graph_test_helpers.h"
#include "data_structures/graph/graph_csr.h"

using Catch::Approx;
// Only GraphCSR-specific tests belong here.
// Generic IGraph behavior is tested once in test_graph_common.cpp.

TEST_CASE("GraphCSR finalization builds CSR arrays and clears temporary edge storage", "[graph][csr][structure]") {
    GraphCSR graph(4);
    graph.addEdge(0, 1, 10.0, 1.0);
    graph.addEdge(0, 2, 20.0, 2.0);
    graph.addEdge(1, 3, 30.0, 3.0);

    REQUIRE(graph.tmp_edges.size() == 6);

    graph.finalize();

    REQUIRE(graph.head.size() == 5);
    REQUIRE(graph.from.size() == 6);
    REQUIRE(graph.to.size() == 6);
    REQUIRE(graph.capacity.size() == 6);
    REQUIRE(graph.distance.size() == 6);

    REQUIRE(graph.tmp_edges.empty());
    REQUIRE(graph.tmp_capacity.empty());
    REQUIRE(graph.tmp_distance.empty());
}

TEST_CASE("GraphCSR bidirectional Dijkstra returns the same path length as ordinary Dijkstra", "[graph][csr][bidirectional-dijkstra]") {
    auto graph = makePathGraph4<GraphCSR>();

    auto ordinary_path = graph.getShortestPath(0, 3);
    auto bidirectional_path = graph.getShortestPathBidirectionalSearch(0, 3);

    REQUIRE(bidirectional_path == ordinary_path);
    REQUIRE(bidirectional_path == std::vector<int>{0, 1, 2, 3});
}

TEST_CASE("GraphCSR bidirectional Dijkstra respects custom edge-distance vectors", "[graph][csr][bidirectional-dijkstra][custom-distances]") {
    GraphCSR graph(3);
    graph.addEdge(0, 1, 10.0, 10.0);
    graph.addEdge(1, 2, 10.0, 10.0);
    graph.addEdge(0, 2, 10.0, 1.0);
    graph.finalize();

    std::vector<double> custom_distances(graph.getNumDirectedEdges(), 100.0);
    custom_distances[graph.getEdgeId(0, 1)] = 1.0;
    custom_distances[graph.getEdgeId(1, 2)] = 1.0;
    custom_distances[graph.getEdgeId(0, 2)] = 10.0;

    auto path = graph.getShortestPathBidirectionalSearch(0, 2, custom_distances);

    REQUIRE(path == std::vector<int>{0, 1, 2});
}

TEST_CASE("GraphCSR exposes path edge ids for a shortest path", "[graph][csr][path-edges]") {
    auto graph = makePathGraph4<GraphCSR>();

    auto path_edges = graph.getPathEdges(0, 3);

    REQUIRE(path_edges.size() == 3);
    REQUIRE(path_edges[0] == graph.getEdgeId(0, 1));
    REQUIRE(path_edges[1] == graph.getEdgeId(1, 2));
    REQUIRE(path_edges[2] == graph.getEdgeId(2, 3));
}
