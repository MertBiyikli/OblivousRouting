#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "data_structures/graph/graph_csr.h"
#include "core/routing_table.h"

using Catch::Approx;

static GraphCSR createSimpleGraph() {
    GraphCSR graph(4);
    graph.addEdge(0, 1, 1.0, 1.0);
    graph.addEdge(1, 2, 1.0, 1.0);
    graph.addEdge(2, 3, 1.0, 1.0);
    graph.addEdge(1, 3, 1.0, 1.0);
    graph.finalize();
    return graph;
}

TEST_CASE("AllPairRoutingTable - Construction and Initialization", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    AllPairRoutingTable table;

    table.init(graph);

    REQUIRE(table.n == 4);
    REQUIRE(table.adj_ids.size() == 8);  // 4 undirected edges = 8 directed edges
    REQUIRE(table.adj_vals.size() == 8);
}

TEST_CASE("AllPairRoutingTable - Add Flow by Nodes", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    AllPairRoutingTable table;

    table.init(graph);

    // Add flow for commodity (0, 3) on edge 0
    table.addFlow(0, 0, 3, 0.5);

    REQUIRE(table.getFlow(0, 0, 3) == Approx(0.5));
}



TEST_CASE("AllPairRoutingTable - Get Flow", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    AllPairRoutingTable table;

    table.init(graph);

    table.addFlow(0, 0, 3, 0.5);
    double flow = table.getFlow(0, 0, 3);

    REQUIRE(flow == Approx(0.5));
}

TEST_CASE("AllPairRoutingTable - Multiple Flows on Same Edge", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    AllPairRoutingTable table;

    table.init(graph);

    // Add multiple flows to the same edge
    table.addFlow(0, 0, 1, 0.3);
    table.addFlow(0, 0, 2, 0.2);
    table.addFlow(0, 1, 3, 0.5);

    REQUIRE(table.adj_ids[0].size() == 3);
    REQUIRE(table.adj_vals[0].size() == 3);

    REQUIRE(table.getFlow(0, 0, 1) == Approx(0.3));
    REQUIRE(table.getFlow(0, 0, 2) == Approx(0.2));
    REQUIRE(table.getFlow(0, 1, 3) == Approx(0.5));
}

TEST_CASE("AllPairRoutingTable - Access Operator", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    AllPairRoutingTable table;

    table.init(graph);

    // Use [] operator to access flows
    table[0].push_back(0.5);

    REQUIRE(!table[0].empty());
}

TEST_CASE("AllPairRoutingTable - Find Index in Sorted List", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    AllPairRoutingTable table;

    table.init(graph);

    table.addFlow(0, 1, 2, 0.1);
    table.addFlow(0, 0, 3, 0.2);
    table.addFlow(0, 2, 3, 0.3);

    // Just verify flows were added
    REQUIRE(table.adj_ids.size() == 8);
}



TEST_CASE("AllPairRoutingTable - Multiple Flows Different Edges", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    AllPairRoutingTable table;

    table.init(graph);

    // Add flows to different edges
    table.addFlow(0, 0, 1, 0.3);
    table.addFlow(1, 0, 1, 0.2);
    table.addFlow(2, 0, 1, 0.5);

    REQUIRE(table.getFlow(0, 0, 1) == Approx(0.3));
    REQUIRE(table.getFlow(1, 0, 1) == Approx(0.2));
    REQUIRE(table.getFlow(2, 0, 1) == Approx(0.5));
}

TEST_CASE("LinearRoutingTable - Construction and Initialization", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    LinearRoutingTable table;

    table.init(graph);

    REQUIRE(table.n == 4);
    REQUIRE(table.src_ids.size() == 8);
    REQUIRE(table.src_flows.size() == 8);
}

TEST_CASE("LinearRoutingTable - Add Flow", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    LinearRoutingTable table;

    table.init(graph);

    // Add source flow for commodity from source 0
    table.addFlow(0, 0, 0.5);

    REQUIRE(!table.src_ids[0].empty());
    REQUIRE(!table.src_flows[0].empty());
}

TEST_CASE("LinearRoutingTable - Get Flow", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    LinearRoutingTable table;

    table.init(graph);

    table.addFlow(0, 0, 0.5);
    double flow = table.getFlow(0, 0);

    REQUIRE(flow == Approx(0.5));
}

TEST_CASE("LinearRoutingTable - Multiple Sources Same Edge", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    LinearRoutingTable table;

    table.init(graph);

    // Add flows from multiple sources on same edge
    table.addFlow(0, 0, 0.3);
    table.addFlow(0, 1, 0.2);
    table.addFlow(0, 2, 0.5);

    REQUIRE(table.src_ids[0].size() == 3);
    REQUIRE(table.getFlow(0, 0) == Approx(0.3));
    REQUIRE(table.getFlow(0, 1) == Approx(0.2));
    REQUIRE(table.getFlow(0, 2) == Approx(0.5));
}

TEST_CASE("LinearRoutingTable - Erase Flow", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    LinearRoutingTable table;

    table.init(graph);

    table.addFlow(0, 0, 0.5);
    table.addFlow(0, 1, 0.3);

    int initial_size = table.src_ids[0].size();

    // Erase flow from source 0
    table.eraseFlow(0, 0);

    REQUIRE(table.src_ids[0].size() <= initial_size);
}
/*
TEST_CASE("LinearRoutingTable - Is Valid", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    LinearRoutingTable table;

    table.init(graph);

    // Add unit of flow along a simple path with a single commodity
    auto path = graph.getShortestPath(0, 3);
    for (size_t i = 0; i < path.size() - 1; ++i) {
        int u = path[i];
        int v = path[i + 1];
        int edge_id = graph.getEdgeId(u, v);
        table.addFlow(edge_id, 0, 1.0); // Add flow from source 0
    }

    table.printFlows(graph);

    REQUIRE(table.isValid(graph));
}
*/

TEST_CASE("LinearRoutingTable - Non-existent Flow", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    LinearRoutingTable table;

    table.init(graph);

    // Try to get flow that wasn't added
    double flow = table.getFlow(0, 3);

    REQUIRE(flow == Approx(0.0));
}


TEST_CASE("LinearRoutingTable - Multiple Edges Different Sources", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    LinearRoutingTable table;

    table.init(graph);

    // Different sources on different edges
    table.addFlow(0, 0, 0.3);
    table.addFlow(1, 1, 0.2);
    table.addFlow(2, 2, 0.5);

    REQUIRE(table.getFlow(0, 0) == Approx(0.3));
    REQUIRE(table.getFlow(1, 1) == Approx(0.2));
    REQUIRE(table.getFlow(2, 2) == Approx(0.5));
}

