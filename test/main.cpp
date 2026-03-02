#include <gtest/gtest.h>
#include <limits>
#include "../src/datastructures/GraphADJList.h"
#include "../src/datastructures/GraphADJMatrix.h"
#include "../src/utils/my_math.h"
#include "../src/utils/hash.h"
#include "../src/solver/routing_table.h"


template<typename T>
T makeTriangle() {
    T g(3);
    g.addEdge(0, 1, 1.0, 1.0);
    g.addEdge(1, 2, 1.0, 1.0);
    g.addEdge(0, 2, 1.0, 1.0);
    g.finalize();
    return g;
}

// ---------------------------------------------------------------------------
// Test 1 – addEdge / getNumNodes / getNumEdges (ADJList)
// ---------------------------------------------------------------------------
TEST(GraphADJList, AddEdgeAndCounts) {
    GraphADJList g(4);
    g.addEdge(0, 1, 2.0, 1.0);
    g.addEdge(1, 2, 3.0, 1.0);
    g.addEdge(2, 3, 4.0, 1.0);
    g.finalize();

    EXPECT_EQ(g.getNumNodes(), 4);
    // 3 undirected edges → 6 directed edges
    EXPECT_EQ(g.getNumDirectedEdges(), 6);
}

// ---------------------------------------------------------------------------
// Test 2 – getEdgeId / edgeEndpoints round-trip (ADJList)
// ---------------------------------------------------------------------------
TEST(GraphADJList, EdgeIdRoundTrip) {
    auto g = makeTriangle<GraphADJList>();

    for (int u = 0; u < 3; ++u) {
        for (int v : g.neighbors(u)) {
            int eid = g.getEdgeId(u, v);
            ASSERT_NE(eid, INVALID_EDGE_ID) << "edge (" << u << "," << v << ") has no id";
            auto [eu, ev] = g.edgeEndpoints(eid);
            EXPECT_EQ(eu, u);
            EXPECT_EQ(ev, v);
        }
    }
}

// ---------------------------------------------------------------------------
// Test 3 – getAntiEdge: anti-edge of (u→v) is (v→u) (ADJList)
// ---------------------------------------------------------------------------
TEST(GraphADJList, AntiEdge) {
    auto g = makeTriangle<GraphADJList>();

    int e01 = g.getEdgeId(0, 1);
    int e10 = g.getEdgeId(1, 0);
    ASSERT_NE(e01, INVALID_EDGE_ID);
    ASSERT_NE(e10, INVALID_EDGE_ID);
    EXPECT_EQ(g.getAntiEdge(e01), e10);
    EXPECT_EQ(g.getAntiEdge(e10), e01);
}

// ---------------------------------------------------------------------------
// Test 4 – getShortestPath: path on a line graph (ADJList)
// ---------------------------------------------------------------------------
TEST(GraphADJList, ShortestPathLine) {
    GraphADJList g(5);
    for (int i = 0; i < 4; ++i)
        g.addEdge(i, i + 1, 1.0, 1.0);
    g.finalize();

    auto path = g.getShortestPath(0, 4);
    ASSERT_EQ(path.size(), 5u);
    for (int i = 0; i < 5; ++i)
        EXPECT_EQ(path[i], i);
}

// ---------------------------------------------------------------------------
// Test 5 – updateEdgeDistance / getEdgeDistance (ADJList)
// ---------------------------------------------------------------------------
TEST(GraphADJList, UpdateEdgeDistance) {
    auto g = makeTriangle<GraphADJList>();

    EXPECT_DOUBLE_EQ(g.getEdgeDistance(0, 1), 1.0);
    bool ok = g.updateEdgeDistance(0, 1, 5.5);
    EXPECT_TRUE(ok);
    EXPECT_DOUBLE_EQ(g.getEdgeDistance(0, 1), 5.5);
}

// ---------------------------------------------------------------------------
// Test 6 – isConnected / disconnected graph (ADJList)
// ---------------------------------------------------------------------------
TEST(GraphADJList, ConnectivityCheck) {
    GraphADJList connected(3);
    connected.addEdge(0, 1, 1.0, 1.0);
    connected.addEdge(1, 2, 1.0, 1.0);
    connected.finalize();
    EXPECT_TRUE(connected.isConnected());

    GraphADJList disconnected(4);
    disconnected.addEdge(0, 1, 1.0, 1.0);
    // nodes 2 and 3 are isolated
    disconnected.finalize();
    EXPECT_FALSE(disconnected.isConnected());
}

// ---------------------------------------------------------------------------
// Test 7 – GraphADJMatrix: addEdge / getEdgeCapacity / getEdgeDistance
// ---------------------------------------------------------------------------
TEST(GraphADJMatrix, EdgeAttributes) {
    GraphADJMatrix g(3);
    g.addEdge(0, 1, 2.5, 3.0);
    g.finalize();

    EXPECT_DOUBLE_EQ(g.getEdgeCapacity(0, 1), 2.5);
    EXPECT_DOUBLE_EQ(g.getEdgeDistance(0, 1), 3.0);
    // undirected: reverse direction should be the same
    EXPECT_DOUBLE_EQ(g.getEdgeCapacity(1, 0), 2.5);
    EXPECT_DOUBLE_EQ(g.getEdgeDistance(1, 0), 3.0);
    // non-existent edge
    EXPECT_DOUBLE_EQ(g.getEdgeCapacity(0, 2), 0.0);
}

// ---------------------------------------------------------------------------
// Test 8 – GraphADJMatrix: edge-id based accessors match node-pair accessors
// ---------------------------------------------------------------------------
TEST(GraphADJMatrix, EdgeIdAccessors) {
    auto g = makeTriangle<GraphADJMatrix>();

    for (int u = 0; u < 3; ++u) {
        for (int v : g.neighbors(u)) {
            int eid = g.getEdgeId(u, v);
            ASSERT_NE(eid, INVALID_EDGE_ID);
            EXPECT_DOUBLE_EQ(g.getEdgeCapacity(eid), g.getEdgeCapacity(u, v));
            EXPECT_DOUBLE_EQ(g.getEdgeDistance(eid), g.getEdgeDistance(u, v));
        }
    }
}

// ---------------------------------------------------------------------------
// Test 9 – GraphADJMatrix: shortest path on a weighted graph
//          0 --10-- 1 --1-- 2
//          |                |
//          +------2---------+
//   Expected shortest path 0→2: [0, 2] with distance 2
// ---------------------------------------------------------------------------
TEST(GraphADJMatrix, ShortestPathWeighted) {
    GraphADJMatrix g(3);
    g.addEdge(0, 1, 1.0, 10.0);
    g.addEdge(1, 2, 1.0, 1.0);
    g.addEdge(0, 2, 1.0, 2.0);
    g.finalize();

    auto path = g.getShortestPath(0, 2);
    ASSERT_FALSE(path.empty());
    EXPECT_EQ(path.front(), 0);
    EXPECT_EQ(path.back(), 2);

    double dist = g.getShortestPathDistance(0, 2);
    EXPECT_DOUBLE_EQ(dist, 2.0);
}

// ---------------------------------------------------------------------------
// Test 10 – NegativeExponent utility
// ---------------------------------------------------------------------------
TEST(NegativeExponent, CorrectValues) {
    EXPECT_DOUBLE_EQ(NegativeExponent(2.0, 3), 1.0 / 8.0);
    EXPECT_DOUBLE_EQ(NegativeExponent(10.0, 0), 1.0);
    EXPECT_DOUBLE_EQ(NegativeExponent(5.0, 1), 0.2);

    EXPECT_THROW(NegativeExponent(0.0, 2), std::invalid_argument);
    EXPECT_THROW(NegativeExponent(2.0, -1), std::invalid_argument);
}

// ===========================================================================
//  MWU Subroutine Tests
//  All tests below use the triangle graph (3 nodes, edges 0-1, 1-2, 0-2)
//  with unit capacity/distance and operate directly on the routing-table
//  data structures that are the core of the MWU algorithm.
// ===========================================================================

// ---------------------------------------------------------------------------
// Helper: build a finalized triangle GraphADJList (reused across MWU tests)
// ---------------------------------------------------------------------------
static GraphADJList makeTriangleList() {
    GraphADJList g(3);
    g.addEdge(0, 1, 1.0, 1.0);
    g.addEdge(1, 2, 1.0, 1.0);
    g.addEdge(0, 2, 1.0, 1.0);
    g.finalize();
    return g;
}

// ---------------------------------------------------------------------------
// MWU Test 1 – getCommodityID: valid and invalid pairs
//   The function encodes an (s,t) pair with s<t as s*n+t.
// ---------------------------------------------------------------------------
TEST(MWU_CommodityID, ValidAndInvalid) {
    const int n = 4;
    // s < t → valid
    EXPECT_EQ(getCommodityID(n, 0, 1), 0 * n + 1);
    EXPECT_EQ(getCommodityID(n, 0, 3), 0 * n + 3);
    EXPECT_EQ(getCommodityID(n, 2, 3), 2 * n + 3);

    // s >= t → invalid
    EXPECT_EQ(getCommodityID(n, 1, 0), INVALID_COMMODITY_ID);
    EXPECT_EQ(getCommodityID(n, 2, 2), INVALID_COMMODITY_ID);
}

// ---------------------------------------------------------------------------
// MWU Test 2 – AllPairRoutingTable::init: sizes match graph
// ---------------------------------------------------------------------------
TEST(MWU_AllPairRoutingTable, InitSizesMatchGraph) {
    auto g = makeTriangleList();
    AllPairRoutingTable table;
    table.init(g);

    EXPECT_EQ((int)table.adj_ids.size(),  g.getNumEdges());
    EXPECT_EQ((int)table.adj_vals.size(), g.getNumEdges());
    EXPECT_EQ((int)table.anti_edge.size(), g.getNumEdges());

    // Every anti-edge must be valid for a fully-connected undirected graph
    for (int e = 0; e < g.getNumEdges(); ++e) {
        EXPECT_NE(table.anti_edge[e], INVALID_EDGE_ID);
    }
}

// ---------------------------------------------------------------------------
// MWU Test 3 – AllPairRoutingTable::addFlow / getFlow: basic round-trip
// ---------------------------------------------------------------------------
TEST(MWU_AllPairRoutingTable, AddAndGetFlow) {
    auto g = makeTriangleList();
    AllPairRoutingTable table;
    table.init(g);

    // Route commodity (0→2) with fraction 0.5 on edge (0,2)
    int e02 = g.getEdgeId(0, 2);
    ASSERT_NE(e02, INVALID_EDGE_ID);

    table.addFlow(e02, 0, 2, 0.5);
    EXPECT_DOUBLE_EQ(table.getFlow(e02, 0, 2), 0.5);

    // Querying the opposite commodity direction should return −flow
    EXPECT_DOUBLE_EQ(table.getFlow(e02, 2, 0), -0.5);

    // Unrelated commodity should return 0
    EXPECT_DOUBLE_EQ(table.getFlow(e02, 0, 1), 0.0);
}

// ---------------------------------------------------------------------------
// MWU Test 4 – AllPairRoutingTable::cancel2Cycle
//   Add flow on e and the same amount on anti(e); they should cancel each other.
// ---------------------------------------------------------------------------
TEST(MWU_AllPairRoutingTable, Cancel2Cycle) {
    auto g = makeTriangleList();
    AllPairRoutingTable table;
    table.init(g);

    int e01 = g.getEdgeId(0, 1);
    int e10 = g.getEdgeId(1, 0);
    ASSERT_NE(e01, INVALID_EDGE_ID);
    ASSERT_NE(e10, INVALID_EDGE_ID);

    // First add 0.7 on the anti-edge (1→0) for commodity (0,1)
    table.addFlow(e10, 0, 1, 0.7);
    EXPECT_DOUBLE_EQ(table.getFlow(e10, 0, 1), 0.7);

    // Now add 0.7 on e01 for the same commodity — the cancel2Cycle routine
    // inside addFlow should wipe out the anti-edge flow completely
    table.addFlow(e01, 0, 1, 0.7);

    EXPECT_NEAR(table.getFlow(e10, 0, 1), 0.0, SOFT_EPS);
    // Net result: no flow on e01 either (they fully cancelled)
    EXPECT_NEAR(table.getFlow(e01, 0, 1), 0.0, SOFT_EPS);
}

// ---------------------------------------------------------------------------
// MWU Test 5 – AllPairRoutingTable: ids stay sorted after multiple insertions
// ---------------------------------------------------------------------------
TEST(MWU_AllPairRoutingTable, SortedInvariant) {
    auto g = makeTriangleList();
    AllPairRoutingTable table;
    table.init(g);

    int e = g.getEdgeId(0, 1);
    ASSERT_NE(e, INVALID_EDGE_ID);

    // Insert three commodities in non-ascending order of commodity_id
    table.addFlow(e, 0, 2, 0.3);  // commodity 0*3+2=2
    table.addFlow(e, 0, 1, 0.4);  // commodity 0*3+1=1
    table.addFlow(e, 1, 2, 0.1);  // commodity 1*3+2=5

    const auto& ids = table.adj_ids[e];
    for (size_t i = 1; i < ids.size(); ++i) {
        EXPECT_LT(ids[i - 1], ids[i]) << "ids not sorted at position " << i;
    }
}

// ---------------------------------------------------------------------------
// MWU Test 6 – LinearRoutingTable::init and addFlow / getFlow
// ---------------------------------------------------------------------------
TEST(MWU_LinearRoutingTable, AddAndGetFlow) {
    auto g = makeTriangleList();
    LinearRoutingTable table;
    table.init(g);

    int e01 = g.getEdgeId(0, 1);
    ASSERT_NE(e01, INVALID_EDGE_ID);

    table.addFlow(e01, 0, 0.6);
    EXPECT_DOUBLE_EQ(table.getFlow(e01, 0), 0.6);

    // Accumulation: add more flow for the same source
    table.addFlow(e01, 0, 0.2);
    EXPECT_DOUBLE_EQ(table.getFlow(e01, 0), 0.8);

    // Different source should be independent
    table.addFlow(e01, 2, 0.3);
    EXPECT_DOUBLE_EQ(table.getFlow(e01, 2), 0.3);
    EXPECT_DOUBLE_EQ(table.getFlow(e01, 0), 0.8); // unchanged
}

// ---------------------------------------------------------------------------
// MWU Test 7 – LinearRoutingTable: src_ids stay sorted after insertions
// ---------------------------------------------------------------------------
TEST(MWU_LinearRoutingTable, SortedInvariant) {
    auto g = makeTriangleList();
    LinearRoutingTable table;
    table.init(g);

    int e = g.getEdgeId(0, 2);
    ASSERT_NE(e, INVALID_EDGE_ID);

    // Insert in reverse order
    table.addFlow(e, 2, 0.5);
    table.addFlow(e, 1, 0.3);
    table.addFlow(e, 0, 0.2);

    const auto& ids = table.src_ids[e];
    for (size_t i = 1; i < ids.size(); ++i) {
        EXPECT_LT(ids[i - 1], ids[i]) << "src_ids not sorted at position " << i;
    }
}

// ---------------------------------------------------------------------------
// MWU Test 8 – LinearRoutingScheme::addFlow / getFlow: linearity w.r.t root
//   We manually set up the linear routing table for a path 0-1-2 with root 0,
//   then verify the LinearRoutingScheme::getFlow(e, s, t) formula:
//     getFlow(e, s, t) = flow_sx - flow_tx  (accounting for anti-edges)
// ---------------------------------------------------------------------------
TEST(MWU_LinearRoutingScheme, GetFlowLinearity) {
    auto g = makeTriangleList();

    LinearRoutingTable table;
    table.init(g);

    // For root x=0: route source 1 → root 0 through edge (1,0)
    // and source 2 → root 0 through path 2-0
    int e10 = g.getEdgeId(1, 0);
    int e20 = g.getEdgeId(2, 0);
    ASSERT_NE(e10, INVALID_EDGE_ID);
    ASSERT_NE(e20, INVALID_EDGE_ID);

    table.addFlow(e10, /*s=*/1, 1.0);   // f_{e10}(1) = 1
    table.addFlow(e20, /*s=*/2, 1.0);   // f_{e20}(2) = 1

    LinearRoutingScheme scheme(g, /*root=*/0, std::move(table));

    // Commodity (1→2): flow on e10 = f(e10,1) - f(e10,2) = 1 - 0 = 1
    // minus contribution from anti-edge e01: f(e01,1)-f(e01,2) = 0-0=0
    // net = 1 - 0 = 1
    double flow_e10_12 = scheme.getFlow(e10, 1, 2);
    EXPECT_NEAR(std::abs(flow_e10_12), 1.0, SOFT_EPS);
}

// ---------------------------------------------------------------------------
// MWU Test 9 – AllPairRoutingScheme::routeDemands: congestion is non-negative
//   and symmetric (the same value is mirrored on the anti-edge).
//
//   We build a complete routing on the triangle:
//     commodity (0,1): all flow on edge (0,1)
//     commodity (0,2): all flow on edge (0,2)
//     commodity (1,2): all flow on edge (1,2)
//   Then route a unit demand for every (s,t) pair and check symmetry.
// ---------------------------------------------------------------------------
TEST(MWU_AllPairRoutingScheme, CongestionSymmetry) {
    auto g = makeTriangleList();
    AllPairRoutingTable table;
    table.init(g);

    // Add unit flow for each commodity on its direct edge
    auto addSymmetric = [&](int u, int v, int s, int t, double f) {
        int e = g.getEdgeId(u, v);
        if (e != INVALID_EDGE_ID) table.addFlow(e, s, t, f);
    };

    addSymmetric(0, 1, 0, 1, 1.0);
    addSymmetric(0, 2, 0, 2, 1.0);
    addSymmetric(1, 2, 1, 2, 1.0);

    AllPairRoutingScheme scheme(g, std::move(table));

    // Uniform unit demand for all pairs
    DemandMap demands;
    for (int s = 0; s < 3; ++s)
        for (int t = 0; t < 3; ++t)
            if (s != t) demands.addDemand(s, t, 1.0);

    std::vector<double> cong(g.getNumEdges(), 0.0);
    scheme.routeDemands(cong, demands);

    // All congestion values must be non-negative
    for (int e = 0; e < g.getNumEdges(); ++e) {
        EXPECT_GE(cong[e], 0.0) << "negative congestion on edge " << e;
    }

    // Congestion must be symmetric: cong[e] == cong[anti(e)]
    for (int e = 0; e < g.getNumEdges(); ++e) {
        int anti_e = g.getAntiEdge(e);
        EXPECT_NEAR(cong[e], cong[anti_e], SOFT_EPS)
            << "congestion not symmetric for edge " << e;
    }
}
