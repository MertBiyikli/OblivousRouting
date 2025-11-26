//
// Created by Mert Biyikli on 25.11.25.
//

#include "graph_base.h"

// --- Test 1: Node count --------------------------------------------------
TYPED_TEST_P(IGraphTest, NodeCount) {
    this->G = TypeParam(5);
    EXPECT_EQ(this->G.n, 5);
    EXPECT_EQ(this->G.vertices.size(), 5);
}

// --- Test 2: Add edges & check neighbors --------------------------------
TYPED_TEST_P(IGraphTest, AddEdgesAndNeighbors) {
    this->buildSimpleTriangle();

    auto N0 = this->G.neighbors(0);
    auto N1 = this->G.neighbors(1);
    auto N2 = this->G.neighbors(2);

    // All nodes have 2 neighbors in a triangle
    EXPECT_EQ(N0.size(), 2);
    EXPECT_EQ(N1.size(), 2);
    EXPECT_EQ(N2.size(), 2);

    // Check expected connectivity using a set
    std::vector<int> nbr0(N0.begin(), N0.end());
    std::vector<int> nbr1(N1.begin(), N1.end());
    std::vector<int> nbr2(N2.begin(), N2.end());

    EXPECT_TRUE((std::set<int>{nbr0.begin(), nbr0.end()} == std::set<int>{1,2}));
    EXPECT_TRUE((std::set<int>{nbr1.begin(), nbr1.end()} == std::set<int>{0,2}));
    EXPECT_TRUE((std::set<int>{nbr2.begin(), nbr2.end()} == std::set<int>{0,1}));
}

// --- Test 3: Edge capacity / distance ------------------------------------
TYPED_TEST_P(IGraphTest, EdgeCapacityDistance) {
    this->buildSimpleTriangle();

    // Check a few edges
    EXPECT_DOUBLE_EQ(this->G.getEdgeCapacity(0,1), 1.0);
    EXPECT_DOUBLE_EQ(this->G.getEdgeDistance(0,1), 1.0);

    // Update and re-check
    this->G.updateEdgeDistance(1,2, 5.0);
    EXPECT_DOUBLE_EQ(this->G.getEdgeDistance(1,2), 5.0);
}

// --- Test 4: Edge-index lookup -------------------------------------------
TYPED_TEST_P(IGraphTest, EdgeIndexLookup) {
    this->buildSimpleTriangle();

    int e01 = this->G.getEdgeId(0,1);
    int e12 = this->G.getEdgeId(1,2);
    int e20 = this->G.getEdgeId(2,0);

    EXPECT_GE(e01, 0);
    EXPECT_GE(e12, 0);
    EXPECT_GE(e20, 0);

    EXPECT_DOUBLE_EQ(this->G.getEdgeCapacity(e01), 1.0);
    EXPECT_DOUBLE_EQ(this->G.getEdgeCapacity(e12), 1.0);
    EXPECT_DOUBLE_EQ(this->G.getEdgeCapacity(e20), 1.0);
}

// --- register all tests --------------------------------------------------
REGISTER_TYPED_TEST_SUITE_P(
    IGraphTest,
    NodeCount,
    AddEdgesAndNeighbors,
    EdgeCapacityDistance,
    EdgeIndexLookup
);
