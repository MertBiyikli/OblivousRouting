//
// Created by Mert Biyikli on 11.05.25.
//

#include "gtest/gtest.h"
#include "../../src/graph.h"
#include <stdexcept>



TEST(GraphTest, NodeAddition) {
    Graph g;
    EXPECT_EQ(g.numNodes(), 0);

    g.addNode();
    g.addNode();
    EXPECT_EQ(g.numNodes(), 2);
}

TEST(GraphTest, ValidEdgeAddition) {
    Graph g;
    g.addNode(); // node 0
    g.addNode(); // node 1

    g.addEdge(0, 1, 2.5);
    EXPECT_EQ(g.numEdges(), 1);

    const auto& neighbors = g.neighbors(0);
    ASSERT_EQ(neighbors.size(), 1);
    EXPECT_EQ(neighbors[0].target, 1);
    EXPECT_DOUBLE_EQ(neighbors[0].capacity, 2.5);
}

TEST(GraphTest, InvalidEdgeThrows) {
    Graph g;
    g.addNode();

    EXPECT_THROW(g.addEdge(0, 1, 3.0), std::out_of_range);  // node 1 doesn't exist
    EXPECT_THROW(g.addEdge(-1, 0, 1.0), std::out_of_range); // invalid index
    EXPECT_THROW(g.addEdge(1, 0, 1.0), std::out_of_range);  // source doesn't exist
}
