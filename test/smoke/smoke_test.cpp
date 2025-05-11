//
// Created by Mert Biyikli on 11.05.25.
//

#include "gtest/gtest.h"
#include "../../src/graph.h"
#include <stdexcept>



TEST(GraphTest, NodeAndEdgeExistence) {
    Graph g(2);

    EXPECT_EQ(g.numNodes(), 2);
}

TEST(GraphTest, Neighbours) {
    Graph g(2);

    g.addEdge(0, 1, 2.5);
    EXPECT_EQ(g.numEdges(), 1);

    ASSERT_EQ(g.neighbors(0).size(), 1);
    EXPECT_EQ(g.neighbors(0)[0].target, 1);
    EXPECT_DOUBLE_EQ(g.neighbors(0)[0].capacity, 2.5);

    ASSERT_EQ(g.neighbors(1).size(), 1);
    EXPECT_EQ(g.neighbors(1)[0].target, 0);
    EXPECT_DOUBLE_EQ(g.neighbors(1)[0].capacity, 2.5);
}

TEST(GraphTest, InvalidEdgeThrows) {
    Graph g(3);


    EXPECT_THROW(g.addEdge(0, 5, 3.0), std::out_of_range);  // node 5 doesn't exist
    EXPECT_THROW(g.addEdge(-1, 0, 1.0), std::out_of_range); // invalid index
}
