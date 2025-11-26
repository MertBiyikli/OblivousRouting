//
// Created by Mert Biyikli on 17.11.25.
//


#include <gtest/gtest.h>
#include "../../src/datastructures/graph_csr.h"

TEST(TinyGraph, Line4BasicGraphSanity) {
    Graph_csr g(4);
    g.addEdge(0,1,1.0);
    g.addEdge(1,2,1.0);
    g.addEdge(2,3,1.0);

    EXPECT_EQ(g.getNumNodes(), 4);
    EXPECT_EQ(g.getNumEdges(), 3);
}

TEST(TinyGraph, Line4GraphNeighbors) {
    Graph_csr g(4);
    g.addEdge(0,1,1.0);
    g.addEdge(1,2,1.0);
    g.addEdge(2,3,1.0);
    g.finalize();

    // Node 1: expects neighbors 0 and 2 (CSR sorted order)
    std::vector<int> nbr1;
    for (auto [u,v] : g.neighbors(1))
        nbr1.push_back(v);

    ASSERT_EQ(nbr1.size(), 2);
    EXPECT_EQ(nbr1[0], 0);
    EXPECT_EQ(nbr1[1], 2);

    // Node 0: expects [1]
    std::vector<int> nbr0;
    for (auto [u,v] : g.neighbors(0))
        nbr0.push_back(v);

    ASSERT_EQ(nbr0.size(), 1);
    EXPECT_EQ(nbr0[0], 1);

    // Node 3: expects [2]
    std::vector<int> nbr3;
    for (auto [u,v] : g.neighbors(3))
        nbr3.push_back(v);

    ASSERT_EQ(nbr3.size(), 1);
    EXPECT_EQ(nbr3[0], 2);
}


TEST(TinyGraph, Line4EdgeUpdates) {
    Graph_csr g(4);

    g.addEdge(0,1,1.0);
    g.addEdge(1,2,1.0);
    g.addEdge(2,3,1.0);

    g.finalize();


    // Update edge distance
    g.updateEdgeDistance(1,2,3.0); // should update the distance

    double dist = g.getEdgeDistance(1, 2);
    EXPECT_DOUBLE_EQ(dist, 3.0);
    double rev_dist = g.getEdgeDistance(2, 1);
    EXPECT_DOUBLE_EQ(rev_dist, 3.0);

}


TEST(TinyGraph, Line4StressTestUpdateMethod) {
    Graph_csr g(4);
    g.addEdge(0,1,1.0);
    g.addEdge(1,2,1.0);
    g.addEdge(2,3,1.0);

    g.finalize();

    double new_distance = 2.0;

    // update a non-existent edge
    EXPECT_THROW(g.updateEdgeDistance(0, 3, new_distance), std::runtime_error);


    // get the values of non-existent edges
    EXPECT_THROW(g.getEdgeDistance(0, 3), std::runtime_error);

    // update with out-of-range nodes
    EXPECT_THROW(g.updateEdgeDistance(1, 4, new_distance), std::out_of_range);


    // update with exising nodes but non existing edge
    EXPECT_THROW(g.updateEdgeDistance(1, 3, new_distance), std::runtime_error);
}


TEST(TinyGraph, Line4ShortestPath) {
    Graph_csr g(4);
    g.addEdge(0,1,1.0);
    g.addEdge(1,2,1.0);
    g.addEdge(2,3,1.0);

    g.finalize();

    auto path_0_3 = g.getShortestPath(0, 3);
    EXPECT_EQ(path_0_3.size(), 4);
    EXPECT_EQ(path_0_3[0], 0);
    EXPECT_EQ(path_0_3[1], 1);
    EXPECT_EQ(path_0_3[2], 2);
    EXPECT_EQ(path_0_3[3], 3);

    auto path_1_2 = g.getShortestPath(1, 2);
    EXPECT_EQ(path_1_2.size(), 2);
    EXPECT_EQ(path_1_2[0], 1);
    EXPECT_EQ(path_1_2[1], 2);
}

TEST(TinyGraph, Line4Diameter) {
    Graph_csr g(4);
    g.addEdge(0,1,1.0);
    g.addEdge(1,2,1.0);
    g.addEdge(2,3,1.0);

    g.finalize();
    EXPECT_EQ(g.GetDiameter(), 3); // distance matrix not initialized

}