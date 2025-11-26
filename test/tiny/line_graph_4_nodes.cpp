//
// Created by Mert Biyikli on 14.11.25.
//

#include <gtest/gtest.h>
#include "../../src/datastructures/graph.h"
#include "../../src/tree_based/frt/raecke_frt_solver.h"

TEST(TinyGraph, Line4BasicGraphSanity) {
    Graph g(4);
    g.addEdge(0,1,1.0);
    g.addEdge(1,2,1.0);
    g.addEdge(2,3,1.0);

    EXPECT_EQ(g.getNumNodes(), 4);
    EXPECT_EQ(g.getNumEdges(), 3);
}

TEST(TinyGraph, Line4GraphNeighbors) {
    Graph g(4);
    g.addEdge(0,1,1.0);
    g.addEdge(1,2,1.0);
    g.addEdge(2,3,1.0);

    auto& neighbors_1 = g.neighbors(1);
    EXPECT_EQ(neighbors_1.size(), 2);
    EXPECT_EQ(neighbors_1[0], 0);
    EXPECT_EQ(neighbors_1[1], 2);

    auto neighbors_0 = g.neighbors(0);
    EXPECT_EQ(neighbors_0.size(), 1);
    EXPECT_EQ(neighbors_0[0], 1);

    auto neighbors_3 = g.neighbors(3);
    EXPECT_EQ(neighbors_3.size(), 1);
    EXPECT_EQ(neighbors_3[0], 2);
}


TEST(TinyGraph, Line4EdgeUpdates) {
    Graph g(4);

    g.addEdge(0,1,1.0);
    g.addEdge(1,2,1.0);
    g.addEdge(2,3,1.0);

    // Update edge capacity
    g.updateEdgeCapacity(1,2,2.5); // should update the capacity

    double cap = g.getEdgeCapacity(1, 2);
    EXPECT_DOUBLE_EQ(cap, 2.5);
    double rev_cap = g.getEdgeCapacity(2, 1);
    EXPECT_DOUBLE_EQ(rev_cap, 2.5);


    // Update edge distance
    g.updateEdgeDistance(1,2,3.0); // should update the distance

    double dist = g.getEdgeDistance(1, 2);
    EXPECT_DOUBLE_EQ(dist, 3.0);
    double rev_dist = g.getEdgeDistance(2, 1);
    EXPECT_DOUBLE_EQ(rev_dist, 3.0);

}


TEST(TinyGraph, Line4StressTestUpdateMethod) {
    Graph g(4);
    g.addEdge(0,1,1.0);
    g.addEdge(1,2,1.0);
    g.addEdge(2,3,1.0);

    double new_capacity = 5.0;
    double new_distance = 2.0;

    // update a non-existent edge
    EXPECT_THROW(g.updateEdgeCapacity(0, 3, new_capacity), std::runtime_error);
    EXPECT_THROW(g.updateEdgeDistance(0, 3, new_distance), std::runtime_error);


    // get the values of non-existent edges
    EXPECT_THROW(g.getEdgeCapacity(0, 3), std::runtime_error);
    EXPECT_THROW(g.getEdgeDistance(0, 3), std::runtime_error);

    // update with out-of-range nodes
    EXPECT_THROW(g.updateEdgeCapacity(5, 2, new_capacity), std::out_of_range);
    EXPECT_THROW(g.updateEdgeDistance(1, 4, new_distance), std::out_of_range);


    // update with exising nodes but non existing edge
    EXPECT_THROW(g.getEdgeCapacity(0, 2), std::runtime_error);
    EXPECT_THROW(g.updateEdgeDistance(1, 3, new_distance), std::runtime_error);
}


TEST(TinyGraph, Line4ShortestPath) {
    Graph g(4);
    g.addEdge(0,1,1.0);
    g.addEdge(1,2,1.0);
    g.addEdge(2,3,1.0);

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
    Graph g(4);
    g.addEdge(0,1,1.0, 2.0);
    g.addEdge(1,2,1.0, 3.0);
    g.addEdge(2,3,1.0, 4.0);

    EXPECT_EQ(g.IsDistanceMatrixComputed(), false);
    EXPECT_THROW(g.GetDiameter(), std::runtime_error); // distance matrix not initialized

    g.createDistanceMatrix();
    EXPECT_EQ(g.IsDistanceMatrixComputed(), true);

    double diameter = g.GetDiameter();
    EXPECT_DOUBLE_EQ(diameter, 9.0); // 2 + 3 + 4 = 9

    // check shortest path distance
    double dist_0_2 = g.getShortestDistance(0, 2);
    EXPECT_DOUBLE_EQ(dist_0_2, 5.0); // 2 + 3 = 5

    double dist_1_3 = g.getShortestDistance(1, 3);
    EXPECT_DOUBLE_EQ(dist_1_3, 7.0); // 3 + 4 = 7


    // Also check for the precomputed shortest path
    auto precomputed_path = g.getPrecomputedShortestPath(0, 3);
    EXPECT_EQ(precomputed_path.size(), 4);

    auto precomputed_path_1_2 = g.getPrecomputedShortestPath(1, 2);
    EXPECT_EQ(precomputed_path_1_2.size(), 2);
}


TEST(TinyGraph, Line4RaeckeTree) {
    Graph g(4);
    g.addEdge(0,1,5);
    g.addEdge(1,2,0.5);
    g.addEdge(2,3,1);

   RaeckeFRTSolver solver;
    solver.solve(g);
    solver.storeFlow();
    solver.printFlow_();

    EXPECT_EQ(true, true);
}
