//
// Created by Mert Biyikli on 11.05.25.
//

#include <gtest/gtest.h>
#include "../../src/graph.h"

TEST(GraphTest, AddAndReorderEdges) {
    Graph g;
    g.addNode();
    g.addNode();
    g.addNode();

    g.addEdge(0, 1, 10.0);
    g.addEdge(0, 2, 5.0);
    g.addEdge(1, 2, 3.0);
    g.addEdge(2, 3, 1.0);



    for (int u = 0; u < g.numNodes(); ++u) {
        for (const auto& e : g.neighbors(u)) {
            EXPECT_GE(e.capacity, 0.0);
        }
    }
}
