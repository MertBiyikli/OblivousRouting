//
// Created by Mert Biyikli on 25.11.25.
//

#ifndef OBLIVIOUSROUTING_GRAPH_BASE_H
#define OBLIVIOUSROUTING_GRAPH_BASE_H

#pragma once
#include "graph_types.h"
#include <gtest/gtest.h>


template <typename GraphType>
class IGraphTest : public ::testing::Test {
public:
    GraphType G;

    void buildSimpleTriangle() {
        // Graph with 3 nodes: 0-1-2-0
        G = GraphType(3);
        G.addEdge(0, 1, 1.0);
        G.addEdge(1, 2, 1.0);
        G.addEdge(2, 0, 1.0);

        // CSR must finalize
        if constexpr (std::is_same_v<GraphType, Graph_csr>) {
            G.finalize();
        }
    }
};

TYPED_TEST_SUITE_P(IGraphTest);



#endif //OBLIVIOUSROUTING_GRAPH_BASE_H