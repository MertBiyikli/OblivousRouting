//
// Created by Mert Biyikli on 08.06.26.
//

#ifndef OBLIVIOUSROUTING_GRAPH_TEST_HELPERS_H
#define OBLIVIOUSROUTING_GRAPH_TEST_HELPERS_H

#pragma once

#include <algorithm>
#include <vector>

#include "data_structures/graph/Igraph.h"

// Shared graph builders for all IGraph implementations.
// Important: every builder returns a finalized graph, because most graph
// operations assume finalize() has been called.

template <typename Graph>
Graph makeEmptyGraph(int number_of_nodes) {
    Graph graph(number_of_nodes);
    graph.finalize();
    return graph;
}

template <typename Graph>
Graph makePathGraph4() {
    Graph graph(4);
    graph.addEdge(0, 1, 10.0, 1.0);
    graph.addEdge(1, 2, 20.0, 1.0);
    graph.addEdge(2, 3, 30.0, 1.0);
    graph.finalize();
    return graph;
}

template <typename Graph>
Graph makeWeightedPathGraph4() {
    Graph graph(4);
    graph.addEdge(0, 1, 10.0, 2.0);
    graph.addEdge(1, 2, 20.0, 3.0);
    graph.addEdge(2, 3, 30.0, 4.0);
    graph.finalize();
    return graph;
}

template <typename Graph>
Graph makeWeightedTriangleGraph() {
    Graph graph(3);
    graph.addEdge(0, 1, 10.0, 1.0);
    graph.addEdge(1, 2, 20.0, 2.0);
    graph.addEdge(0, 2, 30.0, 5.0);
    graph.finalize();
    return graph;
}

template <typename Graph>
Graph makeDiamondGraph() {
    Graph graph(5);
    graph.addEdge(0, 1, 10.0, 1.0);
    graph.addEdge(0, 2, 20.0, 1.0);
    graph.addEdge(1, 4, 30.0, 1.0);
    graph.addEdge(2, 4, 40.0, 1.0);
    graph.addEdge(0, 4, 50.0, 5.0);
    graph.finalize();
    return graph;
}

template <typename Range>
bool containsNeighbor(const Range& neighbors, int node) {
    return std::find(neighbors.begin(), neighbors.end(), node) != neighbors.end();
}

inline bool isValidPathFromTo(const std::vector<int>& path, int source, int target) {
    return !path.empty() && path.front() == source && path.back() == target;
}

inline bool containsNeighbor(const IGraph::NeighborRange& neighbors, int node) {
    for (int neighbor : neighbors) {
        if (neighbor == node) {
            return true;
        }
    }
    return false;
}
#endif //OBLIVIOUSROUTING_GRAPH_TEST_HELPERS_H