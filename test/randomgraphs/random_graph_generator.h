//
// Created by Mert Biyikli on 11.05.25.
//

#ifndef OBLIVOUSROUTING_RANDOM_GRAPH_GENERATOR_H
#define OBLIVOUSROUTING_RANDOM_GRAPH_GENERATOR_H

#include "../../src/graph.h"
#include <random>
#include <set>
#include <utility>

class RandomGraphGenerator {
public:
    static Graph generate(int numNodes, int numEdges, double minCap = 1.0, double maxCap = 10.0) {
        if (numEdges > numNodes * (numNodes - 1) / 2) {
            throw std::invalid_argument("Too many edges for undirected graph");
        }

        Graph g;
        for (int i = 0; i < numNodes; ++i) {
            g.addNode();
        }

        std::set<std::pair<int, int>> existing;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distNode(0, numNodes - 1);
        std::uniform_real_distribution<> distCap(minCap, maxCap);

        while (g.numEdges() < numEdges) {
            int u = distNode(gen);
            int v = distNode(gen);
            if (u == v) continue;

            auto edge = std::make_pair(std::min(u, v), std::max(u, v));

            if (existing.count(edge)) continue;

            double cap = distCap(gen);
            g.addEdge(u, v, cap);
            g.addEdge(v, u, cap); // undirected

            existing.insert(edge);
        }

        return g;
    }
};

#endif //OBLIVOUSROUTING_RANDOM_GRAPH_GENERATOR_H
