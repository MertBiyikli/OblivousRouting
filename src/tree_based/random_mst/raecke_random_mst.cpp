//
// Created by Mert Biyikli on 18.09.25.
//

#include "raecke_random_mst.h"
#include "../frt/raecke_frt.h"

#include <queue>
#include <unordered_set>
#include <limits>
#include <cmath>
#include <iostream>

void RaeckeMST::init(GraphADJ& _g) {
    this->mst.setGraph(_g);
    // Sample MST
    uint64_t seed = dist(rng);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> U(0.0,1.0);
    for(int u = 0; u < _g.getNumNodes(); u++) {
        for (const auto& v : _g.neighbors(u)) {
            double randVal = U(rng);
            _g.updateEdgeDistance(u, v, randVal);
        }
    }
}
/*
double RaeckeMST::iterate(int treeIndex) {
    MSTTree t = getTree(m_graph);
    m_trees.push_back(t);


    computeRLoads(treeIndex, t, m_graph);
    double l = getMaxRload(treeIndex);
    double delta = std::min(1/l, 1-m_lambdaSum);
    m_lambdas.push_back(delta);

    m_graphs.push_back(Graph(m_graph));

    return delta;
}*/


MSTTree RaeckeMST::getTree(GraphADJ &g) {
    // Todo: maybe think of a more lets say clean approach of recomputing the distance
    computeNewDistances(g);
    mst.setGraph(g);
    std::vector<std::pair<int,int>> mst_edges = mst.build_mst();
    MSTTree t = mst.build_tree(g, mst_edges, 0);
    return t;
}

// raecke_random_mst.cpp
void RaeckeMST::computeRLoads(int treeIndex, MSTTree& _t, GraphADJ& g) {
    if (m_idTree2edge2rload.size() <= treeIndex) {
        m_idTree2edge2rload.resize(treeIndex + 1);
    }
    auto& edge2Load = m_idTree2edge2rload[treeIndex];
    edge2Load.clear();

    constexpr double EPS = 1e-12;

    // For each original (undirected) edge, compute rload via the unique MST path
    for (int u = 0; u < g.getNumNodes(); ++u) {
        for (int v : g.neighbors(u)) {
            if (u > v) continue; // handle each undirected edge once

            // Get MST path (sequence of vertices) between u and v
            std::vector<int> path = RandomMST::getMSTPath(u, v, _t.parent);

            // Sum capacities along the MST path
            double path_cap_sum = 0.0;
            for (size_t i = 1; i < path.size(); ++i) {
                int a = path[i-1], b = path[i];
                path_cap_sum += g.getEdgeCapacity(a, b);
            }

            // Denominator: capacity of (u,v), floored to avoid div-by-zero
            double cap_uv = g.getEdgeCapacity(u, v);
            if (cap_uv < EPS) cap_uv = EPS;

            double rload = path_cap_sum / cap_uv;

            // Store both directions for convenience
            edge2Load[{u, v}] = rload;
            edge2Load[{v, u}] = rload;
        }
    }
}



