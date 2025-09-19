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

void RaeckeMST::init(Graph& _g) {
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

double RaeckeMST::iterate(int treeIndex) {
    MSTTree t = getTree(m_graph);
    m_trees.push_back(t);


    computeRLoads(treeIndex, t, m_graph);
    double l = getMaxRload(treeIndex);
    double delta = std::min(1/l, 1-m_lambdaSum);
    m_lambdas.push_back(delta);

    m_graphs.push_back(Graph(m_graph));

    return delta;
}


MSTTree RaeckeMST::getTree(Graph &g) {
    std::vector<std::pair<int,int>> mst_edges = mst.build_mst(g);
    MSTTree t = mst.build_tree(g, mst_edges, 0);

    // Todo: maybe think of a more lets say clean approach of recomputing the distance
    computeNewDistances(g);
    return t;
}

void RaeckeMST::computeRLoads(int treeIndex, MSTTree &_t, Graph& g) {

    if(m_idTree2edge2rload.size() <= treeIndex) {
        m_idTree2edge2rload.resize(treeIndex + 1);
    }



    m_idTree2edge2rload[treeIndex] = std::unordered_map<std::pair<int, int>, double>();
    auto& edge2Load = m_idTree2edge2rload[treeIndex];

    // collect all tree edges that are incidient to one of the endpoints from e
    std::unordered_set<std::pair<int, int>> treeEdges;

    for (int u = 0; u < g.getNumNodes(); u++) {
        for (const auto& v : g.neighbors(u)) {

            treeEdges.clear();
            // compute the rloads for each edge e in the original graph

            // iterate through the parent and children in _t to find the tree edges
            for (auto& child : _t.children[u]) {
                std::pair<int, int> edge = {std::min(u, child), std::max(u, child)};
                treeEdges.insert(edge);
            }

            for (auto& child : _t.children[v]) {
                std::pair<int, int> edge = {std::min(v, child), std::max(v, child)};
                treeEdges.insert(edge);
            }

            // parent edge
            if ( _t.parent[v] != v) {
                std::pair<int, int> edge = {std::min(v, _t.parent[v]), std::max(v, _t.parent[v])};
                treeEdges.insert(edge);
            }
            if ( _t.parent[u] != u) {
                std::pair<int, int> edge = {std::min(u, _t.parent[u]), std::max(u, _t.parent[u])};
                treeEdges.insert(edge);
            }

            double rload = 0;
            // compute the rload for each edge
            for (const auto& edge : treeEdges) {
                rload += g.getEdgeCapacity(edge.first, edge.second);
            }

            rload /= g.getEdgeCapacity(u, v);
            // add to the edge2Load map
            m_idTree2edge2rload[treeIndex][{u, v}] = rload;
            m_idTree2edge2rload[treeIndex][{v, u}] = rload; // undirected graph

        }
    }

}


