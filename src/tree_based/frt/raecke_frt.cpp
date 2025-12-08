
#include "raecke_frt.h"
#include "mcct/mcct_derandomized_weighted_solver.h"

#include <queue>
#include <unordered_set>
#include <limits>
#include <cmath>
#include <iostream>

void RaeckeFRT::init(GraphADJ &g) {
    this->m_mcct.setGraph(g);
}


FRT_Tree RaeckeFRT::getTree(GraphADJ &g) {
    m_mcct.setGraph(g);
    setRequirements(g);
    computeNewDistances(g);
    return m_mcct.getBestTree();
}


void RaeckeFRT::computeRLoads(int treeIndex, FRT_Tree &_t, GraphADJ &copyGraph) {
    std::queue<std::shared_ptr<FRT_Node>> bfsQueue;
    for(auto& node : _t.GetRoot()->getChildren()) {
        bfsQueue.push(node);
    }

    if(m_idTree2edge2rload.size() <= treeIndex) {
        m_idTree2edge2rload.resize(treeIndex + 1);
    }
    m_idTree2edge2rload[treeIndex] = std::unordered_map<std::pair<int, int>, double, PairHash>();
    auto& edge2Load = m_idTree2edge2rload[treeIndex];

    while(!bfsQueue.empty()) {
        auto node = bfsQueue.front();
        bfsQueue.pop();
        auto parent = node->getParent();
        if(parent->getCenter() == node->getCenter()) {
            for(auto& child : node->getChildren()) {
                bfsQueue.push(child);
            }
            continue;
        }

        std::vector<int> nodeVertices;
        for (const int& v: node->GetVertices()) {
            nodeVertices.push_back(v);
        }
        std::vector<int> remaining(copyGraph.getVertices());
        for(const int& v : nodeVertices) {
            auto it = std::find(remaining.begin(), remaining.end(), v);
            if(it != remaining.end()) {
                remaining.erase(it);
            }
        }

        double cut = 0.0;
        for(int u : nodeVertices) {
            for(const auto& v : copyGraph.neighbors(u)) {

                if(std::find(remaining.begin(), remaining.end(), v) != remaining.end()) {
                    cut += copyGraph.getEdgeCapacity(u, v);
                }
            }
        }

        int centerParent = parent->getCenter();
        int currentCenter = node->getCenter();

        auto path = copyGraph.getShortestPath(centerParent, currentCenter);
        for(int i = 0; i<path.size()-1; i++) {
            int u = path[i];
            int v = path[i+1];
            // Normalize the edge distance based on the cut
            double rLoad = edge2Load[{u, v}];
            rLoad += cut / copyGraph.getEdgeCapacity(u, v); // Normalize by capacity
            edge2Load[{u, v}] = rLoad;
            edge2Load[{v,u}] = rLoad; // Also update the reverse arc
        }
        for(auto& child : node->getChildren()) {
            bfsQueue.push(child);
        }
    }
}




void RaeckeFRT::setRequirements(const GraphADJ &g) {
    for (const auto& u : g.getVertices()) {
        for (const auto& v : g.neighbors(u)) {
            double cap = g.getEdgeCapacity(u, v);
            m_mcct.addDemand(u, v, cap);
        }
    }
}
