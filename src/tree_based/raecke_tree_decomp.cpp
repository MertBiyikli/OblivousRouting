
#include "raecke_tree_decomp.h"
#include "mcct/mcct_derandomized_weighted_solver.h"

#include <queue>
#include <unordered_set>
#include <limits>
#include <cmath>
#include <iostream>


Graph RaeckeFRT::getGraph() const {
    return m_graph;
}

void RaeckeFRT::setGraph(const Graph &g) {
    m_graph = g;
}



std::vector<FRT_Tree> RaeckeFRT::getTrees() const {
    return m_trees;
}

std::vector<double> RaeckeFRT::getLambdas() const {
    return m_lambdas;
}

std::vector<Graph> RaeckeFRT::getGraphs() const {
    return m_graphs;
}


void RaeckeFRT::run() {
    auto p = std::make_shared<Graph>(m_graph);
    this->m_mcct.setGraph(p);

    m_lambdaSum = 0.0;
    int treeIndex = 0;

    while(m_lambdaSum < 1) {
        if(debug)
            std::cout << "Computing tree " << treeIndex << std::endl;
        m_lambdaSum += iterate(treeIndex);
        treeIndex++;
    }
}

double RaeckeFRT::iterate(int treeIndex) {
    auto copyGraph = std::make_shared<Graph>(m_graph);

    FRT_Tree t = getTree(copyGraph);
    m_trees.push_back(t);
    m_mcct.reset();

    computeRLoads(treeIndex, t, copyGraph);
    double l = getMaxRload(treeIndex, t);
    double delta = std::min(1/l, 1-m_lambdaSum);
    m_lambdas.push_back(delta);
    // m_lambdas.push_back(delta);
    m_graphs.push_back(*copyGraph);
    return delta;
}

FRT_Tree RaeckeFRT::getTree(std::shared_ptr<Graph> &g) {
    m_mcct.setGraph(g);
    setRequirements(g);
    computeNewDistances(g);
    return m_mcct.getBestTree(debug);
}


void RaeckeFRT::computeRLoads(int treeIndex, FRT_Tree &_t, std::shared_ptr<Graph> &copyGraph) {
    std::queue<std::shared_ptr<FRT_Node>> bfsQueue;
    for(auto& node : _t.GetRoot()->getChildren()) {
        bfsQueue.push(node);
    }

    if(m_idTree2edge2rload.size() <= treeIndex) {
        m_idTree2edge2rload.resize(treeIndex + 1);
    }
    m_idTree2edge2rload[treeIndex] = std::unordered_map<std::pair<int, int>, double>();
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
        std::vector<int> remaining(copyGraph->getVertices());
        for(const int& v : nodeVertices) {
            auto it = std::find(remaining.begin(), remaining.end(), v);
            if(it != remaining.end()) {
                remaining.erase(it);
            }
        }

        double cut = 0.0;
        for(int u : nodeVertices) {
            for(const auto& v : copyGraph->neighbors(u)) {

                if(std::find(remaining.begin(), remaining.end(), v) != remaining.end()) {
                    cut += copyGraph->getEdgeCapacity(u, v);
                }
            }
        }

        int centerParent = parent->getCenter();
        int currentCenter = node->getCenter();

        auto path = copyGraph->getShortestPath(centerParent, currentCenter);
        for(int i = 0; i<path.size()-1; i++) {
            int u = path[i];
            int v = path[i+1];
            // Normalize the edge distance based on the cut
            double rLoad = edge2Load[{u, v}];
            rLoad += cut / copyGraph->getEdgeCapacity(u, v); // Normalize by capacity
            edge2Load[{u, v}] = rLoad;
            edge2Load[{v,u}] = rLoad; // Also update the reverse arc
        }
        for(auto& child : node->getChildren()) {
            bfsQueue.push(child);
        }
    }
}


double RaeckeFRT::getMaxRload(int treeIndex, FRT_Tree &_t) {
    double maxRatio = 0.0;
    for(const auto& [edge, rLoad] : m_idTree2edge2rload[treeIndex]) {
        if (rLoad > maxRatio) {
            maxRatio = rLoad;
        }
    }
    return maxRatio;
}



void RaeckeFRT::setRequirements(const std::shared_ptr<Graph> &g) {
    for (const auto& u : g->getVertices()) {
        for (const auto& v : g->neighbors(u)) {
            double cap = g->getEdgeCapacity(u, v);
            m_mcct.addDemand(u, v, cap);
        }
    }
}

void RaeckeFRT::computeNewDistances(std::shared_ptr<Graph> &g) {
    double totalRLoadsAllEdges = this->getRloadAllEdges(g);
    std::unordered_map<std::pair<int, int>, double> edge2scaledDist;

    for(int u = 0; u<g->getNumNodes(); u++) {
        for (const auto& v : g->neighbors(u)) {
            double totalrLoad = 0.0;
            for (size_t i = 0; i < m_trees.size(); ++i) {
                auto it = this->m_idTree2edge2rload[i].find({u, v});
                double rLoad = (it != this->m_idTree2edge2rload[i].end() ? it->second : 0.0);
                totalrLoad += rLoad * m_lambdas[i];
            }
            double num         = std::exp(totalrLoad) / g->getEdgeCapacity(u, v);
            double newDistance = num / totalRLoadsAllEdges;

            if(std::isinf(newDistance) || std::isinf(num)) {
                throw std::runtime_error("Infinity encountered in newDistance or Exponential calculation for edge: " + std::to_string(u) + " → " + std::to_string(v));
            }

            if(std::isnan(newDistance)) {
                throw std::runtime_error("NaN encountered in newDistance calculation for edge: " + std::to_string(u) + " → " + std::to_string(v));
            }

            edge2scaledDist[{u, v}] = newDistance;
        }
    }

    // Normalize distances
    normalizeDistance(g, edge2scaledDist);
}


void RaeckeFRT::normalizeDistance(std::shared_ptr<Graph> &_g, std::unordered_map<std::pair<int, int>, double> &edge2scaledDist) {
    double minDistance = std::numeric_limits<double>::infinity();
    for(int u = 0; u < _g->getNumNodes(); u++) {
        for (const auto& v : _g->neighbors(u)) {
            double scaledDist = edge2scaledDist[{u, v}];
            if (scaledDist < minDistance) {
                minDistance = scaledDist;
            }
        }
    }

    for(int u = 0; u < _g->getNumNodes(); u++) {
        for (const auto &v: _g->neighbors(u)) {
            double arc2scaledDistValue = edge2scaledDist[{u, v}];
            double newDistance = arc2scaledDistValue / minDistance;

            if (std::isnan(newDistance)) {
                throw std::runtime_error(
                        "NaN encountered in newDistance for edge: " + std::to_string(u) + " → " + std::to_string(v));
            }
            if (newDistance < 1) {
                newDistance = 1.0; // Ensure minimum distance is 1}
            }
            _g->updateEdgeDistance(u, v, newDistance);
        }
    }
}


double RaeckeFRT::getRloadAllEdges(const std::shared_ptr<Graph>& g) const{
    double totalRLoadAllEdges = 0.0;
    for(int u = 0; u < g->getNumNodes(); u++) {
        for (const auto& v : g->neighbors(u)) {
            double totalPerEdge = 0.0;

            std::pair<int, int> edge = {u, v};
            // Sum over all trees i:  rLoad(i,arc) * lambda[i]
            for (size_t i = 0; i < m_trees.size(); ++i) {
                auto it = this->m_idTree2edge2rload[i].find(edge);
                if(it != this->m_idTree2edge2rload[i].end()) {
                    totalPerEdge += it->second * m_lambdas[i];
                } else {
                    // If the edge is not found, we assume the r-load is 0
                    totalPerEdge += 0.0;
                }
            }

            // Add e^(totalPerEdge) to the grand total
            totalRLoadAllEdges += std::exp(totalPerEdge);
        }
    }
    return totalRLoadAllEdges;
}
