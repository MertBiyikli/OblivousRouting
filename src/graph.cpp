//
// Created by Mert Biyikli on 11.05.25.
//


#include "graph.h"

int Graph::numNodes() const {
    return m_iNumNodes;
}

int Graph::numEdges() const {
    return m_iNumEdges;
}

void Graph::addNode() {
    m_adj.push_back(std::vector<Edge>());
    m_iNumNodes++;
}

void Graph::addEdge(int source, int target, double capacity) {
    if (source < 0 || source >= m_iNumNodes || target < 0 || target >= m_iNumNodes) {
        throw std::out_of_range("Invalid node index");
    }
    m_adj[source].push_back({target, capacity});
    m_iNumEdges++;
}
