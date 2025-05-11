//
// Created by Mert Biyikli on 11.05.25.
//


#include "graph.h"
#include <stdexcept>
#include <iostream>

int Graph::numNodes() const {
    return m_iNumNodes;
}

int Graph::numEdges() const {
    return m_iNumEdges;
}

void Graph::addEdge(int source, int target, double capacity) {
    if (source < 0 || source >= m_iNumNodes || target < 0 || target >= m_iNumNodes) {
        throw std::out_of_range("Invalid node index");
    }
    m_adj[source].push_back({target, capacity});
    m_adj[target].push_back({source, capacity});
    m_iNumEdges++;
}


void Graph::print() const {
    for (int i = 0; i < m_iNumNodes; ++i) {
        std::cout << "Node " << i << ": ";
        for (const auto& edge : m_adj[i]) {
            std::cout << "(" << edge.target << ", " << edge.capacity << ") ";
        }
        std::cout << std::endl;
    }
}