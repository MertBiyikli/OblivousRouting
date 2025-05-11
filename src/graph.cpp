//
// Created by Mert Biyikli on 11.05.25.
//


#include "graph.h"
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>


void Graph::InitNodes(int nodes) {
    m_iNumNodes = nodes;
    m_adj.resize(nodes);
}

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

void Graph::readGraph(const std::string &filename) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    int n = 0, m = 0;

    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == 'c') continue;

        std::istringstream iss(line);
        if (line[0] == 'p') {
            std::string tmp;
            iss >> tmp >> tmp >> n >> m;
            break;
        }
    }

    if (n == 0) {
        throw std::runtime_error("Invalid or missing problem line in DIMACS file");
    }

    this->InitNodes(n);

    do {
        if (line.empty() || line[0] == 'c') continue;

        std::istringstream iss(line);
        if (line[0] == 'a') {
            char dummy;
            int u, v;
            double cap = 1.0;
            iss >> dummy >> u >> v;
            if (!(iss >> cap)) cap = 1.0;
            this->addEdge(u - 1, v - 1, cap);  // convert to 0-based
        }
    } while (std::getline(infile, line));


}
