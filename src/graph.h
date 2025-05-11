//
// Created by Mert Biyikli on 11.05.25.
//

#ifndef OBLIVOUSROUTING_GRAPH_H
#define OBLIVOUSROUTING_GRAPH_H

#include <vector>
#include <stdexcept>

struct Edge{
    int target;
    double capacity;
};

class Graph{
public:
    Graph() : m_iNumNodes(0), m_iNumEdges(0) {}
    void addNode();
    void addEdge(int source, int target, double capacity);

    const std::vector<Edge>& neighbors(int node) const {
        return m_adj[node];
    }

    int numNodes() const;
    int numEdges() const;

private:
    std::vector<std::vector<Edge> > m_adj;

    int m_iNumNodes;
    int m_iNumEdges;
};


#endif //OBLIVOUSROUTING_GRAPH_H
