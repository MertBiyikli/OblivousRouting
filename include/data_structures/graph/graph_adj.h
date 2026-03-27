//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_GRAPH_ADJ_H
#define OBLIVIOUSROUTING_GRAPH_ADJ_H

#include "Igraph.h"

class GraphADJList : public IGraph {
public:
    std::vector<std::vector<int>> adjList;
    std::vector<std::vector<int>> edge_ids;
    std::vector<int> edge_id_start;
    std::vector<std::vector<double>> capacity;
    std::vector<std::vector<double>> distance;

    GraphADJList() = default;
    explicit GraphADJList(int n) : IGraph(n) {
        adjList.resize(n);
        capacity.resize(n);
        distance.reserve(n);
    }

    void finalize() override;

    int getNumDirectedEdges() const override;

    void addEdge(int u, int v, double cap, double distance = 1.0) override;

    NeighborRange neighbors(int node) const override;

    // Edge based indexing
    const int getEdgeId(int u, int v) const override;
    std::pair<int,int> getEdgeEndpoints(int e) const override;
    const int getAntiEdge(int e) const override;

    double getEdgeCapacity(int u, int v) const override;
    double getEdgeDistance(int u, int v) const override;
    bool updateEdgeDistance(int u, int v, double distance) override;
    bool updateEdgeCapacity(int u, int v, double capacity);

    double getEdgeCapacity(int e) const override;
    double getEdgeDistance(int e) const override;
    bool updateEdgeDistance(int e, double dist) override;

    void initializeMemberByParser(int maxNodeIdSeen) override;

    // Shortest path computations
    const double getDiameter() const override;
    std::vector<int> getShortestPath(int u, int v) const override;
    std::vector<double> getDistances(int u) const;
    double getShortestDistance(int u, int v) const;
    std::vector<int> getShortestPath(int s, int t, const std::vector<double>& dist_e) const override;





    void readGraph(const std::string &filename);
    virtual void print() const override;


    virtual void printGraphType() const override;

};

#endif //OBLIVIOUSROUTING_GRAPH_ADJ_H