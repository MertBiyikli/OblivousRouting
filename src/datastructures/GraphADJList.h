//
// Created by Mert Biyikli on 25.02.26.
//

#ifndef OBLIVIOUSROUTING_GRAPHADJLIST_H
#define OBLIVIOUSROUTING_GRAPHADJLIST_H

#include <iostream>

#include "IGraph.h"

class GraphADJList : public IGraph {
public:
    std::vector<std::vector<int>> adjList;
    std::vector<std::vector<int>> edge_ids;
    // edge_id_start[u] = the first (smallest) edge id assigned to node u; last element is total edges (size = n+1)
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

    void addEdge(int u, int v, double cap, double distance) override;

    NeighborRange neighbors(int node) const override {
        return NeighborRange{adjList[node].data(), adjList[node].data() + adjList[node].size()};
    }

    // Edge based indexing
    int getEdgeId(int u, int v) const override;
    std::pair<int,int> edgeEndpoints(int e) const override;
    const int getAntiEdge(int e) const override;

    double getEdgeCapacity(int u, int v) const override;
    double getEdgeDistance(int u, int v) const override;
    bool updateEdgeDistance(int u, int v, double distance) override;
    bool updateEdgeCapacity(int u, int v, double capacity);

    double getEdgeCapacity(int e) const override;
    double getEdgeDistance(int e) const override;
    bool updateEdgeDistance(int e, double dist) override;

    void InitializeMemberByParser(int maxNodeIdSeen) override;

    // Shortest path computations
    const double GetDiameter() const override;
    std::vector<int> getShortestPath(int u, int v) const override;
    std::vector<double> getDistances(int u) const;
    double getShortestDistance(int u, int v) const;





    void readGraph(const std::string &filename);
    virtual void print() const override {
        for (size_t i = 0; i < adjList.size(); ++i) {
            std::cout << "Node " << i << ": ";
            for (size_t j = 0; j < adjList[i].size(); ++j) {
                std::cout << "(" << adjList[i][j] << ", " << capacity[i][j] << ") ";
            }
            std::cout << std::endl;
        }
    }


    virtual void printGraphType() const override {
        std::cout << "GraphADJList (Adjacency List) format" << std::endl;
    }

    std::vector<int>
        getShortestPath(int s, int t, const std::vector<double>& dist_e) const override;
};


#endif //OBLIVIOUSROUTING_GRAPHADJLIST_H