//
// Created by Mert Biyikli on 11.05.25.
//

#ifndef OBLIVOUSROUTING_GRAPH_H
#define OBLIVOUSROUTING_GRAPH_H

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <memory>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <queue>
#include "IGraph.h"
#include "../utils/hash.h"



class GraphADJ : public IGraph {
    // int m = 0; // Number of edges
    // std::vector<int> m_vertices; // List of vertices
    std::vector<std::vector<int>> m_adj;
    std::vector<std::vector<double>> m_adj_capacities, m_adj_distances;
    std::vector<std::vector<double>> m_distanceMatrix;
    std::vector<std::vector<int>> m_next; // next[i][j] = first node after i on the shortest path to j


    // this for also being able to use the edge based accessor
    std::vector<std::pair<int,int>> m_edgeList;                 // e -> (u, v)
    std::unordered_map<std::pair<int, int>, int, PairHash> m_uvToEdge;

public:
    GraphADJ() = default;

    explicit GraphADJ(int num_nodes)
        : IGraph(num_nodes),
          m_adj(num_nodes),
          m_adj_capacities(num_nodes),
          m_adj_distances(num_nodes)
    {
        // distanceMatrix / m_next can be built lazily
    }

    GraphADJ(const GraphADJ&)            = default;
    GraphADJ(GraphADJ&&) noexcept        = default;
    GraphADJ& operator=(const GraphADJ&) = default;
    GraphADJ& operator=(GraphADJ&&) noexcept = default;
    ~GraphADJ() override              = default;

    int getNumDirectedEdges() const override {
        return 2*m;
    }

    // derived methods
    IGraph::NeighborRange neighbors(int u) const override;
    void addEdge(int u , int v, double capacity, double distance = 1.0) override;

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
    const double& getShortestDistance(int u, int v) const;
    std::vector<int> getPrecomputedShortestPath(int u, int v) const;
    bool IsDistanceMatrixComputed() const;
    void createDistanceMatrix();




    void print() const;
    void readGraph(const std::string &filename);


    // TODO: need to implement this for the adjacency list based distances
    std::vector<int>
        getShortestPath(int s, int t, const std::vector<double>& dist_e) const override {
        return {};
    };

    void finalize() override {
        // do nothing
        return;
    }


    // void readLFGFile(const std::string& filename, bool withDistances);
};



#endif //OBLIVOUSROUTING_GRAPH_H
