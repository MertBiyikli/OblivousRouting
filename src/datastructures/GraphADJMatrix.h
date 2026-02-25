//
// Created by Mert Biyikli on 25.02.26.
//

#ifndef OBLIVIOUSROUTING_GRAPHADJMATRIX_H
#define OBLIVIOUSROUTING_GRAPHADJMATRIX_H

#include "IGraph.h"
#include <vector>
#include <iostream>

class GraphADJMatrix : public IGraph {
public:
    // Primary storage — true adjacency matrix
    // capacity_matrix[u][v] > 0  ⟺  directed edge u→v exists
    std::vector<std::vector<double>> capacity_matrix;
    std::vector<std::vector<double>> distance_matrix;

    // edge_id_matrix[u][v] = edge id of u→v, INVALID_EDGE_ID if no edge.
    // Assigned in finalize(); O(1) lookup for getEdgeId(u,v).
    std::vector<std::vector<int>> edge_id_matrix;

    // edge_id_start[u] = first edge id assigned to row u  (size = n+1).
    // Used by edgeEndpoints(e): binary-search to find owning row u → O(log n).
    std::vector<int> edge_id_start;

    // neighbors_cache[u] = sorted list of v where capacity_matrix[u][v] != 0.
    // Built in finalize(). Needed so NeighborRange can return a contiguous int*.
    std::vector<std::vector<int>> neighbors_cache;

    GraphADJMatrix() = default;
    explicit GraphADJMatrix(int n) : IGraph(n) {
        capacity_matrix.assign(n, std::vector<double>(n, 0.0));
        distance_matrix.assign(n, std::vector<double>(n, std::numeric_limits<double>::infinity()));
        edge_id_matrix.assign(n, std::vector<int>(n, INVALID_EDGE_ID));
    }

    void finalize() override;
    int getNumDirectedEdges() const override;
    void addEdge(int u, int v, double cap, double dist) override;

    // Iterate over row u: raw pointer to row, length n
    NeighborRange neighbors(int node) const override;

    // Edge-based indexing
    int getEdgeId(int u, int v) const override;
    std::pair<int,int> edgeEndpoints(int e) const override;
    const int getAntiEdge(int e) const override;

    // Node-pair accessors
    double getEdgeCapacity(int u, int v) const override;
    double getEdgeDistance(int u, int v) const override;
    bool updateEdgeDistance(int u, int v, double distance) override;


    // Edge-id accessors (delegate to node-pair versions via edgeEndpoints)
    double getEdgeCapacity(int e) const override;
    double getEdgeDistance(int e) const override;
    bool updateEdgeDistance(int e, double dist) override;

    void InitializeMemberByParser(int maxNodeIdSeen) override;

    const double GetDiameter() const override;
    std::vector<int> getShortestPath(int u, int v) const override;
    std::vector<int> getShortestPath(int s, int t, const std::vector<double>& dist_e) const override;

    void print() const override;
    void printGraphType() const override {
        std::cout << "GraphADJMatrix (Adjacency Matrix) format" << std::endl;
    }
};

#endif //OBLIVIOUSROUTING_GRAPHADJMATRIX_H