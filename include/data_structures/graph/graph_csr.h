//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_GRAPH_CSR_H
#define OBLIVIOUSROUTING_GRAPH_CSR_H

#include "Igraph.h"

class GraphCSR : public IGraph {
public:

    /**
     * The graph members from, to and head represents the graph
     * edges e=(u,v) are stored from[e] = u and to[e] = v,
     * and head[u] gives the starting index of edges from u in the from and to arrays.
     */
    std::vector<int> to;
    std::vector<int> from;
    std::vector<int> head;

    /**
     * The temporary storage members tmp_edges, tmp_capacity and tmp_distance are used to store the edges,
     * capacities and distances before the graph is finalized.
     * Once finalize() is called, these temporary members are used
     * to populate the graph members and then cleared to free up memory.
     */
    std::vector<std::pair<int, int>> tmp_edges;
    std::vector<double> tmp_capacity;
    std::vector<double> tmp_distance;
    mutable std::vector<double> dist_buf;
    mutable std::vector<int> parent_buf;

    /**
     * After finalize() is called, the capacity and distance vectors are populated with the corresponding values for each edge.
     */
    std::vector<double> capacity;
    std::vector<double> distance;


    GraphCSR() = default;

    explicit GraphCSR(int num_nodes)
        : IGraph(num_nodes),
          head(num_nodes + 1, 0)
    {}

    GraphCSR(const GraphCSR&)            = default;
    GraphCSR(GraphCSR&&) noexcept        = default;
    GraphCSR& operator=(const GraphCSR&) = default;
    GraphCSR& operator=(GraphCSR&&) noexcept = default;
    ~GraphCSR() override                  = default;

    void finalize() override;
    void initializeMemberByParser(int maxNodeIdSeen) override;
    IGraph::NeighborRange neighbors(int u) const override;
    void addEdge(int u, int v, double capacity, double distance = 1.0) override;
    const int getEdgeId(int u, int v) const override;

    /**
     *  Edge based helper methods
     */
    std::pair<int, int> getEdgeEndpoints(int e) const override;
    double getEdgeCapacity(int e) const override;
    double getEdgeDistance(int e) const override;
    bool updateEdgeDistance(int e, double dist) override;
    const int getAntiEdge(int e) const override;

    /**
     * Node based helper methods
    */
    double getEdgeDistance(int u, int v) const override;
    double getEdgeCapacity(int u, int v) const override;
    bool updateEdgeDistance(int u, int v, double dist) override;

    /**
     * shortest path computation methods
     */
    std::vector<int> getShortestPath(int source, int target) const override;

    std::vector<int> getPathEdgesFromParent(int src, int tgt) const;

    std::vector<int> getShortestPath(int s, int t, const std::vector<double> &dist) const override;

    double getShortestDistance(int src, int tgt) const;

    const double getDiameter() const override;
    std::vector<int> getPathEdges(int src, int tgt) const;

    /**
     * Debugging methods
     */
    virtual void print() const override;
    virtual void printGraphType() const override;
};

#endif //OBLIVIOUSROUTING_GRAPH_CSR_H