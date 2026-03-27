//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_IGRAPH_H
#define OBLIVIOUSROUTING_IGRAPH_H

#include <vector>
#include <numeric>

#define INVALID_EDGE_ID -1

class IGraph {
public:
    IGraph() = default;
    explicit IGraph(int n) : n(n), m(0) {
        vertices.resize(n);
        std::iota(vertices.begin(), vertices.end(), 0);
        is_processed = false;
    }
    IGraph(const IGraph&) = default;
    IGraph(IGraph&&) noexcept = default;
    IGraph& operator=(const IGraph&) = default;
    IGraph& operator=(IGraph&&) noexcept = default;

    virtual ~IGraph() = default;

    /**
     * shared members by all graph child classes
     */
    int n, m;
    bool is_processed;
    std::vector<int> vertices;


    /**
     * Once the graph is read, finalize should be invoked to
     * settle graph instance for proper usage.
     * Can leave graph in a bad state if not called.
     */
    virtual void finalize() = 0;

    /**
     * This function is meant to be called by the parser before reading the graph,
     * to initialize the members of the graph instance and setting the number of nodes that will be seen.
     */
    virtual void initializeMemberByParser(int maxNodeIdSeen) = 0;

    /**
     * Adds an edge between nodes u and v with the given capacity and distance.
     * Assume that u < n and v < n. So member n should be set before adding edges
     */
    virtual void addEdge(int u, int v, double capacity, double distance = 1) = 0;

    /**
     * Returns the number of nodes and edges in the graph, as well as the list of node IDs.
     */
    virtual int getNumNodes() const {
        return n;
    }

    virtual int getNumUndirectedEdges() const {
        return m/2;
    }

    virtual int getNumDirectedEdges() const {
        return m;
    }

    virtual const std::vector<int>& getVertices() const {
        return vertices;
    }

    /**
     * NeighborRange provides a convenient way to access the neighbors of a node without exposing the underlying data structure.
     */
    struct NeighborRangeBase {
        virtual ~NeighborRangeBase() = default;
        virtual const int* beginPtr() const = 0;
        virtual const int* endPtr() const = 0;
    };

    struct NeighborRange {
        const int* b;
        const int* e;

        const int* begin() const { return b; }
        const int* end()   const { return e; }

        const int size() const { return static_cast<int>(e - b); }
        int operator[](size_t i) const {
            return b[i];
        }
    };
    virtual NeighborRange neighbors(int node) const = 0;

    /**
     * Node based access of edge information
     */
    virtual double getEdgeCapacity(int u, int v) const = 0;
    virtual double getEdgeDistance(int u, int v) const = 0;
    virtual bool updateEdgeDistance(int u, int v, double distance) = 0;
    virtual const int getEdgeId(int u, int v) const = 0;

    /**
     * Edge based access of edge information
     */
    virtual double getEdgeCapacity(int e) const = 0;
    virtual double getEdgeDistance(int e) const = 0;
    virtual bool updateEdgeDistance(int e, double dist) = 0;
    virtual std::pair<int, int> getEdgeEndpoints(int e) const = 0;
    virtual const int getAntiEdge(int e) const = 0;


    /**
     * Set the edge distance back to uniform
     *
     */
    virtual void resetEdgeDistance() {
        for (int v = 0; v < n; ++v) {
            for (size_t i = 0; i < neighbors(v).size(); ++i) {
                updateEdgeDistance(v, neighbors(v)[i], 1.0);
            }
        }
    }

    /**
     * Returns the list of node IDs in the shortest path from source to target, including both endpoints.
     * If no path exists, returns an empty vector. The underlying shortest path algorithm is by default Dijkstra algorithm.
     */
    virtual std::vector<int> getShortestPath(int source, int target) const = 0;
    virtual const double getDiameter() const = 0;
    virtual std::vector<int> getShortestPath(int s, int t, const std::vector<double>& dist) const = 0;
    double getShortestPathDistance(int s, int t) const {
        auto path = getShortestPath(s, t);
        if (path.empty()) {
            return std::numeric_limits<double>::infinity(); // No path exists
        }
        double total_dist = 0.0;
        for (size_t i = 0; i + 1 < path.size(); ++i) {
            total_dist += getEdgeDistance(path[i], path[i + 1]);
        }
        return total_dist;
    }

    /**
     * Debugging method
     */
    virtual void print() const = 0;
    virtual void printGraphType() const = 0;
};

#endif //OBLIVIOUSROUTING_IGRAPH_H