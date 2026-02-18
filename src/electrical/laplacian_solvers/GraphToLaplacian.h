//
// Created by Mert Biyikli on 05.02.26.
//

#ifndef OBLIVIOUSROUTING_GRAPHTOLAPLACIAN_H
#define OBLIVIOUSROUTING_GRAPHTOLAPLACIAN_H
#include <vector>
#include "../../datastructures/IGraph.h"


/*
 * Struct for efficiently bookkeeping the matrix entries synched with the graph edges
 */
struct GraphToLaplacian {
    std::vector<int> to, from, head; // to[e], from[e] gives the endpoints of edge e
    std::vector<double> weights; // weights[e] gives the weight of edge e
    std::vector<int> laplacian_indices; // laplacian_indices[e] gives the index in the Laplacian matrix values array
    std::vector<int> laplacian_indices_for_diagonal_elements; // laplacian_indices_red[e] gives the index in the reduced Laplacian matrix values array

    void init(IGraph& graph) {
        int n = graph.getNumNodes();
        int m = graph.getNumEdges();

        to.resize(m);
        from.resize(m);
        head.resize(n + 1, 0);
        weights.resize(m);
        laplacian_indices.resize(m);
        laplacian_indices_for_diagonal_elements.resize(n);

        for (int e = 0; e < m; e++) {
            from[e] = graph.edgeEndpoints(e).first;
            to[e]   = graph.edgeEndpoints(e).second;
        }

        int curr = 0;
        for (int u = 0; u < n; u++) {
            while (curr < m && from[curr] == u) curr++;
            head[u+1] = curr;
        }
    }

    double getEdgeWeight(int e) const {
        assert(e >= 0 && e < weights.size());
        return weights[e];
    }

    double getEdgeWeight(int u, int v) const {
        assert(u >= 0 && u < head.size() - 1);
        assert(v >= 0 && v < head.size() - 1);
        for (int idx = head[u]; idx < head[u + 1]; ++idx) {
            if (to[idx] == v) {
                return weights[idx];
            }
        }
        throw std::invalid_argument("Edge not found in getEdgeWeight(u,v)");
    }

    void setEdgeWeight(int e, double w) {
        assert(e >= 0 && e < weights.size());
        weights[e] = w;
    }

    void setEdgeWeight(int u, int v, double w) {
        assert(u >= 0 && u < head.size() - 1);
        assert(v >= 0 && v < head.size() - 1);
        for (int idx = head[u]; idx < head[u + 1]; ++idx) {
            if (to[idx] == v) {
                weights[idx] = w;
                return;
            }
        }
        throw std::invalid_argument("Edge not found in setEdgeWeight(u,v)");
    }

    void setLaplacianIndex(int u, int v, int idx) {
        assert(u >= 0 && u < head.size() - 1);
        assert(v >= 0 && v < head.size() - 1);
        if (u == v) {
            laplacian_indices_for_diagonal_elements[u] = idx;
            return;
        }
        for (int e = head[u]; e < head[u + 1]; ++e) {
            if (to[e] == v) {
                laplacian_indices[e] = idx;
                return;
            }
        }
        throw std::invalid_argument("Edge not found in setLaplacianIndex(u,v)");
    }

    int getLaplacianIndex(int u, int v) const {
        assert(u >= 0 && u < head.size() - 1);
        assert(v >= 0 && v < head.size() - 1);
        if (u == v) {
            return laplacian_indices_for_diagonal_elements[u];
        }
        for (int e = head[u]; e < head[u + 1]; ++e) {
            if (to[e] == v) {
                return laplacian_indices[e];
            }
        }
        throw std::invalid_argument("Edge not found in getLaplacianIndex(u,v)");
    }
};

#endif //OBLIVIOUSROUTING_GRAPHTOLAPLACIAN_H