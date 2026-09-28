//
// Created by Mert Biyikli on 20.03.26.
//
#include "algorithms/oblivious/mwu/oracle/electrical/graph_to_laplacian.h"
#include <cassert>
#include <stdexcept>

void GraphToLaplacian::init(const IGraph& graph) {
    int n = graph.getNumNodes();
    int m = graph.getNumDirectedEdges();

    to.resize(m);
    from.resize(m);
    head.resize(n + 1, 0);
    weights.resize(m);
    laplacian_indices.resize(m);
    laplacian_indices_for_diagonal_elements.resize(n);

    for (int e = 0; e < m; e++) {
        from[e] = graph.getEdgeEndpoints(e).first;
        to[e]   = graph.getEdgeEndpoints(e).second;
    }

    int curr = 0;
    for (int u = 0; u < n; u++) {
        while (curr < m && from[curr] == u) curr++;
        head[u+1] = curr;
    }
}

void GraphToLaplacian::init(const optimized::Graph<EdgeData>& graph) {
    int n = graph.getNumNodes();
    int m = graph.getNumDirectedEdges();

    to.resize(m);
    from.resize(m);
    head.resize(n + 1, 0);
    weights.resize(m);
    laplacian_indices.resize(m);
    laplacian_indices_for_diagonal_elements.resize(n);

    int edge_idx = 0;
    for (int u = 0; u < n; u++) {
        for (const auto& e : graph.edgesOf(u)) {
            from[edge_idx] = u;
            to[edge_idx] = e.tail;
            edge_idx++;
        }
    }

    int curr = 0;
    for (int u = 0; u < n; u++) {
        while (curr < m && from[curr] == u) curr++;
        head[u+1] = curr;
    }
}


void GraphToLaplacian::init(const std::vector<std::pair<int, int>>& edges, const std::vector<double>& edge_weights, int n) {
    int m = static_cast<int>(edges.size());

    to.resize(m);
    from.resize(m);
    head.resize(n + 1, 0);
    weights = edge_weights;
    laplacian_indices.resize(m);
    laplacian_indices_for_diagonal_elements.resize(n);

    for (int e = 0; e < m; e++) {
        from[e] = edges[e].first;
        to[e]   = edges[e].second;
    }

    int curr = 0;
    for (int u = 0; u < n; u++) {
        while (curr < m && from[curr] == u) curr++;
        head[u+1] = curr;
    }
}

double GraphToLaplacian::getEdgeWeight(int e) const {
    assert(e >= 0 && e < weights.size());
    return weights[e];
}

double GraphToLaplacian::getEdgeWeight(int u, int v) const {
    assert(u >= 0 && u < head.size() - 1);
    assert(v >= 0 && v < head.size() - 1);
    for (int idx = head[u]; idx < head[u + 1]; ++idx) {
        if (to[idx] == v) {
            return weights[idx];
        }
    }
    assert(false && "Edge not found in getEdgeWeight(u,v)");
    return 0.0;
}

void GraphToLaplacian::setEdgeWeight(int e, double w) {
    assert(e >= 0 && e < weights.size());
    weights[e] = w;
}

Result<void> GraphToLaplacian::setEdgeWeight(int u, int v, double w) {
    assert(u >= 0 && u < head.size() - 1);
    assert(v >= 0 && v < head.size() - 1);
    for (int idx = head[u]; idx < head[u + 1]; ++idx) {
        if (to[idx] == v) {
            weights[idx] = w;
            return {};
        }
    }
    return makeErrorMessage(ErrorCode::InvalidArgument, "Edge not found in laplacian model.");
}

Result<void> GraphToLaplacian::setLaplacianIndex(int u, int v, int idx) {
    assert(u >= 0 && u < head.size() - 1);
    assert(v >= 0 && v < head.size() - 1);
    if (u == v) {
        laplacian_indices_for_diagonal_elements[u] = idx;
        return {};
    }
    for (int e = head[u]; e < head[u + 1]; ++e) {
        if (to[e] == v) {
            laplacian_indices[e] = idx;
            return {};
        }
    }
    return makeErrorMessage(ErrorCode::InvalidArgument, "Edge not found in laplacian model.");
}

Result<int> GraphToLaplacian::getLaplacianIndex(int u, int v) const {
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
    return makeErrorMessage(ErrorCode::InvalidArgument, "Edge not found in laplacian model.");
}