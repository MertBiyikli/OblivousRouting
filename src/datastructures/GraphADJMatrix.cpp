//
// Created by Mert Biyikli on 25.02.26.
//

#include "GraphADJMatrix.h"
#include <cassert>
#include <numeric>
#include <stdexcept>
#include <queue>
#include <algorithm>
#include <limits>

// ---------------------------------------------------------------------------
// addEdge — primary matrix write
// ---------------------------------------------------------------------------
void GraphADJMatrix::addEdge(int u, int v, double cap, double dist) {
    assert(u >= 0 && u < n);
    assert(v >= 0 && v < n);
    capacity_matrix[u][v] = cap;
    distance_matrix[u][v] = dist;

    capacity_matrix[v][u] = cap;
    distance_matrix[v][u] = dist;
}

// ---------------------------------------------------------------------------
// finalize — scan matrix row-by-row, assign edge ids, build edge_id_start
// ---------------------------------------------------------------------------
void GraphADJMatrix::finalize() {
    edge_id_matrix.assign(n, std::vector<int>(n, INVALID_EDGE_ID));
    neighbors_cache.assign(n, {});
    edge_id_start.resize(n + 1);

    int edge_id = 0;
    for (int u = 0; u < n; ++u) {
        edge_id_start[u] = edge_id;
        for (int v = 0; v < n; ++v) {
            if (capacity_matrix[u][v] != 0.0) {
                edge_id_matrix[u][v] = edge_id++;
                neighbors_cache[u].push_back(v);  // v is visited in order, so list is sorted
            }
        }
    }
    edge_id_start[n] = edge_id;
    m = edge_id;
    is_processed = true;
}

// ---------------------------------------------------------------------------
// getNumDirectedEdges
// ---------------------------------------------------------------------------
int GraphADJMatrix::getNumDirectedEdges() const {
    return m;
}

// ---------------------------------------------------------------------------
// neighbors — NeighborRange backed by neighbors_cache built in finalize()
// ---------------------------------------------------------------------------
IGraph::NeighborRange GraphADJMatrix::neighbors(int node) const {
    // neighbors_cache[node] is built in finalize()
    return NeighborRange{neighbors_cache[node].data(),
                         neighbors_cache[node].data() + neighbors_cache[node].size()};
}

// ---------------------------------------------------------------------------
// getEdgeId — O(1) matrix lookup
// ---------------------------------------------------------------------------
int GraphADJMatrix::getEdgeId(int u, int v) const {
    if (u < 0 || u >= n || v < 0 || v >= n) return INVALID_EDGE_ID;
    return edge_id_matrix[u][v];
}

// ---------------------------------------------------------------------------
// edgeEndpoints — O(log n) via binary search on edge_id_start,
//                  then O(n) scan across row u to find column with that id.
// ---------------------------------------------------------------------------
std::pair<int,int> GraphADJMatrix::edgeEndpoints(int e) const {
    if (!is_processed) throw std::runtime_error("edgeEndpoints: graph not finalized");
    if (e < 0 || e >= m)  throw std::out_of_range("edgeEndpoints: edge id out of range");

    // find owning row u
    auto it = std::upper_bound(edge_id_start.begin(), edge_id_start.end(), e);
    int u = static_cast<int>(std::distance(edge_id_start.begin(), it)) - 1;
    if (u < 0 || u >= n) throw std::runtime_error("edgeEndpoints: cannot locate row");

    // within row u, find column v whose edge_id_matrix entry equals e.
    // edge ids are assigned left-to-right across the row, so offset = e - edge_id_start[u]
    // but we need the offset-th neighbour in the row, not the v-th column.
    int offset = e - edge_id_start[u];
    int v = neighbors_cache[u][offset];
    return {u, v};
}

// ---------------------------------------------------------------------------
// getAntiEdge
// ---------------------------------------------------------------------------
const int GraphADJMatrix::getAntiEdge(int e) const {
    auto [u, v] = edgeEndpoints(e);
    int anti = getEdgeId(v, u);
    if (anti == INVALID_EDGE_ID) throw std::runtime_error("getAntiEdge: reverse edge not found");
    return anti;
}

// ---------------------------------------------------------------------------
// capacity / distance — node-pair versions (direct matrix access)
// ---------------------------------------------------------------------------
double GraphADJMatrix::getEdgeCapacity(int u, int v) const {
    if (u < 0 || u >= n || v < 0 || v >= n) throw std::out_of_range("getEdgeCapacity: out of range");
    return capacity_matrix[u][v];
}

double GraphADJMatrix::getEdgeDistance(int u, int v) const {
    if (u < 0 || u >= n || v < 0 || v >= n) throw std::out_of_range("getEdgeDistance: out of range");
    return distance_matrix[u][v];
}

bool GraphADJMatrix::updateEdgeDistance(int u, int v, double distance) {
    if (u < 0 || u >= n || v < 0 || v >= n) return false;
    if (capacity_matrix[u][v] == 0.0) return false;
    distance_matrix[u][v] = distance;
    return true;
}


// ---------------------------------------------------------------------------
// edge-id versions
// ---------------------------------------------------------------------------
double GraphADJMatrix::getEdgeCapacity(int e) const {
    auto [u, v] = edgeEndpoints(e);
    return capacity_matrix[u][v];
}

double GraphADJMatrix::getEdgeDistance(int e) const {
    auto [u, v] = edgeEndpoints(e);
    return distance_matrix[u][v];
}

bool GraphADJMatrix::updateEdgeDistance(int e, double dist) {
    auto [u, v] = edgeEndpoints(e);
    return updateEdgeDistance(u, v, dist);
}

// ---------------------------------------------------------------------------
// InitializeMemberByParser
// ---------------------------------------------------------------------------
void GraphADJMatrix::InitializeMemberByParser(int maxNodeIdSeen) {
    IGraph::n = maxNodeIdSeen + 1;
    IGraph::m = 0;
    capacity_matrix.assign(n, std::vector<double>(n, 0.0));
    distance_matrix.assign(n, std::vector<double>(n, std::numeric_limits<double>::infinity()));
    edge_id_matrix.assign(n, std::vector<int>(n, INVALID_EDGE_ID));
    neighbors_cache.assign(n, {});
    edge_id_start.clear();
    IGraph::vertices.resize(n);
    std::iota(IGraph::vertices.begin(), IGraph::vertices.end(), 0);
}

// ---------------------------------------------------------------------------
// Shortest path — Dijkstra over the distance matrix
// ---------------------------------------------------------------------------
static std::vector<int> dijkstra(
        const std::vector<std::vector<double>>& dist_mat, int s, int t) {
    int n = static_cast<int>(dist_mat.size());
    std::vector<double> dist(n, std::numeric_limits<double>::infinity());
    std::vector<int>    parent(n, -1);
    using PQ = std::pair<double,int>;
    std::priority_queue<PQ, std::vector<PQ>, std::greater<PQ>> pq;
    dist[s] = 0.0;
    pq.emplace(0.0, s);
    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (d > dist[u]) continue;
        if (u == t) break;
        for (int v = 0; v < n; ++v) {
            if (dist_mat[u][v] == 0.0) continue;   // no edge
            double nd = d + dist_mat[u][v];
            if (nd < dist[v]) {
                dist[v] = nd;
                parent[v] = u;
                pq.emplace(nd, v);
            }
        }
    }
    if (dist[t] == std::numeric_limits<double>::infinity()) return {};
    std::vector<int> path;
    for (int cur = t; cur != -1; cur = parent[cur]) path.push_back(cur);
    std::reverse(path.begin(), path.end());
    return path;
}

std::vector<int> GraphADJMatrix::getShortestPath(int u, int v) const {
    if (!is_processed) throw std::runtime_error("getShortestPath: not finalized");
    return dijkstra(distance_matrix, u, v);
}

std::vector<int> GraphADJMatrix::getShortestPath(
        int s, int t, const std::vector<double>& dist_e) const {
    if (!is_processed) throw std::runtime_error("getShortestPath: not finalized");
    // Build a temporary distance matrix from per-edge weights
    std::vector<std::vector<double>> tmp(n, std::vector<double>(n, 0.0));
    for (int u = 0; u < n; ++u) {
        for (int v = 0; v < n; ++v) {
            if (edge_id_matrix[u][v] != INVALID_EDGE_ID) {
                tmp[u][v] = dist_e[edge_id_matrix[u][v]];
            }
        }
    }
    return dijkstra(tmp, s, t);
}



const double GraphADJMatrix::GetDiameter() const {
    double diameter = 0.0;
    for (int u = 0; u < n; ++u) {
        for (int v = u + 1; v < n; ++v) {
            auto path = dijkstra(distance_matrix, u, v);
            if (path.empty()) continue;
            double d = 0.0;
            for (size_t i = 0; i + 1 < path.size(); ++i)
                d += distance_matrix[path[i]][path[i+1]];
            diameter = std::max(diameter, d);
        }
    }
    return diameter;
}

// ---------------------------------------------------------------------------
// print
// ---------------------------------------------------------------------------
void GraphADJMatrix::print() const {
    std::cout << "Adjacency matrix (capacity):\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            std::cout << capacity_matrix[i][j] << ' ';
        std::cout << '\n';
    }
}
