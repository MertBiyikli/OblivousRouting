//
// Created by Mert Biyikli on 05.11.25.
//

#ifndef OBLIVIOUSROUTING_GRAPH_CSR_H
#define OBLIVIOUSROUTING_GRAPH_CSR_H
#include <algorithm>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>
#include <iostream>
#include <sstream>
#include <queue>
#include "IGraph.h"
#include <map>
#include "../utils/priority_queue.h"



/*
 * The idea is to have a graph class that uses CSR format internally
 * to store edges, capacities and distances for faster access.
 */
class GraphCSR : public IGraph {
public:

    // std::vector<std::pair<int, int>> edges;
    std::vector<int> to;         // size m, neighbor v for each edge
    std::vector<int> from;       // size m, source u (optional)
    std::vector<int> head;
    std::vector<double> capacity;
    std::vector<double> distance;
    bool is_processed = false;

    // Temporary storage before finalize()
    std::vector<std::pair<int, int>> tmp_edges;
    std::vector<double> tmp_capacity;
    std::vector<double> tmp_distance;

    mutable std::vector<double> dist_buf;
    mutable std::vector<int> parent_buf;


    GraphCSR() = default;

    explicit GraphCSR(int num_nodes)
        : IGraph(num_nodes),
          from(),
          to(),
          head(num_nodes + 1, 0),
          capacity(),
          distance()
    {}

    GraphCSR(const GraphCSR&)            = default;
    GraphCSR(GraphCSR&&) noexcept        = default;
    GraphCSR& operator=(const GraphCSR&) = default;
    GraphCSR& operator=(GraphCSR&&) noexcept = default;
    ~GraphCSR() override                  = default;


    /*
     * The finalize method is used to preprocess the edges into index sorted vectors.
     * This method should be called after all edges have been added using addEdge().
     * Updating edge distances doesnt require finalize, only adding edges does.
     * If unsure whether finalize has been called, it is safe to call it again as it will
     * do nothing if the graph is already processed.
     */
    void finalize() override {
    if (is_processed) {
        return;
    }
        // 1) sort edges by (u,v)
        std::vector<int> idx(m);
        std::iota(idx.begin(), idx.end(), 0);

        if (tmp_edges.empty()) {
            // fill in the edges if none were added
            for (int e = 0; e < m; e++) {
                int u = from[e];
                int v = to[e];
                tmp_edges.emplace_back(u, v);
                tmp_capacity.push_back(1.0);
                tmp_distance.push_back(1.0);
            }
        }

        std::stable_sort(idx.begin(), idx.end(),
            [&](int a, int b){
                if (tmp_edges[a].first != tmp_edges[b].first)
                    return tmp_edges[a].first < tmp_edges[b].first;
                return tmp_edges[a].second < tmp_edges[b].second;
            });

        // 2) allocate CSR storage
        to.resize(m);
        from.resize(m);
        capacity.resize(m);
        distance.resize(m);
        head.assign(n+1, 0);

        // 3) fill CSR arrays in sorted order
        for (int i = 0; i < m; i++) {
            const int k = idx[i];
            from[i] = tmp_edges[k].first;
            to[i]   = tmp_edges[k].second;
            capacity[i] = tmp_capacity[k];
            distance[i] = tmp_distance[k];
        }

        // 4) build head[]
        int curr = 0;
        for (int u = 0; u < n; u++) {
            while (curr < m && from[curr] == u) curr++;
            head[u+1] = curr;
        }

        // 5) clear temporary storage
        tmp_edges.clear();
        tmp_capacity.clear();
        tmp_distance.clear();

        // 6) check for consistency
        for (int u = 0; u < n; u++) {
            int start = head[u];
            int end   = head[u+1];

            for (int e = start; e < end - 1; e++) {
                if (to[e] == to[e+1]) {
                    is_processed = false;
                    throw std::runtime_error("GraphCSR::finalize: duplicate edges detected");
                }
            }
        }

        is_processed = true;
    }

    int getNumDirectedEdges() const override {
        return m;
    }

    void addEdge(int u, int v, double cap, double dist = 1.0) override {

        if (u<0 || u>=n || v<0 || v>=n)
            throw std::out_of_range("GraphCSR::addEdge: node index");

        // only allow undirected edges to be passed
        if (u > v) {return;};

        // undirected edge => add (u,v) and (v,u)
        tmp_edges.emplace_back(u,v);
        tmp_capacity.push_back(cap);
        tmp_distance.push_back(dist);

        tmp_edges.emplace_back(v,u);
        tmp_capacity.push_back(cap);
        tmp_distance.push_back(dist);
        m+=2;


        is_processed = false;
    }

    const int getAntiEdge(int e) const override {
        const auto& [u, v] = edgeEndpoints(e);
        int anti_e = getEdgeId(v, u);
        if (anti_e == INVALID_EDGE_ID) {
            throw std::runtime_error("getAntiEdge: reverse edge not found");
        }
        return anti_e;
    }


    // returns edge id for (u,v), or -1 if not exists
    int getEdgeId(int u, int v) const override{
        if (!is_processed)
            throw std::runtime_error("GraphCSR::getEdgeId: graph not finalized");

        int start = head[u];
        int end   = head[u+1];

        for (int e = start; e < end; e++) {
            if (to[e] == v) return e;
        }
        return INVALID_EDGE_ID;
    }

    virtual std::pair<int, int> edgeEndpoints(int e) const override {
        if (!is_processed) {
            throw std::runtime_error("edgeEndpoints: graph not finalized");
        }
        if (e < 0 || e >= m) {
            throw std::out_of_range("edgeEndpoints: edge id out of range");
        }
        return {from[e], to[e]};
}


    double getEdgeCapacity(int edge_id) const override {
        if (edge_id >= m || edge_id < 0) {
            throw std::out_of_range("Edge index out of range");
        }
        if (capacity.size() < edge_id) {
            throw std::runtime_error("Edge index exceeds capacity size");
        }
        return capacity[edge_id];
    }

    double getEdgeCapacity(int u, int v) const override {
        if ( u > v) std::swap(u, v);

        int e = getEdgeId(u, v);
        if (e == INVALID_EDGE_ID)
            throw std::runtime_error("Edge not found");
        return capacity[e];
    }

    double getEdgeDistance(int edge_id) const override {
        if (edge_id >= m || edge_id < 0) {
            throw std::out_of_range("Edge index out of range");
        }
        if (distance.size() < edge_id) {
            throw std::runtime_error("Edge index exceeds distance size");
        }
        return distance[edge_id];
    }

    double getEdgeDistance(int u, int v) const override  {
        int source = std::min(u, v);
        int target = std::max(u, v);

        int head_idx = head[source];
        int end_idx = (source + 1 < n) ? head[source + 1] : m;
        for (int e = head_idx; e < end_idx; ++e) {
            const int& a = from[e];
            const int& b = to[e];
            if (a != source) break; // moved past source's edges
            if (b == target) {
                return distance[e];
            }
        }
        throw std::runtime_error("Edge not found");
    }

    bool updateEdgeDistance(int e, double dist) override {
        bool ok = false;
        if (e > m || e < 0) {
            throw std::out_of_range("Edge index out of range");
        }
        if (distance.size() < e) {
            throw std::out_of_range("Edge index exceeds distance size");
        }
        distance[e] = dist;

        // also update the reverse edge if undirected
        const int& a = from[e];
        const int& b = to[e];
        int rev_edge_id = getEdgeId(a, b);
        if (rev_edge_id != INVALID_EDGE_ID) {
            distance[rev_edge_id] = dist;
            ok = true;
        }else {
            throw std::runtime_error("Reverse edge not found");
        }
        return ok;
    }


    // TODO: this function can be optimized && also should handle the case when u > v
    bool updateEdgeDistance(int u, int v, double dist) override {
        int source = u;
        int target = v;

        if (source < 0 || source >= n || target < 0 || target >= n) {
            throw std::out_of_range("Node index out of range");
        }

        if (source > target) std::swap(source, target);

        bool ok = false;

        int e_source_target = getEdgeId(source, target);
        if (e_source_target != INVALID_EDGE_ID) {
            distance[e_source_target] = dist;
            ok = true;
        }else {
            throw std::runtime_error("Edge not found");
        }

        if ( ok ) {
            // also update the reverse edge if undirected
            int rev_edge_id = getEdgeId(target, source);
            if (rev_edge_id != INVALID_EDGE_ID) {
                distance[rev_edge_id] = dist;
                ok &= true;
            } else {
                ok = false;
                throw std::runtime_error("Reverse edge not found");
            }
        }
        return ok;
    }


    std::vector<int> getShortestPath(int src, int tgt = -1) const override {
        const double INF = std::numeric_limits<double>::infinity();

        if (src < 0 || src >= n)
            throw std::out_of_range("GraphCSR::getShortestPath: src or tgt out of range");

        // --- resize reusable buffers only if needed ---
        if ((int)dist_buf.size() != n) {
            dist_buf.resize(n);
            parent_buf.resize(n);
        }

        auto& dist   = dist_buf;
        auto& parent = parent_buf;

        std::fill(dist.begin(),   dist.end(),   INF);
        std::fill(parent.begin(), parent.end(), -1);

        // --- initialize distance of source ---
        dist[src] = 0.0;

        // --- fast heap with handles ---
        MinHeap<double,int> pq(n);
        pq.insert(src, 0.0);

        while (!pq.empty()) {
            int u     = pq.top();
            double du = pq.topKey();
            pq.deleteTop();

            // skip stale entries
            if (du > dist[u]) continue;

            // early exit: we can stop once tgt is extracted
            if (u == tgt) break;

            // iterate CSR edges of u
            int start = head[u];
            int end   = head[u+1];

            for (int e = start; e < end; ++e) {
                int v = to[e];
                double w = distance[e];      // customized weights (parameterized)
                double nd = du + w;

                if (nd < dist[v]) {
                    dist[v] = nd;
                    parent[v] = u;
                    pq.insertOrAdjustKey(v, nd);
                }
            }
        }
        // If we don't have a specific target, skip reconstruction
        if (tgt < 0)
            return {};
        // --- reconstruct path ---
        std::vector<int> path;
        if (dist[tgt] == INF) return path; // unreachable

        for (int v = tgt; v != -1; v = parent[v])
            path.push_back(v);

        std::reverse(path.begin(), path.end());
        return path;
    }


    std::vector<int> getShortestPath(int s, int t, const std::vector<double>& dist) const override;


    // ===============================================================
    // 🔹 Extract edge IDs along the path found by the last Dijkstra run
    // ===============================================================
    std::vector<int> getPathEdgesFromParent(int src, int tgt) const {
        std::vector<int> edges_on_path;
        if (tgt < 0 || tgt >= n) return edges_on_path;
        if (parent_buf.empty()) return edges_on_path;

        for (int v = tgt; v != src && v != -1; v = parent_buf[v]) {
            int u = parent_buf[v];
            if (u == -1) break;

            // Find edge ID for (u,v)
            int start = head[u];
            int end = head[u + 1];
            for (int e = start; e < end; ++e) {
                if (to[e] == v) {
                    edges_on_path.push_back(e);
                    break;
                }
            }
        }

        std::reverse(edges_on_path.begin(), edges_on_path.end());
        return edges_on_path;
    }

    // ===============================================================
    // 🔹 Compute shortest path and return edge indices along it
    // ===============================================================
    std::vector<int> getPathEdges(int src, int tgt) const {
        getShortestPath(src, tgt);               // fills dist_buf and parent_buf
        return getPathEdgesFromParent(src, tgt); // reconstruct via parent links
    }




    double getShortestDistance(int src, int tgt) const {
        const double INF = std::numeric_limits<double>::infinity();

        // run Dijkstra (fills dist_buf)
        getShortestPath(src, tgt);

        // ensure target index is valid
        if (tgt < 0 || tgt >= n) return INF;

        // return the precomputed distance value
        return dist_buf[tgt];
    }


    const double GetDiameter() const override{
        if (n == 0) return 0.0;

        double diameter = 0.0;

        for (int src = 0; src < n; ++src) {
            // We only need the distances, not the path to any particular target.
            // Passing tgt = src just triggers Dijkstra from src (the early-exit won't trigger).
            getShortestPath(src, -1);  // full Dijkstra from src

            // reuse the existing dist_buf from the last run
            const auto &dist = dist_buf;

            // find the farthest reachable node
            double local_max = 0.0;
            for (double d : dist)
                if (d < std::numeric_limits<double>::infinity())
                    local_max = std::max(local_max, d);

            diameter = std::max(diameter, local_max);
        }

        return diameter;
    }



    void print() {
        for (int e = 0; e<m; e++) {
            std::cout << "Edge " << e << " =("<< from[e] << ", "  << to[e] << " cap:" << capacity[e] << " dist: " << distance[e] << ")" << std::endl;
        }

        // print head array
        std::cout << "Head array: ";
        for (int i = 0; i < n; i++) {
            std::cout << head[i] << " ";
        }
    }
    IGraph::NeighborRange neighbors(int u) const override {
        if (!is_processed)
            throw std::runtime_error("GraphCSR::neighbors: graph not finalized");

        int start = head[u];
        int end   = head[u+1];
        return NeighborRange{ &to[start], &to[end] };
    }


    // Note that this function does change the members from the base class
    void InitializeMemberByParser(int maxNodeIdSeen) override{
        IGraph::n = maxNodeIdSeen + 1;
        IGraph::m = 0;
        IGraph::vertices.resize(this->n);
        std::iota(IGraph::vertices.begin(), IGraph::vertices.end(), 0);
    }

};

#endif //OBLIVIOUSROUTING_GRAPH_CSR_H