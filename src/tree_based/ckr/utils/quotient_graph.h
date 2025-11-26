//
// Created by Mert Biyikli on 18.11.25.
//

#ifndef OBLIVIOUSROUTING_QUOTIENT_GRAPH_H
#define OBLIVIOUSROUTING_QUOTIENT_GRAPH_H

#include "../../../datastructures/graph.h"
#include "ultrametric_tree.h"
#include <vector>

namespace MendelScaling {

    template <typename G>
    struct QuotientLevel {
        G Gq;                                 // quotient graph at Δ
        std::vector<int> sigma_compact_of_v;      // size n: original vertex v -> compact qid in [0..k-1]
        std::vector<std::vector<int>> members_of_q; // size k: list of original vertices in each quotient node
    };

    // Build the Δ-level quotient graph and all mappings needed to map back to original vertices.
    // Complexity: O(m log n) amortized over all levels (by the paper’s edge-coverage argument).
    inline QuotientLevel<Graph> build_quotient_graph_with_map(const Graph& G, const UltrametricTree& ultra, double Delta) {
        const int n = G.getNumNodes();

        // 1) Map each original vertex to ancestor σΔ(v)
        std::vector<int> sigma_node(n);
        for (int v = 0; v < n; ++v) {
            sigma_node[v] = ultra.sigmaDelta(v, Delta); // ultra node id
        }

        // 2) Compact those ancestor IDs to 0..k-1
        // std::unordered_map<int,int> idmap;
        // idmap.reserve(n*2);
        int next_id = 0;
        std::vector<int> qid_of_sigma(ultra.N, -1);
        for (int v = 0; v < n; ++v) {
            int s = sigma_node[v];
            int &ref = qid_of_sigma[s];
            if (ref == -1) {
                ref = next_id++;
                // idmap[ref] = v;
            }
        }

        // 3) Original v -> compact quotient id
        std::vector<int> sigma_compact_of_v(n);
        sigma_compact_of_v.reserve(n);
        for (int v = 0; v < n; ++v) {
            sigma_compact_of_v[v] = qid_of_sigma[sigma_node[v]];
        }

        // k is the number of quotient nodes
        const int k = next_id;

        // 4) Members of each quotient node
        std::vector<std::vector<int>> members_of_q(k);
        for (int v = 0; v < n; ++v) {
            members_of_q[sigma_compact_of_v[v]].push_back(v);
        }

        // 5) Build quotient edges with min inter-edge weight
        Graph Gq(k);
        std::unordered_map<long long,std::pair<double, double>> min_w;
        min_w.reserve((size_t)G.getNumEdges());
        auto key = [](int a,int b){ return ((long long)a<<32) | (unsigned)b; };

        for (int u = 0; u < n; ++u) {
            int cu = sigma_compact_of_v[u];
            for (int v : G.neighbors(u)) {
                if (u >= v) continue; // undirected, process each edge once
                int cv = sigma_compact_of_v[v];
                if (cu == cv) continue;
                int a = (cu < cv) ? cu : cv;
                int b = (cu < cv) ? cv : cu;
                double cap = G.getEdgeCapacity(u, v);
                double w = G.getEdgeDistance(u, v);
                long long K = key(a,b);
                auto it = min_w.find(K);
                if (it == min_w.end() || w < it->second.second)
                    min_w[K] = {cap, w};
            }
        }

        // TODO: here we need to ensure that we add the edge distances correctly and not the capacities
        for (auto& [K, cap_and_weight] : min_w) {
            int a = (int)(K >> 32);
            int b = (int)(K & 0xffffffff);
            Gq.addEdge(a,b, cap_and_weight.first, cap_and_weight.second);
        }

        return {.Gq = std::move(Gq),
                .sigma_compact_of_v = std::move(sigma_compact_of_v),
                .members_of_q = std::move(members_of_q)};
    }


    // Build the Δ-level quotient graph and all mappings needed to map back to original vertices.
    // Complexity: O(m log n) amortized over all levels (by the paper’s edge-coverage argument).
    inline QuotientLevel<Graph_csr> build_quotient_graph_with_map(const Graph_csr& G, const UltrametricTree& ultra, double Delta) {
        const int n = G.getNumNodes();

        // 1) Map each original vertex to ancestor σΔ(v)
        std::vector<int> sigma_node(n);
        for (int v = 0; v < n; ++v) {
            sigma_node[v] = ultra.sigmaDelta(v, Delta); // ultra node id
        }

        // 2) Compact those ancestor IDs to 0..k-1
        // std::unordered_map<int,int> idmap;
        // idmap.reserve(n*2);
        int next_id = 0;
        std::vector<int> qid_of_sigma(ultra.N, -1);
        for (int v = 0; v < n; ++v) {
            int s = sigma_node[v];
            int &ref = qid_of_sigma[s];
            if (ref == -1) {
                ref = next_id++;
                // idmap[ref] = v;
            }
        }

        // 3) Original v -> compact quotient id
        std::vector<int> sigma_compact_of_v(n);
        sigma_compact_of_v.reserve(n);
        for (int v = 0; v < n; ++v) {
            sigma_compact_of_v[v] = qid_of_sigma[sigma_node[v]];
        }

        // k is the number of quotient nodes
        const int k = next_id;

        // 4) Members of each quotient node
        std::vector<std::vector<int>> members_of_q(k);
        for (int v = 0; v < n; ++v) {
            members_of_q[sigma_compact_of_v[v]].push_back(v);
        }

        // 5) Build quotient edges with min inter-edge weight
        Graph_csr Gq(k);
        std::unordered_map<long long,std::pair<double, double>> min_w;
        min_w.reserve((size_t)G.getNumEdges());
        auto key = [](int a,int b){ return ((long long)a<<32) | (unsigned)b; };

        for (int u = 0; u < n; ++u) {
            int cu = sigma_compact_of_v[u];
            for (auto&  v : G.neighbors(u)) {
                if (u >= v) continue; // undirected, process each edge once
                int cv = sigma_compact_of_v[v];
                if (cu == cv) continue;
                int a = (cu < cv) ? cu : cv;
                int b = (cu < cv) ? cv : cu;
                double cap = G.getEdgeCapacity(u, v);
                double w = G.getEdgeDistance(u, v);
                long long K = key(a,b);
                auto it = min_w.find(K);
                if (it == min_w.end() || w < it->second.second)
                    min_w[K] = {cap, w};
            }
        }

        // TODO: here we need to ensure that we add the edge distances correctly and not the capacities
        for (auto& [K, cap_and_weight] : min_w) {
            int a = (int)(K >> 32);
            int b = (int)(K & 0xffffffff);
            Gq.addEdge(a,b, cap_and_weight.first, cap_and_weight.second);
        }

        Gq.finalize();

        return {.Gq = std::move(Gq),
                .sigma_compact_of_v = std::move(sigma_compact_of_v),
                .members_of_q = std::move(members_of_q)};
    }
}

#endif //OBLIVIOUSROUTING_QUOTIENT_GRAPH_H