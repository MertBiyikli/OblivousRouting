//
// Created by Mert Biyikli on 20.03.26.
//

#include "../../include/data_structures/mendel_scaling/quotient_graph.h"
#include "../../include/data_structures/graph/graph_csr.h"
#include <algorithm>
#include <unordered_map>

void QuotientConstruction::preprocessEdges(const IGraph& G) {
     assert(G.getNumNodes() > 0);
     original_n = G.getNumNodes();
     weights.clear();
     weights.reserve(G.getNumUndirectedEdges());
     for (int u = 0; u < original_n; ++u) {
         for (auto& v : G.neighbors(u)) {
             if (u < v) {
                 double w = G.getEdgeDistance(u, v);
                 int e_id = G.getEdgeId(u, v);
                 weights.emplace_back(e_id, w);
             }
         }
     }


     std::sort(weights.begin(), weights.end(),
               [](auto& a, auto& b){ return a.second < b.second; });
 }


QuotientLevel QuotientConstruction::constructQuotientGraph(const UltrametricTree& ultra, double Delta, IGraph& G) {
    assert(original_n > 0);

    // Map each original vertex to ancestor σΔ(v)
    std::vector<int> sigma_node(original_n);
    for (int v = 0; v < original_n; ++v) {
        sigma_node[v] = ultra.sigmaDelta(v, Delta); // ultra node id
    }

    // Compact those ancestor IDs to 0..k-1
    int next_id = 0;
    std::vector<int> qid_of_sigma(ultra.N, -1);
    for (int v = 0; v < original_n; ++v) {
        int s = sigma_node[v];
        int &ref = qid_of_sigma[s];
        if (ref == -1) {
            ref = next_id++;
        }
    }

    // Original v -> compact quotient id
    std::vector<int> sigma_compact_of_v(original_n);
    sigma_compact_of_v.reserve(original_n);
    for (int v = 0; v < original_n; ++v) {
        sigma_compact_of_v[v] = qid_of_sigma[sigma_node[v]];
    }

    // k is the number of quotient nodes
    const int k = next_id;

    // Members of each quotient node
    std::vector<std::vector<int>> members_of_q(k);
    for (int v = 0; v < original_n; ++v) {
        members_of_q[sigma_compact_of_v[v]].push_back(v);
    }


    // 5) Build quotient edges with SLIDING WINDOW method:

    // find the interval of edges that are within the sliding window of size Delta
    // interval: [(left, right) = edges with weight in [Delta, Delta /2n]
    std::unique_ptr<IGraph> Gq_ = std::make_unique<GraphCSR>(k);

    double w_low = Delta / (2.0 * original_n);
    double w_high = Delta/2.0;
    size_t left = 0, right = 0;

    while (right < weights.size() && weights[right].second <= w_high) {
        ++right;
    }
    while (left < weights.size() && weights[left].second < w_low) {
        ++left;
    }

    // two passes to cover all edges
    std::unordered_map<long long,std::pair<double, double>> min_w;
    min_w.reserve((size_t)(right - left + 1));
    auto key = [](int a,int b){ return ((long long)a<<32) | (unsigned)b; };

    for (size_t i = left; i < right; ++i) {
        //auto& [w, u, v] = edges[i];
        const auto& [e_id, w2] = weights[i];
        const auto& [v1, v2] = G.getEdgeEndpoints(e_id);
        int cu = sigma_compact_of_v[v1];
        int cv = sigma_compact_of_v[v2];
        if (cu == cv) continue;
        int a = (cu < cv) ? cu : cv;
        int b = (cu < cv) ? cv : cu;
        double cap = G.getEdgeCapacity(v1, v2);
        long long K = key(a,b);
        auto it = min_w.find(K);
        if (it == min_w.end() || w2 < it->second.second)
            min_w[K] = {cap, w2};
    }

    for (auto& [K, cap_and_weight] : min_w) {
        int a = (int)(K >> 32);
        int b = (int)(K & 0xffffffff);
        Gq_->addEdge(a,b, cap_and_weight.first, cap_and_weight.second);
    }

    Gq_->finalize();
    return {.Gq = std::move(Gq_),
            .sigma_compact_of_v = std::move(sigma_compact_of_v),
            .members_of_q = std::move(members_of_q)};
}