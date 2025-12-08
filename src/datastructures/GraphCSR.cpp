#include <vector>

#include "GraphCSR.h"

std::vector<int>
GraphCSR::getShortestPath(int src, int tgt, const std::vector<double>& dist_e) const {
    using P = std::pair<double, int>;
    const double INF = std::numeric_limits<double>::infinity();

    if (src < 0 || src >= n || tgt < 0 || tgt >= n) {
        throw std::out_of_range("GraphCSR::getShortestPath(IDistance): node index out of range");
    }

    // Ensure reusable buffers have correct size
    if (static_cast<int>(dist_buf.size()) != n) {
        dist_buf.resize(n);
        parent_buf.resize(n);
    }

    std::fill(dist_buf.begin(),   dist_buf.end(),   INF);
    std::fill(parent_buf.begin(), parent_buf.end(), -1);

    auto& d     = dist_buf;
    auto& prev  = parent_buf;

    d[src] = 0.0;
    std::priority_queue<P, std::vector<P>, std::greater<>> pq;
    pq.emplace(0.0, src);

    while (!pq.empty()) {
        auto [du, u] = pq.top();
        pq.pop();

        if (du > d[u]) continue;
        if (u == tgt) break;

        // iterate over CSR neighbors of u
        for (int e = head[u]; e < head[u + 1]; ++e) {
            int v = to[e];

            double w = dist_e[e];
            if (w < 0.0) continue; // optional: forbid negative weights

            double nd = du + w;
            if (nd < d[v]) {
                d[v]    = nd;
                prev[v] = u;
                pq.emplace(nd, v);
            }
        }
    }

    std::vector<int> path;
    if (d[tgt] == INF) {
        return path; // no path
    }

    for (int v = tgt; v != -1; v = prev[v]) {
        path.push_back(v);
    }
    std::reverse(path.begin(), path.end());
    return path;
}
