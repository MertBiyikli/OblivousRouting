//
// Created by Mert Biyikli on 18.11.25.
//
#include "ultrametric_tree.h"
#include <queue>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include "../../random_mst/mst.h"

using namespace MendelScaling;

void UltrametricTree::buildFromMST(
    int n_,
    const std::vector<std::pair<int,int>>& mst_edges,
    const std::vector<double>& mst_w
) {
    n = n_;
    const int m = (int)mst_edges.size();

    // 1) Init leaves
    T.clear();
    T.resize(n);
    for (int v = 0; v < n; ++v) {
        T[v].Gamma = 0.0;
        T[v].parent = -1;
        T[v].child.clear();
    }

    // 2) Standard DSU on [0..n-1]
    std::vector<int> dsu(n);
    std::iota(dsu.begin(), dsu.end(), 0);

    auto find = [&](int x) {
        while (dsu[x] != x) {
            dsu[x] = dsu[dsu[x]];
            x = dsu[x];
        }
        return x;
    };

    auto unite_roots = [&](int a, int b) {
        a = find(a); b = find(b);
        if (a == b) return a;
        dsu[a] = b;
        return b;
    };

    // 3) Map DSU root -> current tree node representing that component
    std::vector<int> comp_node(n);
    for (int v = 0; v < n; ++v) comp_node[v] = v;

    // 4) sort MST edges by weight
    std::vector<int> ord(m);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int i, int j){
        return mst_w[i] < mst_w[j];
    });

    int cur = n; // next internal node id

    for (int idx : ord) {
        int u = mst_edges[idx].first;
        int v = mst_edges[idx].second;
        double w = mst_w[idx];

        int ru = find(u);
        int rv = find(v);
        if (ru == rv) continue; // already same component

        // Create internal node for merging ru and rv
        if ((int)T.size() <= cur) T.resize(cur + 1);

        // This is key: Look into the paper from Mendel and Har-Peled: Fast Construction of Nets in Low Dimensional Metrics and Their
        //Applications
        // They multiply the distance by (n-1) to ensure the ultrametric property holds
        T[cur].Gamma = w * (double)(n-1);
        T[cur].child.clear();

        int nodeU = comp_node[ru];
        int nodeV = comp_node[rv];

        T[cur].child.push_back(nodeU);
        T[cur].child.push_back(nodeV);
        T[nodeU].parent = cur;
        T[nodeV].parent = cur;

        // Union DSU roots and set comp_node[new_root] = cur
        int new_root = unite_roots(ru, rv);
        comp_node[new_root] = cur;

        ++cur;
    }

    N = cur;

    // 5) Root is the tree node representing the final DSU component
    int global_root = find(0);
    root = comp_node[global_root];

    // 6) preprocess lifting
    preprocessLifting();
}

    void UltrametricTree::preprocessLifting() {
        gamma.assign(N, 0.0);
        for (int i = 0; i < N; ++i) gamma[i] = T[i].Gamma;

        depth.assign(N, 0);
        std::queue<int> q;
        if (root != -1) q.push(root), depth[root] = 0;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : T[u].child) {
                depth[v] = depth[u] + 1;
                q.push(v);
            }
        }
        int LOG = 1;
        while ((1<<LOG) <= std::max(1, N)) ++LOG;
        up.assign(LOG, std::vector<int>(N, -1));
        for (int u = 0; u < N; ++u) {
            up[0][u] = T[u].parent;
        }
        for (int k = 1; k < LOG; ++k) {
            for (int u = 0; u < N; ++u) {
                int mid = up[k-1][u];
                up[k][u] = (mid == -1 ? -1 : up[k-1][mid]);
            }
        }
    }

    // Return the *highest* ancestor of leaf v whose Γ <= threshold.
    // If leaf’s Γ(leaf)=0 is already > threshold (never happens), it returns the leaf itself.
    int UltrametricTree::sigmaDelta(int v_leaf, double Delta) const {
        if (N == 0 || v_leaf < 0 || v_leaf >= n) return v_leaf;
        const double thr = Delta / (2.0 * (double) n);
        // std::cout << "threshold: " << thr << std::endl;

        int u = v_leaf;
        if (gamma[u] > thr) return u; // safety
        // climb as long as parent exists and Γ(parent) <= thr
        for (int k = (int)up.size()-1; k >= 0; --k) {
            int p = up[k][u];
            if (p != -1 && gamma[p] <= thr) {
                u = p;
            }
        }
        return u;
    }