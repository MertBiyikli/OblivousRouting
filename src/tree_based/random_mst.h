//
// Created by Mert Biyikli on 11.02.26.
//

#ifndef OBLIVIOUSROUTING_RANDOM_MST_H
#define OBLIVIOUSROUTING_RANDOM_MST_H

#include "tree_oracle.h"
#include "mst.h"

class TreeMST : public TreeOracle {
public:
    explicit TreeMST(IGraph& g):TreeOracle(g) {}

    std::shared_ptr<HSTNode> getTree(std::vector<double>& distances) override {
            updateDistances(distances);
            RandomMST mst_algo(graph);
            auto mst_edges  = mst_algo.computeMST();
            auto t = mst_algo.buildHST( mst_edges,0);
            return t;
    }

    // Flat representation: build the same Räcke/MST tree but pack into HST
    // using HSTBuilder directly — no level-partition path needed.
    HST getFlatTree(std::vector<double>& distances) override {
        updateDistances(distances);
        RandomMST mst_algo(graph);
        auto mst_edges = mst_algo.computeMST();

        // Build the keyed edge list sorted by distance (same order as buildHST)
        // We re-use mst_algo's internal keyed list after computeMST populated it.
        // Re-sort by weight ascending (Kruskal order).
        std::vector<std::tuple<double,int,int>> keyed;
        keyed.reserve(mst_edges.size());
        for (auto [u, v] : mst_edges) {
            double w = graph.getEdgeDistance(graph.getEdgeId(u, v));
            keyed.emplace_back(w, u, v);
        }
        std::sort(keyed.begin(), keyed.end());

        const int N = n;
        HSTBuilder builder(N);  // creates leaf nodes 0..N-1 with members={v}

        // cluster[v] holds the current builder-index representing the
        // DSU component whose representative is v.
        std::vector<int> cluster(N);
        for (int v = 0; v < N; ++v) cluster[v] = builder.leafOf(v);

        DSU dsu(N);

        for (auto& [w, u, v] : keyed) {
            int pu = dsu.find(u);
            int pv = dsu.find(v);
            if (pu == pv) continue;

            // Pick the center as the representative of the heavier component
            // (mirrors buildHST: parent->center = cluster[pu]->center)
            int center = builder.node(cluster[pu]).center;
            int parent_idx = builder.addNode(center);

            builder.attach(parent_idx, cluster[pu]);
            builder.attach(parent_idx, cluster[pv]);
            builder.sortMembers(parent_idx);

            dsu.unite(pu, pv);
            cluster[dsu.find(pu)] = parent_idx;
        }

        int root_idx = cluster[dsu.find(0)];
        return builder.finalise(root_idx);
    }

    void computeLevelPartition(IGraph &g, HSTLevel &level, const std::vector<int> &x_perm, double delta) override {
        // not used in this oracle
        return;
    }
};

#endif //OBLIVIOUSROUTING_RANDOM_MST_H