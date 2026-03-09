//
// Created by Mert Biyikli on 11.02.26.
//

#ifndef OBLIVIOUSROUTING_RANDOM_MST_H
#define OBLIVIOUSROUTING_RANDOM_MST_H

#include "tree_oracle.h"
#include "mst.h"

template<typename T>
class TreeMST : public TreeOracle<T> {
public:
    explicit TreeMST(IGraph& g) : TreeOracle<T>(g) {}

    // Override getTree: bypass the level-partition path entirely and
    // build the MST-based Räcke tree directly, dispatching on T.
    T getTree(std::vector<double>& distances) override {
        this->updateDistances(distances);
        RandomMST mst_algo(this->graph);
        auto mst_edges = mst_algo.computeMST();

        if constexpr (std::is_same_v<T, std::shared_ptr<HSTNode>>) {
            return mst_algo.buildHST(mst_edges, 0);
        } else {
            // Flat FlatHST build
            std::vector<std::tuple<double,int,int>> keyed;
            keyed.reserve(mst_edges.size());
            for (auto [u, v] : mst_edges) {
                double w = this->graph.getEdgeDistance(this->graph.getEdgeId(u, v));
                keyed.emplace_back(w, u, v);
            }
            std::sort(keyed.begin(), keyed.end());

            const int N = this->n;
            HSTBuilder builder(N);
            std::vector<int> cluster(N);
            for (int v = 0; v < N; ++v) cluster[v] = builder.leafOf(v);
            DSU dsu(N);

            for (auto& [w, u, v] : keyed) {
                int pu = dsu.find(u);
                int pv = dsu.find(v);
                if (pu == pv) continue;
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
    }

    void computeLevelPartition(IGraph& /*g*/, HSTLevel& /*level*/,
                               const std::vector<int>& /*x_perm*/, double /*delta*/) override {
        // not used — getTree handles everything directly
    }
};

#endif //OBLIVIOUSROUTING_RANDOM_MST_H
