//
// Created by Mert Biyikli on 24.03.26.
//

#ifndef OBLIVIOUSROUTING_MST_ORACLE_H
#define OBLIVIOUSROUTING_MST_ORACLE_H

#include "../tree_oracle.h"
#include "../../../../../data_structures/union_find/union_find.h"



template<typename T>
class TreeMST : public TreeOracle<T> {
    int n;
    int root;
    std::vector<std::pair<int, int>> mst_edges;

public:
    explicit TreeMST(IGraph& g) : TreeOracle<T>(g) {
        n = g.getNumNodes();
        root =0; // by default
    }

    // Override getTree: bypass the level-partition path entirely and
    // build the MST-based decomposition tree directly, dispatching on T.
    T getTree(std::vector<double>& distances) override {
        this->updateDistances(distances);
        MST mst_algo(this->graph);
        mst_edges = mst_algo.computeMST();
        if constexpr (std::is_same_v<T, std::shared_ptr<HSTNode>>) {
            return buildPointerHST();
        } else {
            return buildFlatHST();
        }
    }

    std::shared_ptr<HSTNode> buildPointerHST() {
        std::vector<std::tuple<double,int,int>> keyed;
        keyed.reserve(mst_edges.size());
        for (auto [u, v] : mst_edges) {
            double w = this->graph.getEdgeDistance(this->graph.getEdgeId(u, v));
            keyed.emplace_back(w, u, v);
        }
        std::vector<std::shared_ptr<HSTNode>> cluster(n);
        for (int v = 0; v < n; ++v) {
            cluster[v] = std::make_shared<HSTNode>(v);
            cluster[v]->id = v;
            cluster[v]->members = {v};
            cluster[v]->center = v;
        }

        // build Hierarchical Decomposition Tree where each spanning tree edge induces the partitions
        DSU dsu(n);

        for (auto& [w, u, v] : keyed) {

            int pu = dsu.find(u);
            int pv = dsu.find(v);
            if (pu == pv) continue;

            auto parent = std::make_shared<HSTNode>(std::numeric_limits<int>::max());
            parent->children = {cluster[pu], cluster[pv]};
            cluster[pu]->parent = parent;
            cluster[pv]->parent = parent;

            // merge members
            parent->members = cluster[pu]->members;
            parent->members.insert(
                parent->members.end(),
                cluster[pv]->members.begin(),
                cluster[pv]->members.end()
            );

            dsu.unite(pu, pv);
            int new_rep = dsu.find(pu);
            parent->center = cluster[pu]->center; // arbitrary
            cluster[new_rep] = parent;
        }
        int rep = dsu.find(root);


        this->scales.push_back(calculatePointerHSTHeight(cluster[rep]));
        assert(cluster[dsu.find(root)] != nullptr);
        assert(cluster[dsu.find(root)]->children.size() > 0);
        return cluster[rep];
    }

    FlatHST buildFlatHST() {
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

        auto t = builder.finalise(root_idx);
        this->scales.push_back(calculateFlatHSTHeight(t));
        return t;
    }

    void computeLevelPartition(IGraph& /*g*/, HSTLevel& /*level*/,
                               const std::vector<int>& /*x_perm*/, double /*delta*/) override {
        // not used — getTree handles everything directly
    }
};

#endif //OBLIVIOUSROUTING_MST_ORACLE_H


