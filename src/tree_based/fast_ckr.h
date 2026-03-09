//
// Created by Mert Biyikli on 09.02.26.
//

#ifndef OBLIVIOUSROUTING_FAST_CKR_H
#define OBLIVIOUSROUTING_FAST_CKR_H

#include "tree_oracle.h"
#include "../utils/priority_queue.h"

template<typename T>
class FastCKR : public TreeOracle<T> {
public:
    explicit FastCKR(IGraph& g) : TreeOracle<T>(g) {}
    explicit FastCKR(IGraph& g, bool mendelscaling) : TreeOracle<T>(g, mendelscaling) {}

    void computeLevelPartition(IGraph& _g, HSTLevel& level, const std::vector<int>& x_perm, double delta) override {
        const int n = _g.getNumNodes();

        // Sample R in [Δ/4, Δ/2] ---
        std::mt19937_64 gen(std::random_device{}());
        std::uniform_real_distribution<double> dist(delta / 4.0, delta / 2.0);
        const double R = dist(gen);

        level.R = R;
        level.owner.assign(n, -1);
        level.centers.clear();

        std::vector<double> estimated_distances(n, std::numeric_limits<double>::infinity());
        std::vector<int> P(n, 0);


        MinHeap<double, int> Q;
        Q.resize(x_perm.size() * 4);

        for (int i = 0; i<x_perm.size(); i++) {
            int source = x_perm[i];
            if (level.owner[source] != -1) continue; // already assigned

            // if not assigned, make it a new center
            level.centers.push_back(source);
            P[source] = -1; // center has no predecessor
            double& dist_source = estimated_distances[source];
            if (dist_source > 0.0) {
                dist_source = 0.0;
                Q.insertOrAdjustKey(source, 0.0);
            }

            while (!Q.empty()
                && Q.topKey() <= R) {
                int w = Q.top();
                double dist_w = Q.topKey();
                Q.deleteTop();

                if (level.owner[w] != -1) continue; // already assigned

                level.owner[w] = source; // assign owner

                // label assignment
                if (P[w]== 0) P[w] = i+1;

                if (_g.neighbors(w).size() == 0) continue;
                for (const auto& u : _g.neighbors(w)) {
                    if (u >= level.owner.size()) continue;
                    if (level.owner[u]!=-1) continue;

                    double new_dist = dist_w + _g.getEdgeDistance(w, u);
                    if (new_dist < estimated_distances[u]) {
                        estimated_distances[u] = new_dist;
                        Q.insertOrAdjustKey(u, new_dist);
                    }
                }
            }

        }

        assert(level.centers.size() <= level.owner.size());
    }
};

#endif //OBLIVIOUSROUTING_FAST_CKR_H