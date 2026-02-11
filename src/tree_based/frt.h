//
// Created by Mert Biyikli on 09.02.26.
//

#ifndef OBLIVIOUSROUTING_FRT_H
#define OBLIVIOUSROUTING_FRT_H

#include "tree_oracle.h"

class FRT : public TreeOracle {
public:
    explicit FRT(IGraph& g):TreeOracle(g) {}

    virtual void computeLevelPartition(IGraph& g, HSTLevel& level, const std::vector<int>& x_perm, double delta) override {

        level.owner.resize(g.getNumNodes());
        for (const auto& v : g.getVertices()) {
            level.owner[v] = v; // initialize owner to itself
        }
        level.R = delta;

        // compute for each node in the permutation the ball of radius delta around it and assign to the cluster of that node
        std::unordered_set<int> assigned;
        for (const auto& v : x_perm) {
            if (assigned.find(v) != assigned.end()) {
                continue;
            }

            // v is the new center for the next cluster
            level.centers.push_back(v);
            assigned.insert(v);

            for (const auto& u : g.getVertices()) {
                if (assigned.find(u) != assigned.end()) {
                    continue;
                }

                double dist = g.getShortestPathDistance(v, u);
                if (dist <= delta) {
                    // assign u to the ball of v
                    level.owner[u] = v;
                    assigned.insert(u);
                }
            }
        }
        assert(level.centers.size() <= level.owner.size());
    }
};

#endif //OBLIVIOUSROUTING_FRT_H