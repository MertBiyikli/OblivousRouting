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

    void computeLevelPartition(IGraph &g, HSTLevel &level, const std::vector<int> &x_perm, double delta) override {
        // not used in this oracle
        return;
    }
};

#endif //OBLIVIOUSROUTING_RANDOM_MST_H