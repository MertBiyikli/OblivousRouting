//
// Created by Mert Biyikli on 12.12.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_ORACLE_RANDOM_H
#define OBLIVIOUSROUTING_RAECKE_ORACLE_RANDOM_H

#include "../raecke_oracle.h"
#include "../raecke_tree.h"
#include "mst.h"

class RandomOracle : public RaeckeOracle {
public:
    RandomMST mst_algo;
    explicit RandomOracle(IGraph& g) : RaeckeOracle(g), mst_algo(g) {}
    std::shared_ptr<ITreeNode> getTree(std::vector<double>& distances) override {
        mst_algo.updateEdgeDistances(distances);
        auto mst_edges  = mst_algo.build_mst();
        auto t = mst_algo.buildRaeckeTree( mst_edges,0);
        return t;

    }

};

#endif //OBLIVIOUSROUTING_RAECKE_ORACLE_RANDOM_H