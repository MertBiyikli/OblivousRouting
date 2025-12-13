//
// Created by Mert Biyikli on 08.12.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_ORACLE_CKR_H
#define OBLIVIOUSROUTING_RAECKE_ORACLE_CKR_H

#include "../raecke_oracle.h"
#include "../raecke_tree.h"
#include "optimized/efficient_oracle_ckr.h"

class CKROracle : public RaeckeOracle {
public:
    EfficientCKR ckr_algo;
    explicit CKROracle(IGraph& g) : RaeckeOracle(g), ckr_algo(g) {}
    std::shared_ptr<ITreeNode> getTree(std::vector<double>& distances) override {
        // first update the edge distances
        ckr_algo.updateEdgeDistances(distances);
        auto tree = ckr_algo.getTree();
        return tree;
    }
};


#endif //OBLIVIOUSROUTING_RAECKE_ORACLE_CKR_H