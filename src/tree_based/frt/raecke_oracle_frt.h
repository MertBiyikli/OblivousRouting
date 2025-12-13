//
// Created by Mert Biyikli on 11.12.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_ORACLE_FRT_H
#define OBLIVIOUSROUTING_RAECKE_ORACLE_FRT_H

#include "../raecke_oracle.h"
#include "../raecke_tree.h"
#include "efficient_frt.h"

class FRTOracle : public RaeckeOracle {
public:
    EfficientFRT frt_algo;
    explicit FRTOracle(IGraph& g) : RaeckeOracle(g), frt_algo(g) {}
    std::shared_ptr<ITreeNode> getTree(std::vector<double>& distances) override {
        // first update the edge distances
        frt_algo.updateEdgeDistances(distances);
        auto tree = frt_algo.getTree();

        return tree;
    }

};

#endif //OBLIVIOUSROUTING_RAECKE_ORACLE_FRT_H