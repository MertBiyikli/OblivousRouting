//
// Created by Mert Biyikli on 08.12.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_ORACLE_H
#define OBLIVIOUSROUTING_RAECKE_ORACLE_H
#include "raecke_oracle_iteration.h"
#include "../datastructures/IGraph.h"
#include "raecke_tree.h"

class RaeckeOracle {
public:
    const IGraph& graph;   // non-owning reference
    explicit RaeckeOracle(const IGraph& g) : graph(g) {}
    virtual ~RaeckeOracle() = default;

    virtual std::shared_ptr<ITreeNode> getTree(std::vector<double>& distances) = 0;
};

#endif //OBLIVIOUSROUTING_RAECKE_ORACLE_H