//
// Created by Mert Biyikli on 27.11.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_ORACLE_ITERATION_H
#define OBLIVIOUSROUTING_RAECKE_ORACLE_ITERATION_H


#include <memory>
#include <cassert>
#include "raecke_tree.h"
#include "ckr/ckr_tree_decomposer.h"


class OracleTreeIteration {
private:
    double lambda = 0.0;
    const std::shared_ptr<ITreeNode> tree;          // owns the tree of this level
    std::vector<double> distance;

public:
    OracleTreeIteration() = default;

    OracleTreeIteration(double lambda,
                        std::shared_ptr<ITreeNode> tree,
                        const std::vector<double>& distance)
        : lambda(lambda),
          tree((tree)),
          distance(distance)
    {}

    // ---- lambda ----
    void setLambda(double l) noexcept { lambda = l; }
    double getLambda() const noexcept { return lambda; }

    // ---- tree ----
    std::shared_ptr<ITreeNode> getTree() const {
        assert(tree != nullptr);
        return tree;
    }

    // ---- tree ----
    /*
    void setTree(std::shared_ptr<ITreeNode> t) noexcept {
        assert(t != nullptr);
        tree = t;
    }*/
/*
    const ITreeType*& getTree() const {
        assert(tree != nullptr);
        return tree;
    }*/

    std::shared_ptr<ITreeNode> getTree() {
        assert(tree != nullptr);
        return tree;
    }


    // ---- distance ----
    void setDistance(const std::vector<double>& d) noexcept {
        distance = d;
    }

    const std::vector<double>& getDistance() const {
        return distance;
    }

};




#endif //OBLIVIOUSROUTING_RAECKE_ORACLE_ITERATION_H