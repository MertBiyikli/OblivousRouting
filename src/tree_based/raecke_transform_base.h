//
// Created by Mert Biyikli on 19.09.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_TRANSFORM_BASE_H
#define OBLIVIOUSROUTING_RAECKE_TRANSFORM_BASE_H

#include "../datastructures/graph.h"
#include <unordered_map>
#include <utility>
#include "../utils/hash.h"



template<typename Tree>
class RaeckeTransformBase {
protected:
    EdgeDemandMap arc2demand2cumulativeFraction;

    virtual void normalizeOldSolutionBasedOnNewLambda(double lambda) {
        for (auto& [arc, dmap] : arc2demand2cumulativeFraction) {
            for (auto& [d, frac] : dmap) {
                frac *= (1.0 - lambda);
            }
        }
    }

public:
    virtual ~RaeckeTransformBase() = default;
    const EdgeDemandMap& getRoutingTable() const { return arc2demand2cumulativeFraction; }

    virtual EdgeDemandMap& addTree(
        Tree& tree,   // generic pointer, concrete subclasses cast appropriately
        double lambda,
        Graph& g) = 0;
};

#endif //OBLIVIOUSROUTING_RAECKE_TRANSFORM_BASE_H