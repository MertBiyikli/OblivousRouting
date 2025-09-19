//
// Created by Mert Biyikli on 18.09.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_MST_TRANSFORM_H
#define OBLIVIOUSROUTING_RAECKE_MST_TRANSFORM_H

#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <utility>
#include <iostream>
#include "../../utils/hash.h"
#include "raecke_random_mst.h"   // for MSTTree + RandomMST
#include "../raecke_transform_base.h"

class RaeckeMSTTransform : public RaeckeTransformBase<MSTTree>{
public:
    EdgeDemandMap& addTree(
            MSTTree& t,
            double         lambda,
            Graph&   graph)
    {
        normalizeOldSolutionBasedOnNewLambda(lambda);

        // For every pair (src,dst) of terminals, route on unique MST path
        auto vertices = graph.getVertices();
        for (int src : vertices) {
            for (int dst : vertices) {
                if (src == dst) continue;
                if (dst > src) continue; // avoid double-counting

                auto path = RandomMST::getMSTPath(src, dst, t.parent);
                for (size_t i=0; i+1<path.size(); ++i) {
                    int u = path[i], v = path[i+1];
                    std::pair<int,int> arcFwd = {u,v};
                    std::pair<int,int> arcRev = {v,u};
                    std::pair<int,int> demand = {src,dst};
                    std::pair<int,int> revDemand = {dst,src};

                    auto& fwdMap = arc2demand2cumulativeFraction[arcFwd];
                    auto& revMap = arc2demand2cumulativeFraction[arcRev];

                    fwdMap[demand] += lambda;
                    revMap[revDemand] += lambda;
                }
            }
        }

        return arc2demand2cumulativeFraction;
    }

private:
    void normalizeOldSolutionBasedOnNewLambda(double lambda) {
        for (auto& [arc, dmap] : arc2demand2cumulativeFraction) {
            for (auto& [d, frac] : dmap) {
                frac *= (1.0 - lambda);
            }
        }
    }

};
#endif //OBLIVIOUSROUTING_RAECKE_MST_TRANSFORM_H