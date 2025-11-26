//
// Created by Mert Biyikli on 29.09.25.
//

#ifndef OBLIVIOUSROUTING_LINEAR_OBLIVIOUS_ROUTING_RATIO_H
#define OBLIVIOUSROUTING_LINEAR_OBLIVIOUS_ROUTING_RATIO_H

#include <map>
#include <vector>
#include <tuple>
#include "../../src/utils/hash.h"
#include "../../src/datastructures/graph.h"


class LinearObliviousRatio{
    bool debug = false;
    Graph g;
    // This is the oblivious routing strategy for which the oblivious ratio should be computed for
    std::unordered_map<std::pair<int,int>, std::unordered_map<std::pair<int,int>, double>> routing;
    // Worst case LP result
    double max_congestion = 0;
    std::unordered_map<std::pair<int, int>, double> demand_vector;

public:
    void init(Graph& g, const std::unordered_map<std::pair<int,int>, std::unordered_map<std::pair<int,int>, double>>& routing);
    double solve();

};

#endif //OBLIVIOUSROUTING_LINEAR_OBLIVIOUS_ROUTING_RATIO_H