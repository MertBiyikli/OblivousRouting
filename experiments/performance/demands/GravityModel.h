//
// Created by Mert Biyikli on 12.08.25.
//

#ifndef OBLIVIOUSROUTING_GRAVITYMODEL_H
#define OBLIVIOUSROUTING_GRAVITYMODEL_H
#include "../../src/tree_based/raecke_transform.h"

class GravityModel {
public:
    bool debug = false;
    GravityModel() = default;
    // Generate a gravity model demand matrix
    DemandMap generate(Graph& g, std::vector<std::pair<int, int>>& demands, double margin = 1.0);
};

#endif //OBLIVIOUSROUTING_GRAVITYMODEL_H
