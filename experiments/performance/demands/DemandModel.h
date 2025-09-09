//
// Created by Mert Biyikli on 08.09.25.
//

#ifndef OBLIVIOUSROUTING_DEMANDMODEL_H
#define OBLIVIOUSROUTING_DEMANDMODEL_H

#include "../../src/utils/hash.h"
#include "../../src/graph.h"

class DemandModel {
public:
    bool debug = false;
    DemandModel() = default;
    virtual ~DemandModel() = default;
    // Generate a gravity model demand matrix
    virtual DemandMap generate(Graph& g, std::vector<std::pair<int, int>>& demands, double margin = 1.0) = 0;
};


#endif //OBLIVIOUSROUTING_DEMANDMODEL_H