//
// Created by Mert Biyikli on 08.09.25.
//

#ifndef OBLIVIOUSROUTING_BINMOALMODEL_H
#define OBLIVIOUSROUTING_BINMOALMODEL_H

#include <vector>
#include "DemandModel.h"

class BimodalModel : public DemandModel {
public:
    bool debug = false;
    BimodalModel() = default;
    // Generate a gravity model demand matrix
    virtual DemandMap generate(IGraph& g, std::vector<std::pair<int, int>>& demands, double margin = 1.0);
};

#endif //OBLIVIOUSROUTING_BINMOALMODEL_H