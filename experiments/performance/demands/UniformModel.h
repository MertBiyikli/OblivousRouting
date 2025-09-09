//
// Created by Mert Biyikli on 08.09.25.
//

#ifndef OBLIVIOUSROUTING_UNIFORMMODEL_H
#define OBLIVIOUSROUTING_UNIFORMMODEL_H

#include <vector>
#include "DemandModel.h"

class UniformModel : public DemandModel {
public:
    bool debug = false;
    UniformModel() = default;
    // Generate a gravity model demand matrix
    virtual DemandMap generate(Graph& g, std::vector<std::pair<int, int>>& demands, double margin = 1.0);
};

#endif //OBLIVIOUSROUTING_UNIFORMMODEL_H