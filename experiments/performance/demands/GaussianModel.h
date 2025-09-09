//
// Created by Mert Biyikli on 08.09.25.
//

#ifndef OBLIVIOUSROUTING_GAUSSIANMODEL_H
#define OBLIVIOUSROUTING_GAUSSIANMODEL_H

#include <vector>
#include "DemandModel.h"

class GaussianModel : public DemandModel {
public:
    bool debug = false;
    GaussianModel() = default;
    // Generate a gravity model demand matrix
    virtual DemandMap generate(Graph& g, std::vector<std::pair<int, int>>& demands, double margin = 1.0);
};

#endif //OBLIVIOUSROUTING_GAUSSIANMODEL_H