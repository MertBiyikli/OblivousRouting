//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_DEMANDMAP_H
#define OBLIVIOUSROUTING_DEMANDMAP_H

#include <vector>
#include "../data_structures/graph/Igraph.h"

class demands {
    public:
    std::vector<int> source, target;
    std::vector<double> demand_values;

    void addDemand(int s, int t, double demand);

    size_t size() const;
    std::pair<int, int> getDemandPair(size_t idx) const;
    double getDemandValue(size_t idx) const;
};

// base classe for the demand models
class DemandModel {
public:
    DemandModel() = default;
    virtual ~DemandModel() = default;
    virtual demands generate(IGraph& g, std::vector<std::pair<int, int>>& _demands, double margin = 1.0) = 0;
};


class BimodalModel : public DemandModel {
public:
    BimodalModel() = default;
    virtual demands generate(IGraph& g, std::vector<std::pair<int, int>>& demands, double margin = 1.0) override;
};

class UniformModel : public DemandModel {
public:
    UniformModel() = default;
    virtual demands generate(IGraph& g, std::vector<std::pair<int, int>>& demands, double margin = 1.0) override;
};

class GravityModel : public DemandModel {
public:
    GravityModel() = default;
    virtual demands generate(IGraph& g, std::vector<std::pair<int, int>>& demands, double margin = 1.0) override;
};

class GaussianModel : public DemandModel {
public:
    GaussianModel() = default;
    virtual demands generate(IGraph& g, std::vector<std::pair<int, int>>& demands, double margin = 1.0) override;
};


#endif //OBLIVIOUSROUTING_DEMANDMAP_H