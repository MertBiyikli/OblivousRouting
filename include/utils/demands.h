//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_DEMANDMAP_H
#define OBLIVIOUSROUTING_DEMANDMAP_H

#include "../data_structures/graph/Igraph.h"
#include "../data_structures/graph/graph.h"
#include "hash.h"

#include <optional>
#include <unordered_map>
#include <vector>

class demands {
  public:
    std::vector<int> source, target;
    std::vector<double> demand_values;
    std::unordered_map<std::pair<int, int>, double, PairHash> demand_map;

    void addDemand(int s, int t, double demand);

    size_t size() const;
    std::pair<int, int> getDemandPair(size_t idx) const;
    double getDemandValue(size_t idx) const;
    std::optional<double> getDemandValue(int s, int t) const;
};

// base class for the demand models
class DemandModel {
  public:
    int seed = 42;
    DemandModel() = default;
    virtual ~DemandModel() = default;
    virtual Result<demands> generate(optimized::Graph<EdgeData>& g, std::vector<std::pair<int, int>>& _demands, double margin = 1.0) = 0;
};

struct CombinedDemand {
    int source;
    int target;
    double value;
};

inline std::vector<CombinedDemand> combineDemandPairs(const demands& demandMap, int numNodes) {
    std::unordered_map<std::uint64_t, double> aggregated;
    aggregated.reserve(demandMap.size());

    for (std::size_t i = 0; i < demandMap.size(); ++i) {
        const int s = demandMap.source[i];
        const int t = demandMap.target[i];
        const double value = demandMap.demand_values[i];

        if (s == t || value <= 0.0) {
            continue;
        }

        const int u = std::min(s, t);
        const int v = std::max(s, t);

        const auto key = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(u)) << 32) | static_cast<std::uint32_t>(v);

        aggregated[key] += value;
    }

    std::vector<CombinedDemand> result;
    result.reserve(aggregated.size());

    for (const auto& [key, value] : aggregated) {
        result.push_back({static_cast<int>(key >> 32), static_cast<int>(key & 0xFFFFFFFFULL), value});
    }

    return result;
}

class BimodalModel : public DemandModel {
  public:
    BimodalModel() = default;
    virtual Result<demands> generate(optimized::Graph<EdgeData>& g, std::vector<std::pair<int, int>>& demands, double margin = 1.0) override;
};

class UniformModel : public DemandModel {
  public:
    UniformModel() = default;
    virtual Result<demands> generate(optimized::Graph<EdgeData>& g, std::vector<std::pair<int, int>>& demands, double margin = 1.0) override;
};

class GravityModel : public DemandModel {
  public:
    GravityModel() = default;
    virtual Result<demands> generate(optimized::Graph<EdgeData>& g, std::vector<std::pair<int, int>>& demands, double margin = 1.0) override;
};

class GaussianModel : public DemandModel {
  public:
    GaussianModel() = default;
    virtual Result<demands> generate(optimized::Graph<EdgeData>& g, std::vector<std::pair<int, int>>& demands, double margin = 1.0) override;
};

#endif // OBLIVIOUSROUTING_DEMANDMAP_H