//
// Created by Mert Biyikli on 25.03.26.
//

#ifndef OBLIVIOUSROUTING_DEMAND_IO_H
#define OBLIVIOUSROUTING_DEMAND_IO_H

#include <memory>
#include <map>
#include <functional>
#include "../utils/demands.h"
#include "../algorithms/lp/lp_mcf.h"

enum class DemandModelType {
    GRAVITY,
    BIMODAL,
    GAUSSIAN,
    UNIFORM,
    NONE
};

static const std::map<std::string, DemandModelType> DEMAND_MAP{
    {"gravity", DemandModelType::GRAVITY}, {"gravity_model", DemandModelType::GRAVITY},
    {"bimodal", DemandModelType::BIMODAL}, {"bimodal_model", DemandModelType::BIMODAL},
    {"gaussian", DemandModelType::GAUSSIAN}, {"gaussian_model", DemandModelType::GAUSSIAN},
    {"uniform", DemandModelType::UNIFORM}, {"uniform_model", DemandModelType::UNIFORM}
};

using DemandModelFactory = std::function<std::unique_ptr<DemandModel>()>;

static const std::map<DemandModelType, std::pair<std::string, DemandModelFactory>> DEMAND_MODELS{
    { DemandModelType::GRAVITY,
      {"gravity", []() { return std::make_unique<GravityModel>(); }} },

    { DemandModelType::BIMODAL,
      {"bimodal", []() { return std::make_unique<BimodalModel>(); }} },

    { DemandModelType::GAUSSIAN,
      {"gaussian", []() { return std::make_unique<GaussianModel>(); }} },

    { DemandModelType::UNIFORM,
      {"uniform", []() { return std::make_unique<UniformModel>(); }} }
};

inline std::string demandModelName(DemandModelType type) {
    auto it = DEMAND_MODELS.find(type);
    return (it != DEMAND_MODELS.end()) ? it->second.first : "<unknown>";
}

inline std::unique_ptr<DemandModel> makeDemandModel(DemandModelType type) {
    auto it = DEMAND_MODELS.find(type);
    if (it == DEMAND_MODELS.end())
        throw std::runtime_error("Unknown demand model type.");
    return it->second.second();
}

inline void printStatsForDemandModel(const std::string& model_name,
                                     std::pair<double, double> result) {
    if (result.first > 0.0 && result.second > 0.0) {
        std::cout << "Ratio off the optimal offline solution ["
                  << model_name << "] demand model: "
                  << (result.second / result.first) * 100.0 << "% "
                  << "(" << result.first << " / " << result.second << ")\n";
    } else {
        std::cout << "Invalid congestion values for demand model [" << model_name << "].\n";
    }
}

inline double computeRoutingSchemeCongestion(IGraph& _g,
                                             const std::unique_ptr<RoutingScheme>& routing_scheme,
                                             const demands& demand_map) {
    std::vector<double> congestion_per_edge(_g.getNumDirectedEdges(), 0.0);
    routing_scheme->routeDemands(congestion_per_edge, demand_map);
    double max_cong = routing_scheme->getMaxCongestion(congestion_per_edge);
    for (const auto& cong : congestion_per_edge)
        if (cong > max_cong) max_cong = cong;
    return max_cong;
}

inline double computeOfflineOptimalCongestion(IGraph& _g, const demands& demand_map) {
    CMMF_Solver mccf(_g);
    mccf.AddDemandMap(demand_map);
    auto offline_scheme = mccf.solve();
    return mccf.getCongestionForPassedDemandMap();
}



#endif //OBLIVIOUSROUTING_DEMAND_IO_H