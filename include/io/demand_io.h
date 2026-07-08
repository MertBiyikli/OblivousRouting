//
// Created by Mert Biyikli on 25.03.26.
//

#ifndef OBLIVIOUSROUTING_DEMAND_IO_H
#define OBLIVIOUSROUTING_DEMAND_IO_H

#include <memory>
#include <map>
#include <functional>
#include "core/types.h"
#include "core/config.h"
#include "../utils/demands.h"
#include "../algorithms/lp/lp_mcf.h"



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




// Generate all-pairs demand list
inline std::vector<std::pair<int,int>> generateAllDemandPairs(IGraph& g) {
    std::vector<std::pair<int,int>> result;
    result.reserve(static_cast<size_t>(g.getNumNodes()) * (g.getNumNodes() - 1));
    for (int v : g.getVertices())
        for (int u : g.getVertices())
            if (v != u) result.push_back({v, u});
    return result;
}

// Demand model handling
inline Result<void> HandleDemandModels(const std::optional<Config>& cfg, IGraph& g,
                               std::function<void(const std::string&, const demands&)> callback) {
    if (!cfg || cfg->demand_models.empty()) {
        return makeErrorMessage(ErrorCode::InvalidDemand, "No demand models specified in the configuration.");
    }

    auto pairs = generateAllDemandPairs(g);
    for (DemandModelType type : cfg->demand_models) {
        auto model = makeDemandModel(type);
        auto dmap = model->generate(g, pairs);
        if (!dmap) {
            return getError(dmap);
        }
        callback(demandModelName(type), dmap.value());
    }
    return {};
}

inline demands GetSingleDemandModel(const std::optional<Config>& cfg, IGraph& g) {
    if (!cfg || cfg->demand_models.empty()) return demands{};
    auto pairs = generateAllDemandPairs(g);
    auto model = makeDemandModel(cfg->demand_models.front());
    auto generated = model->generate(g, pairs);
    return generated ? generated.value() : demands{};
}


inline Result<void> offlineOptimal(std::unique_ptr<IGraph>& g, Config& cfg) {
    if (cfg.evaluate_demand_models) {
        auto handle_demand = HandleDemandModels(cfg, *g,
            [&](const std::string& model_name, const demands& dmap) {
                cfg.demand_maps[model_name] = dmap;
                auto offline_opt = computeOfflineOptimalCongestion(*g, dmap);
                // Void callback: only set the result if computation succeeded
                if (offline_opt) {
                    cfg.offline_opt_per_model[model_name] = offline_opt.value();
                }
                // If offline_opt fails, silently skip (don't propagate error from void lambda)
            });
        if (!handle_demand) {
            return getError(handle_demand);
        }
    }
    return {};
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





#endif //OBLIVIOUSROUTING_DEMAND_IO_H