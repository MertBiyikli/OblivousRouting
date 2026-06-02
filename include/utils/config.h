//
// Created by Mert Biyikli on 02.06.26.
//

#ifndef OBLIVIOUSROUTING_CONFIG_H
#define OBLIVIOUSROUTING_CONFIG_H

#include "../io/solver_io.h"



class RoutingResult {
    public:
    SolverType type;
    std::unique_ptr<RoutingScheme> scheme;
    double congestion;
    double oblivious_ratio;
    double total_runtime;
    int mwu_iterations;
};


class RoutingEngine
{
public:
    std::optional<RoutingResult> solve(
        IGraph& graph,
        const Config& cfg,
        const SolverType& type) {

        RoutingResult result;
        result.type = type;
        std::cout << "\n=== Running solver: " << getSolverName(type) << " ===\n";

        auto solver_opt = makeSolver(type, graph);
        if (!solver_opt) {
            std::cerr << "Failed to create solver of type " << static_cast<int>(type) << "\n";
        }
        auto& solver = *solver_opt;

        auto t0 = timeNow();
        result.scheme = solver->solve();
        auto t1 = timeNow();
        result.total_runtime = duration(t1-t0);

        std::cout << "Total time: " << result.total_runtime << " micro seconds\n";


        // Print time statistics if available
        if (auto mwu = dynamic_cast<MWUFramework*>(solver.get())) {
            mwu->printTimeStats();
        }

        // Compute oblivious ratio for linear schemes
        if (auto linear_scheme = dynamic_cast<LinearRoutingScheme*>(result.scheme.get())) {
            result.oblivious_ratio = linear_scheme->computeObliviousRatio();
            std::cout << "Oblivious ratio: " << result.oblivious_ratio << "\n";
        }

        // Evaluate demand models if provided
        if (cfg.evaluate_demand_models) {
            for (const auto& [model_name, offline_cong] : cfg.offline_opt_per_model) {
                auto it = cfg.demand_maps.find(model_name);
                if (it == cfg.demand_maps.end()) {
                    std::cerr << "Missing demand map for model: " << model_name << "\n";
                    return std::nullopt;
                }

                const demands& dmap = it->second;

                if (!result.scheme) {
                    std::cerr << "Solver returned null routing scheme\n";
                    return std::nullopt;
                }
                double scheme_cong = computeRoutingSchemeCongestion(graph, result.scheme, dmap);

                printStatsForDemandModel(model_name, {offline_cong, scheme_cong});
            }
        }
        return result;
    }
};

#endif //OBLIVIOUSROUTING_CONFIG_H