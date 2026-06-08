//
// Created by Mert Biyikli on 02.06.26.
//

#ifndef OBLIVIOUSROUTING_CONFIG_H
#define OBLIVIOUSROUTING_CONFIG_H

#include "../io/solver_io.h"
#include "../io/parse_argurment_io.h"

enum class OutPutFormat {
    TEXT,
    JASON
};

enum class ResultStatus {
    OK,
    ERROR_INVALID_SOLVER,
    ERROR_INVALID_ROUTING_SCHEME
};

class RoutingResult {
    public:
    ResultStatus status;
    SolverType type;
    std::unique_ptr<RoutingScheme> scheme;
    double congestion;
    double oblivious_ratio;
    double total_runtime;
    int mwu_iterations;

    void storeAsFile(const std::string& str, const OutPutFormat& format) {
        std::ofstream file;
        file.open(str);

        if (!file.is_open()) {
            std::cerr << "[ERROR] Failed to open file for writing: " << str << "\n";
            return;
        }

        if (!scheme->isValid()) {
            std::cerr << "[ERROR] Routing scheme is broken." << std::endl;
        }

        // store result in
        switch (format) {
            case OutPutFormat::TEXT:
                file << "Solver: " << getSolverName(type) << "\n";
                file << "Total runtime (micro seconds): " << total_runtime << "\n";
                file << "Oblivious ratio: " << oblivious_ratio << "\n";
                break;

            case OutPutFormat::JASON:
                file << "{\n";
                file << "  \"solver\": \"" << getSolverName(type) << "\",\n";
                file << "  \"total_runtime_microseconds\": " << total_runtime << ",\n";
                file << "  \"oblivious_ratio\": " << oblivious_ratio << "\n";
                file << "}\n";
                break;

            default:
                std::cerr << "[ERROR] Unknown output format.\n";
        }
        file.close();
    }
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
        //std::cout << "\n=== Running solver: " << getSolverName(type) << " ===\n";

        auto solver_opt = makeSolver(type, graph);
        if (!solver_opt) {
            std::cerr << "[ERROR] Failed to create solver of type " << static_cast<int>(type) << "\n";
            result.status = ResultStatus::ERROR_INVALID_SOLVER;
        }
        auto& solver = *solver_opt;
        auto t0 = timeNow();
        result.scheme = solver->solve();
        auto t1 = timeNow();
        result.total_runtime = duration(t1-t0);


        //result.scheme->printRoutingTable();


        // Print time statistics if available
        if (auto mwu = dynamic_cast<MWUFramework*>(solver.get())) {
            //mwu->printTimeStats();
        }

        // Compute oblivious ratio for linear schemes
        if (auto linear_scheme = dynamic_cast<LinearRoutingScheme*>(result.scheme.get())) {
            result.oblivious_ratio = linear_scheme->computeObliviousRatio();
        }

        // Evaluate demand models if provided
        if (cfg.evaluate_demand_models) {
            for (const auto& [model_name, offline_cong] : cfg.offline_opt_per_model) {

                auto it = cfg.demand_maps.find(model_name);
                if (it == cfg.demand_maps.end()) {
                    std::cerr << "[ERROR]: Missing demand for evaluating demand model. " << model_name << "\n";

                    return std::nullopt;
                }

                const demands& dmap = it->second;

                if (!result.scheme) {
                    std::cerr << "[ERROR]: Solver returned null routing scheme\n";
                    return std::nullopt;
                }
                result.congestion = computeRoutingSchemeCongestion(graph, result.scheme, dmap);
                //printStatsForDemandModel(model_name, {offline_cong, scheme_cong});
            }
        }
        result.status = ResultStatus::OK;
        return result;
    }
};

#endif //OBLIVIOUSROUTING_CONFIG_H