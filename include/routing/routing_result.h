//
// Created by Mert Biyikli on 22.06.26.
//

#ifndef OBLIVIOUSROUTING_ROUTING_RESULT_H
#define OBLIVIOUSROUTING_ROUTING_RESULT_H

#include <unordered_map>
#include <vector>
#include <memory>
#include "../io/demand_io.h"
#include "routing_table.h"
#include "core/types.h"

enum class OutputFormat {
    TEXT,
    JSON
};

enum class ResultStatus {
    OK,
    ERROR_INVALID_SOLVER,
    ERROR_INVALID_ROUTING_SCHEME,
    ERROR_MISSING_DEMAND_MODELS,
    ERROR_EXCEPTION
};

struct DemandEvaluationResult {
    DemandModelType demand_type{};
    double congestion = -1.0;

    double runtime_microseconds = -1.0;
};

class MWUMetrics{
public:
    MWUMetrics() {
        iteration_count = 0;
        solve_time = 0;
        transformation_time = 0;
        mwu_weight_update_time = 0;
    }

    virtual ~MWUMetrics() = default;


    std::vector<double> oracle_running_times;
    int iteration_count;
    double solve_time;
    double transformation_time;
    double mwu_weight_update_time;
    double load_computation_time{};

    [[nodiscard]] bool empty() const {
        return iteration_count == 0
            && solve_time == 0.0
            && transformation_time == 0.0
            && mwu_weight_update_time == 0.0
            && load_computation_time == 0.0
            && oracle_running_times.empty();
    }

    [[nodiscard]] double averageOracleTime() const {
        if (oracle_running_times.empty()) {
            return -1.0;
        }

        double sum = 0.0;
        for (double t : oracle_running_times) {
            sum += t;
        }

        return sum / static_cast<double>(oracle_running_times.size());
    }

    [[nodiscard]] int getIterationCount() const {
        return iteration_count;
    }

};


struct IRoutingResult {

    ResultStatus status = ResultStatus::OK; // By default
    std::string solver_name;
    std::string routing_base;
    std::string graph_name;
    SolverType type;
    double oblivious_ratio = 0.0;


    // Runtime
    double total_runtime_microseconds = 0.0;
    double preprocessing_runtime_microseconds = 0.0;
    double solve_runtime_microseconds = 0.0;

    // Optional: only filled for normal solvers or last scheme if needed
    std::unique_ptr<RoutingScheme> scheme;
    std::size_t candidate_paths = 0;
    double average_paths_per_pair = 0.0;

    MWUMetrics mwu_metrics;

    std::vector<DemandEvaluationResult> demand_evaluations;
    // Semi-oblivious: one demand-specific scheme per demand model
    std::unordered_map<std::string, std::unique_ptr<RoutingScheme>> demand_schemes;

};

struct SemiObliviousRoutingResult :public IRoutingResult{

    DemandModelType demand_type{};
    double congestion = -1.0;
    std::string path_selection_strategy;
};




inline void printTimeStats(const MWUMetrics& metrics) {
    std::cout << "Solve time: " << metrics.solve_time << " micro seconds\n";
    std::cout << "Transformation time: " << metrics.transformation_time << " micro seconds\n";
    std::cout << "MWU iterations: " << metrics.iteration_count << "\n";
    std::cout << "MWU load computation: " << metrics.load_computation_time << " micro seconds\n";
    std::cout << "Average oracle time: " << metrics.averageOracleTime() << " micro seconds\n";
    std::cout << "Total MWU weight update time: " << metrics.mwu_weight_update_time << " micro seconds\n";
}


inline void printSemiObliviousResult(
    const SemiObliviousRoutingResult& r
) {
    std::cout << "Routing base: " << r.path_selection_strategy << std::endl;
    std::cout << "Demand [" << demandModelName(r.demand_type) << "]\n";
    std::cout << "  Congestion: " << r.congestion << '\n';
    std::cout << "  Runtime: " << r.total_runtime_microseconds << " us\n";
    std::cout << "  Candidate paths: " << r.candidate_paths << '\n';
    std::cout << "  Avg paths/pair: " << r.average_paths_per_pair << '\n';
}



#endif //OBLIVIOUSROUTING_ROUTING_RESULT_H