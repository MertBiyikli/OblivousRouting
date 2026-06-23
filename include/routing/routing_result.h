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

struct RoutingRunResult {
    ResultStatus status = ResultStatus::OK;
    SolverType type{};

    std::string solver_name;
    std::string routing_base; // "Tree", "Electrical", or empty for normal solvers

    double total_runtime_microseconds = 0.0;
    double preprocessing_runtime_microseconds = 0.0;
    double solve_runtime_microseconds = 0.0;

    double oblivious_ratio = -1.0;
    int mwu_iterations = -1;

    std::size_t candidate_paths = 0;
    double average_paths_per_pair = 0.0;

    std::vector<DemandEvaluationResult> demand_evaluations;

    // Optional: only filled for normal solvers or last scheme if needed
    std::unique_ptr<RoutingScheme> scheme;

    // Semi-oblivious: one demand-specific scheme per demand model
    std::unordered_map<std::string, std::unique_ptr<RoutingScheme>> demand_schemes;
};


#endif //OBLIVIOUSROUTING_ROUTING_RESULT_H