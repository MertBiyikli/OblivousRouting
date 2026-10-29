//
// Created by Mert on 29.10.26.
//

#ifndef E_ROUTING_PREPARED_ANALYSIS_H
#define E_ROUTING_PREPARED_ANALYSIS_H

#include "scenario_results.h"

// Read-only by convention once returned by ScenarioEngine::prepare().
// Never stores graph-bound routing schemes: only their evaluated metrics.
struct PreparedSolver {
    SolverType type;
    NetworkData baseline;
    double baseline_solver_seconds = 0.0;
};

struct PreparedAnalysis {
    std::uint64_t graph_fingerprint = 0;
    demands fixed_demands;
    double total_demand_volume = 0.0;
    ExecutionMode mode = ExecutionMode::FAST;
    std::vector<PreparedSolver> solvers;
};


#endif // E_ROUTING_PREPARED_ANALYSIS_H
