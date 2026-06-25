//
// Created by Mert Biyikli on 24.06.26.
//

#include "algorithms/oblivious/oblivious_routing_runner.h"

IRoutingResult ObliviousSolverRunner::run(IGraph& graph,const Config& cfg,SolverType type) const {
    IRoutingResult result;
    result.type = type;
    result.graph_name = cfg.filename;
    result.solver_name = getSolverName(type);

    auto solverOpt = makeSolver(type, graph);
    if (!solverOpt) {
        std::cerr << "[ERROR] Failed to create solver of type "
                  << static_cast<int>(type) << '\n';

        result.status = ResultStatus::ERROR_INVALID_SOLVER;
        return result;
    }

    auto& solver = *solverOpt;

    const auto t0 = timeNow();
    result.scheme = solver->solve();
    result.solve_runtime_microseconds = duration(timeNow() - t0);
    result.total_runtime_microseconds = result.solve_runtime_microseconds;

    if (!result.scheme) {
        std::cerr << "[ERROR] Solver returned null routing scheme\n";
        result.status = ResultStatus::ERROR_INVALID_ROUTING_SCHEME;
        return result;
    }

    appendMetricsIfAvailable(solver, result);
    appendObjectiveIfAvailable(solver, *result.scheme, result);

    DemandEvaluator::evaluate(graph, result.scheme, cfg, result);

    if (result.status == ResultStatus::OK) {
        result.status = ResultStatus::OK;
    }

    return result;
}

