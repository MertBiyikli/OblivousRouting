//
// Created by Mert Biyikli on 24.06.26.
//

#include "algorithms/oblivious/oblivious_routing_runner.h"

Result<IRoutingResult> ObliviousSolverRunner::run(optimized::Graph<EdgeData>& graph,const Config& cfg,SolverType type) const {
    IRoutingResult result;
    result.type = type;
    result.graph_name = cfg.filename;
    result.nodes = graph.getNumNodes();
    result.edges = graph.getNumUndirectedEdges();
    result.solver_name = getSolverName(type);

    auto solverOpt = makeSolver(type, graph);
    if (!solverOpt
        || !(*solverOpt)) {
        result.status = ResultStatus::ERROR_INVALID_SOLVER;
        return makeErrorMessage(ErrorCode::InvalidSolver, getSolverName(type)+" solver was not found.");
    }

    auto& solver = *solverOpt;

    const auto t0 = timeNow();
    auto scheme = solver->solve();
    if (!scheme) {
        return getError(scheme);
    }
    if ( !(*scheme)) {
        result.status = ResultStatus::ERROR_INVALID_ROUTING_SCHEME;
        return makeErrorMessage(ErrorCode::SolverFailed, "The computed scheme is null.");
    }else {
        result.scheme = std::move(scheme.value());
    }


    result.solve_runtime_microseconds = duration(timeNow() - t0);
    result.total_runtime_microseconds = result.solve_runtime_microseconds;


    appendMetricsIfAvailable(solver, result);
    appendObjectiveIfAvailable(solver, *result.scheme, result);

    if (!cfg.evaluate_demand_models) {
        return result;
    }

    auto evaluate = DemandEvaluator::evaluate(graph, result.scheme, cfg, result);
    if (!evaluate) {
        return getError(evaluate);
    }

    return result;
}



