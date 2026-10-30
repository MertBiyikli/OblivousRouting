//
// Created by Mert Biyikli on 24.06.26.
//

#include "algorithms/oblivious/oblivious_routing_runner.h"

#include "analysis/failures/link_exposure_result_adapter.h"
#include "routing/demand_preparation.h"

Result<IRoutingResult> ObliviousSolverRunner::run(optimized::Graph<EdgeData>& graph, const Config& cfg, SolverType type) const {
    IRoutingResult result;
    result.type = type;
    result.graph_name = cfg.filename;
    result.nodes = graph.getNumNodes();
    result.edges = graph.getNumUndirectedEdges();
    result.solver_name = getSolverName(type);
    result.total_runtime_microseconds = 0;

    auto solverOpt = makeSolver(type, graph);
    if (!solverOpt || !(*solverOpt)) {
        result.status = ResultStatus::ERROR_INVALID_SOLVER;
        return makeErrorMessage(ErrorCode::InvalidSolver, getSolverName(type) + " solver was not found.");
    }

    auto& solver = *solverOpt;

    const auto t0 = timeNow();
    auto scheme = solver->solve();
    if (!scheme) {
        return getError(scheme);
    }
    if (!(*scheme)) {
        result.status = ResultStatus::ERROR_INVALID_ROUTING_SCHEME;
        return makeErrorMessage(ErrorCode::SolverFailed, "The computed scheme is null.");
    } else {
        result.scheme = std::move(scheme.value());
    }

    result.solve_runtime_microseconds = duration(timeNow() - t0);
    result.total_runtime_microseconds += result.solve_runtime_microseconds;
    result.total_runtime_microseconds += result.preprocessing_runtime_microseconds;

    appendMetricsIfAvailable(solver, result);
    appendObjectiveIfAvailable(solver, *result.scheme, result);

    if (!cfg.evaluate_demand_models) {
        return result;
    }

    // Existing solver construction remains unchanged.

    if (cfg.evaluate_demand_models) {
        for (const auto demandType : cfg.demand_models) {
            auto demand = prepareDemands(graph, demandType, cfg.demand_scale);

            if (!demand) {
                return getError(demand);
            }

            const auto start = std::chrono::steady_clock::now();

            auto evaluated = DemandEvaluator::evaluate(graph, result, demand.value());

            if (!evaluated) {
                return getError(evaluated);
            }

            const auto end = std::chrono::steady_clock::now();

            DemandEvaluationResult summary{};
            summary.demand_type = demandType;
            summary.congestion = evaluated->max_congestion;

            summary.runtime_microseconds =
                std::chrono::duration<double, std::micro>(
                    end - start
                ).count();

            summary.network_data = std::move(evaluated.value());

            if (cfg.failure_analysis) {
                const auto failure_start = std::chrono::steady_clock::now();

                auto exposure = LinkExposureAnalyzer::analyze(
                    graph,
                    *result.scheme,
                    demand.value()
                );

                if (!exposure)
                    return getError(exposure);

                const auto failure_end = std::chrono::steady_clock::now();

                const double runtime =
                    std::chrono::duration<double, std::micro>(
                        failure_end - failure_start
                    ).count();

                populateLinkExposureMetrics(
                    summary,
                    exposure.value(),
                    runtime
                );
            }

            result.demand_evaluations.push_back(std::move(summary));
        }
    }

    return result;
}
