//
// Created by Mert on 29.10.26.
//

#ifndef E_ROUTING_RESILIENCE_ANALYZER_H
#define E_ROUTING_RESILIENCE_ANALYZER_H

#include "analysis/failures/link_failures.h"
#include "analysis/failures/network_connectivity.h"
#include "analysis/graph_rollback_guard.h"
#include "analysis/scenario_engine.h"
#include <algorithm>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

enum class FailureStatus { SUCCESS, DISCONNECTED, EXECUTION_FAILED };

struct FailureOutcome {
    int edge_id = -1; // Stable canonical physical edge ID
    int source = -1;
    int target = -1;
    FailureStatus status = FailureStatus::EXECUTION_FAILED;
    std::optional<NetworkData> reoptimized;
    std::optional<ComparisonData> comparison;
};

struct N1Result {
    std::uint64_t graph_fingerprint = 0;
    SolverType solver;
    NetworkData baseline;
    std::vector<FailureOutcome> outcomes; // All active physical edges, even unsuccessful ones
};

// Sequential first version. Never share one mutable graph across threads.
// Restrict to one solver so a failure can be attributed unambiguously.
inline Result<N1Result> analyzeN1(optimized::Graph<EdgeData>& graph, const Config& cfg, ExecutionMode mode = ExecutionMode::FAST) {
    const auto initial_fingerprint = optimized::GraphFingerprint::compute(graph);
    auto prepared_result = ScenarioEngine::prepare(graph, cfg, mode);

    if (!prepared_result)
        return getError(prepared_result);
    const PreparedAnalysis& prepared = prepared_result.value();

    if (prepared.solvers.size() != 1)
        return makeErrorMessage(ErrorCode::InvalidArgument, "N-1 v1 supports exactly one solver per analysis.");

    if (optimized::GraphFingerprint::compute(graph) != initial_fingerprint)
        return makeErrorMessage(ErrorCode::LogicError, "Graph changed during N-1 preparation.");

    N1Result result;
    result.graph_fingerprint = initial_fingerprint;
    result.solver = prepared.solvers.front().type;
    result.baseline = prepared.solvers.front().baseline;

    // Canonical directed IDs are even in this graph representation and remain stable.
    // Collect IDs and endpoints *before* applying any failure.
    const int edge_bound = graph.globalEdgeCount() * 2;
    for (int e = 0; e < edge_bound; e += 2) {
        if (!graph.edgeAlive(e))
            continue;
        const auto [source, target] = graph.getEdgeEndpoints(e);
        FailureOutcome outcome;
        outcome.edge_id = e;
        outcome.source = source;
        outcome.target = target;
        result.outcomes.push_back(std::move(outcome));
    }

    for (auto& outcome : result.outcomes) {
        LinkFailure failure;
        auto registered = failure.linkFailure(outcome.edge_id);
        if (!registered)
            return getError(registered);

        // Classify demand disconnection separately, using the same predicate as
        // ScenarioEngine, while ensuring even early errors restore the graph.
        bool disconnected = false;
        {
            GraphRollbackGuard guard(graph);
            auto applied = failure.apply(graph);
            if (!applied)
                return getError(applied);
            auto connected = verifyDemandConnectivity(graph, prepared.fixed_demands);
            disconnected = !connected;
        }

        if (optimized::GraphFingerprint::compute(graph) != initial_fingerprint)
            return makeErrorMessage(ErrorCode::LogicError, "Graph not restored after connectivity check.");

        if (disconnected) {
            outcome.status = FailureStatus::DISCONNECTED;
            continue;
        }

        auto executed = ScenarioEngine::execute(graph, failure, prepared);
        if (optimized::GraphFingerprint::compute(graph) != initial_fingerprint)
            return makeErrorMessage(ErrorCode::LogicError, "Graph not restored after scenario execution.");
        if (!executed) {
            outcome.status = FailureStatus::EXECUTION_FAILED;
            continue; // Preserve the event and finish evaluating the other links.
        }
        const auto& solver_results = executed.value().solver_results;
        if (solver_results.size() != 1 || !solver_results.front().modified_Data || !solver_results.front().comparison)
            return makeErrorMessage(ErrorCode::LogicError, "Incomplete successful N-1 scenario result.");

        outcome.status = FailureStatus::SUCCESS;
        outcome.reoptimized = *solver_results.front().modified_Data;
        outcome.comparison = *solver_results.front().comparison;
    }
    return result;
}

#endif // E_ROUTING_RESILIENCE_ANALYZER_H
