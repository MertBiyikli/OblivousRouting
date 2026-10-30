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

inline Result<N1Result> analyzeN1(optimized::Graph<EdgeData>& graph, const PreparedAnalysis& prepared) {
    const auto initial_fingerprint = optimized::GraphFingerprint::compute(graph);

    // The prepared analysis must correspond to the
    // effective workspace graph.
    if (prepared.graph_fingerprint != initial_fingerprint) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Prepared N-1 analysis does not match workspace graph.");
    }

    if (prepared.solvers.size() != 1) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "N-1 supports exactly one solver per analysis.");
    }

    N1Result result;
    result.graph_fingerprint = initial_fingerprint;
    result.solver = prepared.solvers.front().type;
    result.baseline = prepared.solvers.front().baseline;

    // Enumerate active physical links before executing failures.
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

    // Each failure is independent and must restore
    // the effective workspace graph.
    for (auto& outcome : result.outcomes) {
        LinkFailure failure;

        auto registered = failure.linkFailure(outcome.edge_id);
        if (!registered)
            return getError(registered);

        bool disconnected = false;

        {
            GraphRollbackGuard guard(graph);

            auto applied = failure.apply(graph);
            if (!applied)
                return getError(applied);

            auto connected = verifyDemandConnectivity(graph, prepared.fixed_demands);

            disconnected = !connected;
        }

        if (optimized::GraphFingerprint::compute(graph) != initial_fingerprint) {
            return makeErrorMessage(ErrorCode::LogicError, "Graph not restored after connectivity check.");
        }

        if (disconnected) {
            outcome.status = FailureStatus::DISCONNECTED;
            continue;
        }

        auto executed = ScenarioEngine::execute(graph, failure, prepared);

        if (optimized::GraphFingerprint::compute(graph) != initial_fingerprint) {
            return makeErrorMessage(ErrorCode::LogicError, "Graph not restored after N-1 execution.");
        }

        if (!executed) {
            outcome.status = FailureStatus::EXECUTION_FAILED;
            continue;
        }

        const auto& solvers = executed.value().solver_results;

        if (solvers.size() != 1 || !solvers.front().modified_Data || !solvers.front().comparison) {
            return makeErrorMessage(ErrorCode::LogicError, "Incomplete successful N-1 result.");
        }

        outcome.status = FailureStatus::SUCCESS;

        outcome.reoptimized = *solvers.front().modified_Data;

        outcome.comparison = *solvers.front().comparison;
    }

    return result;
}

// Sequential first version. Never share one mutable graph across threads.
// Restrict to one solver so a failure can be attributed unambiguously.
inline Result<N1Result> analyzeN1(optimized::Graph<EdgeData>& graph, const Config& cfg, ExecutionMode mode) {
    auto prepared = ScenarioEngine::prepare(graph, cfg, mode);

    if (!prepared)
        return getError(prepared);

    return analyzeN1(graph, prepared.value());
}

#endif // E_ROUTING_RESILIENCE_ANALYZER_H
