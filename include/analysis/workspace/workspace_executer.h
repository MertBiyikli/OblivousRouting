//
// Created by Mert on 30.10.26.
//

#ifndef E_ROUTING_WORKSPACE_EXECUTER_H
#define E_ROUTING_WORKSPACE_EXECUTER_H

#include "analysis/analysis_request.h"
#include "analysis/capacity_upgrade/capacity_upgrade.h"
#include "analysis/failures/link_failures.h"
#include "analysis/failures/network_connectivity.h"
#include "analysis/graph_rollback_guard.h"
#include "analysis/scenario_engine.h"

class WorkspaceExecutor {
  public:
    static Result<ScenarioResult> linkFailure(optimized::Graph<EdgeData>& graph, const PreparedAnalysis& prepared, int edge_id) {
        LinkFailure failure;

        auto configured = failure.linkFailure(edge_id);
        if (!configured)
            return getError(configured);

        auto verified = failure.verify(graph);
        if (!verified)
            return getError(verified);

        return ScenarioEngine::execute(graph, failure, prepared);
    }

    static Result<ScenarioResult> capacityUpgrade(optimized::Graph<EdgeData>& graph, const PreparedAnalysis& prepared, const std::vector<CapacityUpgradeEvent>& events) {
        if (events.empty()) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "At least one capacity upgrade is required.");
        }

        CapacityUpgrade upgrade;

        for (const auto& event : events) {
            auto configured = upgrade.upgrade(event.edge_id, event.factor);

            if (!configured)
                return getError(configured);
        }

        auto verified = upgrade.verify(graph);
        if (!verified)
            return getError(verified);

        return ScenarioEngine::execute(graph, upgrade, prepared);
    }

    static Result<ScenarioResult> capacityUpgrade(optimized::Graph<EdgeData>& graph, const PreparedAnalysis& prepared, int edge_id, double factor) {
        return capacityUpgrade(graph, prepared, std::vector{CapacityUpgradeEvent{edge_id, factor}});
    }

    static Result<ScenarioResult> analyzeLinkFailure(optimized::Graph<EdgeData>& graph, const PreparedAnalysis& prepared, int edge_id) {
        const auto initial_fingerprint = optimized::GraphFingerprint::compute(graph);

        if (initial_fingerprint != prepared.graph_fingerprint) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Prepared analysis does not match workspace graph.");
        }

        // Validate the baseline demands and graph before
        // interpreting any modified-network connectivity error.
        auto baseline_connected = verifyDemandConnectivity(graph, prepared.fixed_demands);

        if (!baseline_connected)
            return getError(baseline_connected);

        if (prepared.solvers.empty()) {
            return makeErrorMessage(ErrorCode::InvalidSolver, "Prepared analysis has no baseline solver.");
        }

        LinkFailure failure;

        auto registered = failure.linkFailure(edge_id);
        if (!registered)
            return getError(registered);

        auto valid = failure.verify(graph);
        if (!valid)
            return getError(valid);

        bool disconnected = false;
        std::uint64_t modified_fingerprint = 0;

        {
            GraphRollbackGuard guard(graph);

            auto applied = failure.apply(graph);
            if (!applied)
                return getError(applied);

            modified_fingerprint = optimized::GraphFingerprint::compute(graph);

            auto connectivity = verifyDemandConnectivity(graph, prepared.fixed_demands);

            // The baseline demand vectors were already validated.
            // Removing an edge cannot invalidate their endpoints
            // or values, so this failure indicates disconnection.
            disconnected = !connectivity;
        }

        if (optimized::GraphFingerprint::compute(graph) != initial_fingerprint) {
            return makeErrorMessage(ErrorCode::LogicError, "Graph not restored after connectivity preflight.");
        }

        if (!disconnected) {
            return linkFailure(graph, prepared, edge_id);
        }

        ScenarioResult output;

        output.kind = ScenarioKind::LINK_FAILURE;
        output.execution_mode = prepared.mode;
        output.original_fingerprint = initial_fingerprint;
        output.modified_fingerprint = modified_fingerprint;

        output.status = ScenarioExecutionStatus::Disconnected;
        output.diagnostic = "Link failure disconnects nodes with positive traffic demand.";

        const int canonical_edge = std::min(edge_id, graph.reverse(edge_id).id);

        output.link_failures.push_back(canonical_edge);

        // Preserve all existing baseline solver results.
        // Do not invent modified congestion measurements.
        for (const auto& solver : prepared.solvers) {
            SolverAnalysisResult item;

            item.solver_algorithms.push_back(getSolverTypeName(solver.type));

            item.baseline_solver_seconds = solver.baseline_solver_seconds;

            item.baseline_Data = solver.baseline;

            output.solver_results.push_back(std::move(item));
        }

        return output;
    }
};

#endif // E_ROUTING_WORKSPACE_EXECUTER_H
