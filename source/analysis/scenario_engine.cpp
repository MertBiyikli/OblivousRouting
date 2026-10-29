//
// Created by Mert on 27.10.26.
//

#include "analysis/scenario_engine.h"
#include "analysis/capacity_upgrade/capacity_upgrade.h"
#include "analysis/commodity_flow_diagnostics.h"
#include "analysis/failures/link_failures.h"
#include "analysis/failures/network_connectivity.h"
#include "analysis/graph_rollback_guard.h"
#include "analysis/scenario_demands.h"

Result<ScenarioResult> ScenarioEngine::runScenario(optimized::Graph<EdgeData>& graph, const Scenario& scenario, const Config& cfg, ExecutionMode mode, double demand_scale) {
    using Graph = optimized::Graph<EdgeData>;

    // Stage 1: Validate scenario.
    auto valid = scenario.verify(graph);
    if (!valid)
        return getError(valid);

    if (!(scenario.kind() == ScenarioKind::CAPACITY_UPGRADE) && !(scenario.kind() == ScenarioKind::LINK_FAILURE)) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Unsupported scenario type.");
    }

    // Stage 2: Generate one fixed demand matrix.
    // demand_scale is retained in the signature for
    // compatibility; cfg.demand_scale is the source
    // of truth in computeTotalDemandVolume().
    (void)demand_scale;

    auto generated = computeTotalDemandVolume(graph, cfg);

    if (!generated)
        return getError(generated);

    auto& [fixed_demands, total_volume] = generated.value();

    ScenarioResult output;

    output.kind = scenario.kind();
    output.execution_mode = mode;

    output.original_fingerprint = optimized::GraphFingerprint::compute(graph);
    if (scenario.kind() == ScenarioKind::LINK_FAILURE) {

        const auto& failure = static_cast<const LinkFailure&>(scenario);

        for (const int edge_id : failure.failed_links) {

            const int canonical = std::min(edge_id, graph.reverse(edge_id).id);

            output.link_failures.push_back(canonical);
        }

        std::sort(output.link_failures.begin(), output.link_failures.end());

        output.link_failures.erase(std::unique(output.link_failures.begin(), output.link_failures.end()), output.link_failures.end());
    }

    if (scenario.kind() == ScenarioKind::CAPACITY_UPGRADE) {
        const auto& upgrade = static_cast<const CapacityUpgrade&>(scenario);

        for (const auto [edge_id, factor] : upgrade.upgrades) {

            const int canonical = std::min(edge_id, graph.reverse(edge_id).id);

            const double original_capacity = graph.edgeData(edge_id).capacity;
            const double upgraded_capacity = original_capacity * factor;

            output.capacity_upgrades.push_back({edge_id, original_capacity, upgraded_capacity, factor});
        }
    }

    // Stage 3: Execute each solver independently.
    for (const auto solver_type : getSolverFromConfigAndMode(cfg,mode)) {

        SolverAnalysisResult item;

        item.solver_algorithms.push_back(getSolverTypeName(solver_type));

        NetworkData baseline;

        // Baseline routing must be evaluated while
        // the original graph is still intact.
        {
            std::optional<IRoutingResult> baseline_routing;

            auto solved = solveInto(solver_type, graph, baseline_routing, item.baseline_solver_seconds);

            if (!solved)
                return getError(solved);

            if (!baseline_routing->scheme || !baseline_routing->scheme->isValid()) {
                return makeErrorMessage(ErrorCode::InvalidRouting, "Invalid baseline routing scheme.");
            }

            auto evaluated = evaluateRouting(graph, *baseline_routing, fixed_demands, total_volume);

            if (!evaluated)
                return getError(evaluated);

            baseline = std::move(evaluated.value());
        }
        const auto& solver = output.solver_results.front();

        // Stage 4: Apply modifications transactionally.
        {
            GraphRollbackGuard guard(graph);

            auto applied = scenario.apply(graph);

            if (!applied)
                return getError(applied);

            if (scenario.kind() == ScenarioKind::LINK_FAILURE) {

                auto connectivity = verifyDemandConnectivity(graph, fixed_demands);

                if (!connectivity)
                    return getError(connectivity);
            }

            output.modified_fingerprint = optimized::GraphFingerprint::compute(graph);

            // Stage 5: Fixed-routing evaluation.
            // Only valid for capacity-only modifications.
            std::optional<NetworkData> fixed;

            if (scenario.kind() == ScenarioKind::CAPACITY_UPGRADE) {

                auto evaluated = evaluateFixedRoutingUpgrade(baseline, graph);

                if (!evaluated)
                    return getError(evaluated);

                fixed.emplace(std::move(evaluated.value()));
            }

            NetworkData modified;

            // Stage 6: Solve modified topology/capacities.
            {
                std::optional<IRoutingResult> modified_routing;

                auto solved = solveInto(solver_type, graph, modified_routing, item.modified_solver_seconds);

                if (!solved || !modified_routing)
                    return getError(solved);

                if (!modified_routing->scheme || !modified_routing->scheme->isValid()) {
                    return makeErrorMessage(ErrorCode::InvalidRouting, "Invalid modified routing scheme.");
                }

                auto evaluated = evaluateRouting(graph, *modified_routing, fixed_demands, total_volume);

                if (!evaluated)
                    return getError(evaluated);

                modified = std::move(evaluated.value());
            }

            // Stage 7: Compare materialized results.

            // Valid for both capacity upgrades and link failures.
            item.comparison = compareNetworkData(baseline, modified);

            if (scenario.kind() == ScenarioKind::CAPACITY_UPGRADE) {

                if (!fixed) {
                    return makeErrorMessage(ErrorCode::LogicError, "Missing fixed-routing capacity evaluation.");
                }

                item.fixed_routing_comparison = compareNetworkData(baseline, *fixed);

                auto performance = deriveScenarioPerformanceMetrics(baseline, *fixed, modified);

                if (!performance)
                    return getError(performance);

                item.performance_metrics = std::move(performance.value());

                item.fixed_routing_Data = std::move(*fixed);
            }

            // For link failures, no fixed-routing result exists yet.
            // We also avoid deriveScenarioPerformanceMetrics(),
            // which currently requires matching physical-link sets.

            // These exist for both scenario types.
            item.baseline_Data = std::move(baseline);
            item.modified_Data = std::move(modified);

            // The routing objects are already destroyed.
            // The graph is restored when guard exits.
        }

        output.solver_results.push_back(std::move(item));
    }

    const auto& solver = output.solver_results.front();

    std::cout << "Baseline congestion: " << solver.baseline_Data->max_congestion << '\n'
              << "Modified congestion: " << solver.modified_Data->max_congestion << '\n'
              << "Original fingerprint: " << output.original_fingerprint << '\n'
              << "Modified fingerprint: " << output.modified_fingerprint << '\n';

    return output;
}

Result<std::pair<demands, double>> ScenarioEngine::computeTotalDemandVolume(const optimized::Graph<EdgeData>& graph, const Config& cfg) {
    if (cfg.demand_models.empty()) {
        return makeErrorMessage(ErrorCode::InvalidDemand, "Demand models are empty.");
    }
    double total_demand_volume = 0.0; // {total demand volume, max congestion
    // Generate demands once using the ORIGINAL graph.
    auto generated = generateScenarioDemands(graph, cfg.demand_models.front(), cfg.demand_scale);
    if (!generated) {
        return getError(generated);
    }

    auto _demands = generated.value();

    for (const double value : _demands.demand_values) {
        if (!std::isfinite(value) || value < 0.0) {
            return makeErrorMessage(ErrorCode::InvalidDemand, "Demand value must be non-negative.");
        }
        total_demand_volume += value;
    }

    return std::make_pair(_demands, total_demand_volume);
}