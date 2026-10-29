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
#include "core/errors.h"

Result<PreparedAnalysis> ScenarioEngine::prepare(optimized::Graph<EdgeData>& graph, const Config& cfg, ExecutionMode mode) {
    auto generated = computeTotalDemandVolume(graph, cfg);
    if (!generated)
        return getError(generated);

    PreparedAnalysis prepared;
    prepared.graph_fingerprint = optimized::GraphFingerprint::compute(graph);
    prepared.fixed_demands = std::move(generated.value().first);
    prepared.total_demand_volume = generated.value().second;
    prepared.mode = mode;

    const auto solver_types = getSolverFromConfigAndMode(cfg, mode);
    if (solver_types.empty())
        return makeErrorMessage(ErrorCode::InvalidSolver, "No solvers selected for scenario analysis.");

    for (const auto type : solver_types) {
        std::optional<IRoutingResult> routing;
        double seconds = 0.0;
        auto solved = solveInto(type, graph, routing, seconds);
        if (!solved)
            return getError(solved);
        if (!routing || !routing->scheme || !routing->scheme->isValid())
            return makeErrorMessage(ErrorCode::InvalidRouting, "Invalid baseline routing scheme.");

        auto evaluated = evaluateRouting(graph, *routing, prepared.fixed_demands, prepared.total_demand_volume);
        if (!evaluated)
            return getError(evaluated);

        PreparedSolver item;
        item.type = type;
        item.baseline = std::move(evaluated.value());
        item.baseline_solver_seconds = seconds;
        prepared.solvers.push_back(std::move(item));
    }
    return prepared;
}

Result<ScenarioResult> ScenarioEngine::execute(optimized::Graph<EdgeData>& graph, const Scenario& scenario, const PreparedAnalysis& prepared) {
    if (optimized::GraphFingerprint::compute(graph) != prepared.graph_fingerprint)
        return makeErrorMessage(ErrorCode::InvalidArgument, "Prepared analysis does not match the current graph.");

    auto valid = scenario.verify(graph);
    if (!valid)
        return getError(valid);

    if (scenario.kind() != ScenarioKind::CAPACITY_UPGRADE && scenario.kind() != ScenarioKind::LINK_FAILURE)
        return makeErrorMessage(ErrorCode::InvalidArgument, "Unsupported scenario type.");
    if (prepared.solvers.empty())
        return makeErrorMessage(ErrorCode::InvalidSolver, "Prepared analysis has no solver baselines.");

    ScenarioResult output;
    output.kind = scenario.kind();
    output.execution_mode = prepared.mode;
    output.original_fingerprint = prepared.graph_fingerprint;

    if (scenario.kind() == ScenarioKind::LINK_FAILURE) {
        const auto& failure = static_cast<const LinkFailure&>(scenario);
        for (const int edge_id : failure.failed_links)
            output.link_failures.push_back(std::min(edge_id, graph.reverse(edge_id).id));
        std::sort(output.link_failures.begin(), output.link_failures.end());
        output.link_failures.erase(std::unique(output.link_failures.begin(), output.link_failures.end()), output.link_failures.end());
    } else {
        const auto& upgrade = static_cast<const CapacityUpgrade&>(scenario);
        for (const auto& [edge_id, factor] : upgrade.upgrades) {
            const int canonical = std::min(edge_id, graph.reverse(edge_id).id);
            const double capacity = graph.edgeData(canonical).capacity;
            output.capacity_upgrades.push_back({canonical, capacity, capacity * factor, factor});
        }
    }

    // Graph rollback covers all solver executions and every early return.
    GraphRollbackGuard guard(graph);
    auto applied = scenario.apply(graph);
    if (!applied)
        return getError(applied);

    if (scenario.kind() == ScenarioKind::LINK_FAILURE) {
        auto connected = verifyDemandConnectivity(graph, prepared.fixed_demands);
        if (!connected)
            return getError(connected);
    }
    output.modified_fingerprint = optimized::GraphFingerprint::compute(graph);

    for (const auto& prepared_solver : prepared.solvers) {
        SolverAnalysisResult item;
        item.solver_algorithms.push_back(getSolverTypeName(prepared_solver.type));
        item.baseline_solver_seconds = prepared_solver.baseline_solver_seconds;
        const NetworkData& baseline = prepared_solver.baseline;

        std::optional<NetworkData> fixed;
        if (scenario.kind() == ScenarioKind::CAPACITY_UPGRADE) {
            auto evaluated = evaluateFixedRoutingUpgrade(baseline, graph);
            if (!evaluated)
                return getError(evaluated);
            fixed.emplace(std::move(evaluated.value()));
        }

        std::optional<IRoutingResult> modified_routing;
        auto solved = solveInto(prepared_solver.type, graph, modified_routing, item.modified_solver_seconds);
        if (!solved)
            return getError(solved);
        if (!modified_routing || !modified_routing->scheme || !modified_routing->scheme->isValid())
            return makeErrorMessage(ErrorCode::InvalidRouting, "Invalid or missing modified routing scheme.");

        auto evaluated = evaluateRouting(graph, *modified_routing, prepared.fixed_demands, prepared.total_demand_volume);
        if (!evaluated)
            return getError(evaluated);
        NetworkData modified = std::move(evaluated.value());

        item.comparison = compareNetworkData(baseline, modified);
        if (fixed) {
            item.fixed_routing_comparison = compareNetworkData(baseline, *fixed);
            auto performance = deriveScenarioPerformanceMetrics(baseline, *fixed, modified);
            if (!performance)
                return getError(performance);
            item.performance_metrics = std::move(performance.value());
            item.fixed_routing_Data = std::move(*fixed);
        }
        item.baseline_Data = baseline;
        item.modified_Data = std::move(modified);
        output.solver_results.push_back(std::move(item));
    }
    return output;
}

Result<ScenarioResult> ScenarioEngine::runScenario(optimized::Graph<EdgeData>& graph, const Scenario& scenario, const Config& cfg, ExecutionMode mode, double demand_scale) {
    // cfg.demand_scale is the source of truth; preserve the historic signature.
    (void)demand_scale;
    auto prepared = prepare(graph, cfg, mode);
    if (!prepared)
        return getError(prepared);
    auto output = execute(graph, scenario, prepared.value());
    if (!output)
        return getError(output);
    const auto& first = output.value().solver_results.front();
    std::cout << "Baseline congestion: " << first.baseline_Data->max_congestion << '\n'
              << "Modified congestion: " << first.modified_Data->max_congestion << '\n'
              << "Original fingerprint: " << output.value().original_fingerprint << '\n'
              << "Modified fingerprint: " << output.value().modified_fingerprint << '\n';
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

Result<TrafficGrowthResult> ScenarioEngine::evaluateTrafficGrowth(const optimized::Graph<EdgeData>& graph, const IRoutingResult& routing, const demands& baseline_demands,
                                                                  const TrafficGrowth& scenario) {
    // 1. Validate graph and existing routing scheme.
    if (!routing.scheme) {
        return makeErrorMessage(ErrorCode::InvalidRouting, "Traffic analysis requires an existing routing scheme.");
    }

    auto verified = scenario.verify(graph);
    if (!verified) {
        return getError(verified);
    }

    const auto fingerprint = optimized::GraphFingerprint::compute(graph);

    // 2. Normalize baseline demands.
    auto baseline_vector = DemandVector::fromDemands(baseline_demands, graph.getNumNodes());

    if (!baseline_vector) {
        return getError(baseline_vector);
    }

    // 3. Apply traffic growth to the demands, not the graph.
    auto modified_vector = scenario.transform(baseline_vector.value());

    if (!modified_vector) {
        return getError(modified_vector);
    }

    auto baseline_engine_demands = baseline_vector->toDemands();

    if (!baseline_engine_demands) {
        return getError(baseline_engine_demands);
    }

    auto modified_engine_demands = modified_vector->toDemands();

    if (!modified_engine_demands) {
        return getError(modified_engine_demands);
    }

    // 4. Evaluate the SAME routing scheme under both demands.
    // No solver is constructed or executed.
    auto baseline = evaluateRouting(graph, routing, baseline_engine_demands.value(), baseline_vector->totalVolume());

    if (!baseline) {
        return getError(baseline);
    }

    auto modified = evaluateRouting(graph, routing, modified_engine_demands.value(), modified_vector->totalVolume());

    if (!modified) {
        return getError(modified);
    }

    // 5. Ensure that the graph has not changed.
    if (optimized::GraphFingerprint::compute(graph) != fingerprint) {
        return makeErrorMessage(ErrorCode::InvalidGraph, "Graph changed during traffic evaluation.");
    }

    // 6. Produce demand-aware comparison metrics.
    TrafficGrowthResult result;

    result.graph_fingerprint = fingerprint;
    result.baseline = std::move(baseline.value());
    result.modified = std::move(modified.value());

    const double original_volume = result.baseline.total_demand_volume;

    const double modified_volume = result.modified.total_demand_volume;

    result.absolute_demand_change = modified_volume - original_volume;

    if (original_volume > 0.0) {
        result.demand_scale = modified_volume / original_volume;

        result.relative_demand_change = result.absolute_demand_change / original_volume;
    }

    result.absolute_congestion_change = result.modified.max_congestion - result.baseline.max_congestion;

    if (result.baseline.max_congestion > 0.0) {
        result.relative_congestion_change = result.absolute_congestion_change / result.baseline.max_congestion;
    }

    return result;
}
