//
// Created by Mert on 22.10.26.
//

#include "analysis/scenario_analysis.h"
#include "routing/routing_engine.h"

Result<void> solveInto(SolverType type, optimized::Graph<EdgeData>& working_graph, std::optional<IRoutingResult>& destination, double& seconds) {
    auto solver_opt = makeSolver(type, working_graph);
    if (!solver_opt || !(*solver_opt)) {
        return makeErrorMessage(ErrorCode::InvalidSolver, getSolverName(type) + " solver was not found.");
    }

    const auto started = std::chrono::steady_clock::now();
    auto scheme = (*solver_opt)->solve();
    seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();

    if (!scheme)
        return getError(scheme);
    if (!(*scheme)) {
        return makeErrorMessage(ErrorCode::SolverFailed, "The computed scheme is null.");
    }
    IRoutingResult routing;
    routing.scheme = std::move(scheme.value());
    destination.emplace(std::move(routing));
    return {};
}

static Result<NetworkData> evaluateFixedRoutingUpgrade(const NetworkData& baseline, const optimized::Graph<EdgeData>& upgraded_graph) {
    NetworkData result;
    result.total_demand_volume = baseline.total_demand_volume;
    result.links.reserve(baseline.links.size());

    const int m = upgraded_graph.getNumDirectedEdges();
    double utilization_sum = 0.0;

    for (const auto& link : baseline.links) {
        const int e = link.edge_id;

        if (e < 0 || e >= m) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid edge ID in baseline Data.");
        }

        const int reverse = upgraded_graph.reverse(e).id;

        if (reverse < 0 || reverse >= m || upgraded_graph.reverse(reverse).id != e) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid reverse edge mapping.");
        }

        const double forward_capacity = upgraded_graph.edgeData(e).capacity;

        const double reverse_capacity = upgraded_graph.edgeData(reverse).capacity;

        if (!std::isfinite(forward_capacity) || forward_capacity <= 0.0 || !std::isfinite(reverse_capacity) || reverse_capacity <= 0.0 || !std::isfinite(link.load) || link.load < 0.0) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid capacity or physical-link load.");
        }

        // Same routing and demands => identical physical load.
        const double load = link.load;

        const double utilization = std::max(load / forward_capacity, load / reverse_capacity);

        if (!std::isfinite(utilization)) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid upgraded link utilization.");
        }

        result.links.push_back({e, forward_capacity, load, utilization});

        result.max_congestion = std::max(result.max_congestion, utilization);

        utilization_sum += utilization;
    }

    if (!result.links.empty()) {
        result.average_congestion = utilization_sum / result.links.size();
    }

    return result;
}

Result<ScenarioResult> CapacityUpgrade::execute(const ScenarioContext& context, ExecutionMode mode) {
    if (context.fingerprint != fingerprint || optimized::GraphFingerprint::compute(context.graph) != fingerprint) {

        return makeErrorMessage(ErrorCode::InvalidArgument, "Scenario context does not match baseline graph.");
    }

    auto valid = verify();
    if (!valid) {
        return getError(valid);
    }

    // Construct upgraded graph before running solvers.
    auto modified = graph;

    auto applied = apply(modified);
    if (!applied) {
        return getError(applied);
    }

    std::cout << "\n[SCENARIO DIAGNOSTICS]\n";

    std::cout << "Original fingerprint: " << optimized::GraphFingerprint::compute(graph) << '\n';

    std::cout << "Modified fingerprint: " << optimized::GraphFingerprint::compute(modified) << '\n';

    for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
        const double before = graph.edgeData(e).capacity;
        const double after = modified.edgeData(e).capacity;

        if (before != after) {
            std::cout << "Edge " << e << " | Original capacity: " << before << " | Modified capacity: " << after << " | Factor: " << after / before << '\n';
        }
    }

    ScenarioResult output;
    output.kind = ScenarioKind::CAPACITY_UPGRADE;
    output.execution_mode = mode;
    output.original_fingerprint = fingerprint;
    output.modified_fingerprint = optimized::GraphFingerprint::compute(modified);

    for (const auto& [edge_id, factor] : upgrades) {
        output.capacity_upgrades.push_back({edge_id, graph.edgeData(edge_id).capacity, modified.edgeData(edge_id).capacity, factor});
    }

    std::sort(output.capacity_upgrades.begin(), output.capacity_upgrades.end(), [](const auto& a, const auto& b) { return a.edge_id < b.edge_id; });

    for (const SolverType type : getSolverBasedOnExecutionMode(mode)) {

        SolverAnalysisResult item;
        item.solver_algorithms.push_back(getSolverTypeName(type));

        item.baseline_graph = std::make_shared<optimized::Graph<EdgeData>>(graph);

        item.modified_graph = std::make_shared<optimized::Graph<EdgeData>>(modified);

        auto first = solveInto(type, *item.baseline_graph, item.baseline_routing, item.baseline_solver_seconds);

        if (!first) {
            return getError(first);
        }

        auto second = solveInto(type, *item.modified_graph, item.modified_routing, item.modified_solver_seconds);

        if (!second) {
            return getError(second);
        }

        if (!item.baseline_routing->scheme->isValid() || !item.modified_routing->scheme->isValid()) {
            return makeErrorMessage(ErrorCode::RuntimeError, "Routing scheme is not valid.");
        }

        output.solver_results.emplace_back(std::move(item));
    }

    return output;
}

Result<void> CapacityUpgrade::upgrade(const Config& cfg) {
    if (cfg.capacity_upgrades.empty()) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Capacity upgrade already empty.");
    }

    for (const auto& [edge_id, factor] : cfg.capacity_upgrades) {
        auto fail = upgrade(edge_id, factor);
        if (!fail) {
            return getError(fail);
        }
    }
    return {};
}

Result<void> CapacityUpgrade::upgrade(int edge_id, double scale) {
    auto valid = validEdge(edge_id);
    if (!valid)
        return valid;

    if (!std::isfinite(scale) || scale < 1.0) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Capacity upgrade factor must be finite and >= 1.");
    }
    const int canonical_edge_id = canonicalEdge(edge_id);
    if (!upgrades.emplace(canonical_edge_id, scale).second) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Capacity upgrade already exists for this physical link.");
    }
    return {};
}

Result<void> CapacityUpgrade::apply(optimized::Graph<EdgeData>& target) const {
    for (const auto& [edge_id, scale] : upgrades) {
        const int reverse_id = target.reverse(edge_id).id;
        auto forward = target.edgeData(edge_id);
        auto reverse = target.edgeData(reverse_id);
        // Validate BOTH representations before mutating this pair.
        const double upgraded_forward = forward.capacity * scale;
        const double upgraded_reverse = reverse.capacity * scale;
        if (!std::isfinite(upgraded_forward) || upgraded_forward <= 0.0 || !std::isfinite(upgraded_reverse) || upgraded_reverse <= 0.0) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Capacity upgrade produced invalid link capacity.");
        }
        forward.capacity = upgraded_forward;
        reverse.capacity = upgraded_reverse;
        target.setEdgeData(edge_id, forward);
        if (reverse_id != edge_id)
            target.setEdgeData(reverse_id, reverse);
    }
    return {};
}

Result<void> CapacityUpgrade::verify() {
    for (const auto& [edge_id, scale] : upgrades) {
        const int reverse_id = graph.reverse(edge_id).id;
        const auto forward = graph.edgeData(edge_id);
        const auto reverse = graph.edgeData(reverse_id);
        if (!std::isfinite(forward.capacity) || forward.capacity <= 0.0 || !std::isfinite(reverse.capacity) || reverse.capacity <= 0.0) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Capacity upgrade produced invalid link capacity.");
        }
        const double expected_forward = forward.capacity / scale;
        const double expected_reverse = reverse.capacity / scale;
        if (!std::isfinite(expected_forward) || expected_forward <= 0.0 || !std::isfinite(expected_reverse) || expected_reverse <= 0.0) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Capacity upgrade produced invalid link capacity.");
        }
    }
    return {};
}

Result<ScenarioResult> CapacityUpgrade::executeEvaluated(const ScenarioContext& context, ExecutionMode mode, const RoutingMetricsEvaluator& evaluator) {
    if (!evaluator) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "A metrics evaluator is required.");
    }

    // Reuse the same verified baseline/upgrade solver workflow.
    auto computed = execute(context, mode);
    if (!computed)
        return getError(computed);

    // Construct the upgraded graph again without touching scenario.graph.
    auto modified = graph;
    auto applied = apply(modified);
    if (!applied)
        return getError(applied);

    auto& result = computed.value();
    for (auto& solver_result : result.solver_results) {
        if (!solver_result.baseline_routing || !solver_result.modified_routing) {
            return makeErrorMessage(ErrorCode::SolverFailed, "Missing routing result for metric evaluation.");
        }
        // The same evaluator is used for both states. Callers must capture a
        // single, precomputed demand matrix (never regenerate it per graph).
        auto before = evaluator(*solver_result.baseline_graph, *solver_result.baseline_routing);

        auto after = evaluator(*solver_result.modified_graph, *solver_result.modified_routing);

        if (!before)
            return getError(before);

        if (!after)
            return getError(after);

        auto valid_Datas = [](const NetworkData& m) {
            return std::isfinite(m.max_congestion) && m.max_congestion >= 0.0 && std::isfinite(m.average_congestion) && m.average_congestion >= 0.0 && std::isfinite(m.total_demand_volume) &&
                   m.total_demand_volume >= 0.0;
        };
        if (!valid_Datas(before.value()) || !valid_Datas(after.value())) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Evaluator produced invalid network metrics.");
        }
        const double x = before.value().total_demand_volume;
        const double y = after.value().total_demand_volume;
        const double tolerance = 1e-9 * std::max({1.0, std::abs(x), std::abs(y)});
        if (std::abs(x - y) > tolerance) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Baseline and upgraded demand volumes differ. Use fixed demands.");
        }

        // 1. Evaluate the capacity upgrade without rerouting.
        // Reuse the baseline physical-link loads.
        auto fixed = evaluateFixedRoutingUpgrade(before.value(), *solver_result.modified_graph);

        if (!fixed) {
            return getError(fixed);
        }

        // 2. Baseline vs capacity upgrade with fixed routing.
        solver_result.fixed_routing_comparison = compareNetworkData(before.value(), fixed.value());

        // 3. Baseline vs capacity upgrade with recomputation.
        solver_result.comparison = compareNetworkData(before.value(), after.value());

        auto performance = deriveScenarioPerformanceMetrics(before.value(), fixed.value(), after.value());

        if (!performance)
            return getError(performance);

        solver_result.performance_metrics = std::move(performance.value());

        solver_result.baseline_Data = std::move(before.value());
        solver_result.fixed_routing_Data = std::move(fixed.value());
        solver_result.modified_Data = std::move(after.value());


    }

    return computed;
}

Result<ScenarioResult> CapacityUpgrade::executeWithDemands(const ScenarioContext& context, ExecutionMode mode, const demands& fixed_demands, double total_demand_volume) {
    if (!std::isfinite(total_demand_volume) || total_demand_volume < 0.0) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid demand volume.");
    }
    // Both evaluations use the exact same demand object; the solver's scheme
    // supplies its own graph-dependent routing behavior.

    auto evaluate = [&fixed_demands, total_demand_volume](const optimized::Graph<EdgeData>& graph, const IRoutingResult& routing) -> Result<NetworkData> {
        if (!routing.scheme) {
            return makeErrorMessage(ErrorCode::SolverFailed, "Routing scheme is null.");
        }

        std::vector<double> congestion;
        routing.scheme->routeDemands(congestion, fixed_demands);

        const int m = graph.getNumDirectedEdges();

        if (congestion.size() != static_cast<std::size_t>(m)) {
            return makeErrorMessage(ErrorCode::InvalidRouting, "Unexpected congestion vector size.");
        }

        NetworkData result;
        result.total_demand_volume = total_demand_volume;
        result.links.reserve(graph.getNumUndirectedEdges());

        double utilization_sum = 0.0;

        for (int e = 0; e < m; ++e) {
            const int reverse = graph.reverse(e).id;

            if (e > reverse) {
                continue;
            }

            const double forward_capacity = graph.edgeData(e).capacity;

            const double reverse_capacity = graph.edgeData(reverse).capacity;

            const double forward_utilization = congestion[e];

            const double reverse_utilization = congestion[reverse];

            if (!std::isfinite(forward_capacity) || !std::isfinite(reverse_capacity) || forward_capacity <= 0.0 || reverse_capacity <= 0.0 || !std::isfinite(forward_utilization) ||
                !std::isfinite(reverse_utilization) || forward_utilization < 0.0 || reverse_utilization < 0.0) {

                return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid physical-link metrics.");
            }

            const double forward_load = forward_utilization * forward_capacity;

            const double reverse_load = reverse_utilization * reverse_capacity;

            // Both directed representations must describe
            // the same physical-link load.
            const double tolerance = 1e-8 * std::max({1.0, std::abs(forward_load), std::abs(reverse_load)});

            if (std::abs(forward_load - reverse_load) > tolerance) {

                return makeErrorMessage(ErrorCode::InvalidRouting, "Inconsistent reverse-edge loads.");
            }

            const double load = std::max(forward_load, reverse_load);

            const double utilization = std::max(forward_utilization, reverse_utilization);

            result.links.push_back({e, forward_capacity, load, utilization});

            result.max_congestion = std::max(result.max_congestion, utilization);

            utilization_sum += utilization;
        }

        if (!result.links.empty()) {
            result.average_congestion = utilization_sum / result.links.size();
        }

        return result;
    };

    return executeEvaluated(context, mode, evaluate);
}

Result<void> LinkFailure::verify() {
    std::unordered_set<int> seen;
    for (const auto& edge_id : failed_links) {
        if (!seen.insert(canonicalEdge(edge_id)).second) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Duplicate physical link failure");
        }
    }

    return {};
}

Result<void> LinkFailure::linkFailure(int edge_id) {
    auto valid = validEdge(edge_id);
    if (!valid) {
        return valid;
    }

    const int canonical_edge_id = canonicalEdge(edge_id);
    failed_links.insert(canonical_edge_id);
    return {};
}

Result<void> TrafficGrowth::setGlobalScale(double scale) {
    if (scale < 0 || !std::isfinite(scale)) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Global scale must be finite and non-negative");
    }

    globalScale = scale;
    return {};
}

Result<void> TrafficGrowth::growPerCommodity(int source, int target, double scale) {
    if (!validNode(source) || !validNode(target)) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid commodity endpoint");
    }

    if (source == target) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Commodity endpoints must differ");
    }

    if ((scale < 0 || !std::isfinite(scale))) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Commodity scale must be finite and positive");
    }

    if (!commodities.insert({{source, target}, scale}).second) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Commodity scale already configured");
    }
    return {};
}

Result<void> TrafficGrowth::verify() {
    if ((globalScale < 0 || !std::isfinite(globalScale))) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid global traffic scale");
    }

    // Detect duplicate ordered source-target pairs.
    std::set<std::pair<int, int>> seen;

    for (const auto& [source_target, scale] : commodities) {
        const auto& [source, target] = source_target;
        if (!validNode(source) || !validNode(target)) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid commodity endpoint");
        }

        if (source == target || !(scale > 0 && std::isfinite(scale))) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid commodity scale");
        }

        if (!seen.emplace(source, target).second) {
            return makeErrorMessage(ErrorCode::LogicError, "Duplicate commodity scale");
        }
    }

    return {};
}

Result<ScenarioResult> Scenario::execute(const ScenarioContext& context, ExecutionMode mode) {
    if (context.fingerprint != fingerprint || optimized::GraphFingerprint::compute(context.graph) != fingerprint) {

        return makeErrorMessage(ErrorCode::InvalidArgument, "Scenario context does not match scenario graph.");
    }

    ScenarioResult output;
    output.execution_mode = mode;
    output.original_fingerprint = fingerprint;
    output.modified_fingerprint = fingerprint;

    for (const SolverType type : getSolverBasedOnExecutionMode(mode)) {

        SolverAnalysisResult item;
        item.solver_algorithms.emplace_back(getSolverTypeName(type));
        // Keep graph alive as long as its routing scheme exists.
        item.modified_graph = std::make_shared<optimized::Graph<EdgeData>>(graph);

        auto solved = solveInto(type, *item.modified_graph, item.modified_routing, item.modified_solver_seconds);

        if (!solved) {
            return getError(solved);
        }

        output.solver_results.emplace_back(std::move(item));
    }

    return output;
}
