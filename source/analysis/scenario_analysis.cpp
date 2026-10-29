//
// Created by Mert on 22.10.26.
//

#include "analysis/scenario_analysis.h"
#include "routing/routing_engine.h"



/*
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


Result<ScenarioResult> Scenario::executeEvaluated(const ScenarioContext& context, ExecutionMode mode, const RoutingMetricsEvaluator& evaluator) {
    if (!evaluator) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "A metrics evaluator is required.");
    }

    // Reuse the same verified baseline/upgrade solver workflow.
    auto computed = execute(context, mode);
    if (!computed)
        return getError(computed);




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

Result<ScenarioResult> Scenario::executeWithDemands(const ScenarioContext& context, ExecutionMode mode, const demands& fixed_demands, double total_demand_volume) {
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
*/