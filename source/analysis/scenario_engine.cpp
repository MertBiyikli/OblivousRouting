//
// Created by Mert on 27.10.26.
//

#include "analysis/scenario_engine.h"
#include "analysis/commodity_flow_diagnostics.h"
#include "analysis/scenario_demands.h"


Result<ScenarioResult> ScenarioEngine::runScenario(optimized::Graph<EdgeData>& graph, Scenario& abstract_scenario, const Config& cfg, ExecutionMode mode, double demand_scale) {
    // Generate demands from the original network.
    auto demand_res = computeTotalDemandVolume(graph, cfg);

    if (!demand_res) {
        return getError(demand_res);
    }

    auto& [fixed_demands, total_demand_volume] = demand_res.value();

    // The scenario must describe the original graph.
    // Do not change graph capacities here.
    ScenarioContext context(graph);

    if (auto* capacity_upgrade = dynamic_cast<CapacityUpgrade*>(&abstract_scenario)) {

        // Register capacity changes only.
        auto registered = capacity_upgrade->upgrade(cfg);

        if (!registered) {
            return getError(registered);
        }

        auto baselineExport = ScenarioLPExporter::save(graph, fixed_demands, "results/performance/baseline_lp_input.json");

        if (!baselineExport)
            return getError(baselineExport);

        auto upgraded_graph = graph;

        auto applied = capacity_upgrade->apply(upgraded_graph);
        if (!applied)
            return getError(applied);

        auto upgradedExport = ScenarioLPExporter::save(upgraded_graph, fixed_demands, "results/performance/upgraded_lp_input.json");

        if (!upgradedExport)
            return getError(upgradedExport);

        // executeWithDemands() will internally:
        // 1. Solve the original graph.
        // 2. Copy and upgrade the graph.
        // 3. Solve the modified graph.
        // 4. Evaluate both with fixed demands.
        auto analysis = capacity_upgrade->executeWithDemands(context, mode, fixed_demands, total_demand_volume);

        if (!analysis) {
            return getError(analysis);
        }

        // Temporary debugging diagnostics.
        for (const auto& solver : analysis->solver_results) {
            if (solver.solver_algorithms.empty() || solver.solver_algorithms.front() != "electrical_sketching") {
                continue;
            }

            if (!solver.baseline_routing || !solver.modified_routing || !solver.baseline_routing->scheme || !solver.modified_routing->scheme) {
                continue;
            }

            analyzeCommodityFlows(*solver.baseline_routing->scheme, *solver.modified_routing->scheme, fixed_demands, 8, 15);

            const int source = 2;
            const int target = 6;

            std::cout << "\n[COMMODITY ROUTING COMPARISON]\n";
            std::cout << "Commodity: " << source << " -> " << target << '\n';

            for (int e = 0; e < solver.baseline_graph->getNumDirectedEdges(); ++e) {

                const int reverse = solver.baseline_graph->reverse(e).id;

                if (e > reverse)
                    continue;

                const double before = solver.baseline_routing->scheme->getFlow(e, source, target);

                const double after = solver.modified_routing->scheme->getFlow(e, source, target);

                std::cout << "Edge " << e << " | Before: " << before << " | After: " << after << " | Delta: " << after - before << '\n';

            }
            if (solver.performance_metrics) {
                const auto& m = *solver.performance_metrics;

                std::cout << "\n[PERFORMANCE METRICS]\n";

                if (m.capacity_only_improvement)
                    std::cout << "Capacity improvement: " << *m.capacity_only_improvement * 100.0 << "%\n";

                if (m.reoptimization_gain)
                    std::cout << "Additional rerouting gain: " << *m.reoptimization_gain * 100.0 << "%\n";

                if (m.total_absolute_load_redistribution)
                    std::cout << "Load redistribution: " << *m.total_absolute_load_redistribution << '\n';

                if (m.changed_physical_links)
                    std::cout << "Changed physical links: " << *m.changed_physical_links << '\n';
            }
        }

        return std::move(analysis.value());
    }

    return makeErrorMessage(ErrorCode::InvalidArgument, "Unsupported scenario type.");
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