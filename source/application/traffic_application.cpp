//
// Created by Mert on 02.11.26.
//

#include "application/traffic_application.h"

#include "analysis/capacity_upgrade/capacity_upgrade.h"
#include "analysis/scenario_engine.h"
#include "analysis/traffic_change/traffic_growth.h"
#include "analysis/traffic_change/traffic_growth_result.h"
#include "analysis/workspace/workspace_preparer.h"
#include "io/parse_argument_io.h"
#include "routing/routing_result.h"
#include "routing/storage/routing_scheme_verifier.h"
#include "routing/storage/routing_table_importer.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <tuple>
#include <utility>

int TrafficApplication::run(const Command& command) {
    if (command.action != Action::AnalyzeTraffic) {
        std::cerr << "Invalid command for TrafficApplication.\n";
        return 1;
    }

    auto cfg = command.config;

    // 1. Load the network.
    auto loadedGraph = load_graph_optimized(cfg);

    if (!loadedGraph) {
        std::cerr << "Error loading graph: " << loadedGraph.error().message << '\n';
        return 1;
    }

    optimized::Graph<EdgeData> graph = *loadedGraph.value();

    // 2. Import the existing routing scheme.
    const auto import_start = std::chrono::steady_clock::now();
    auto imported = RoutingTableImporter::load(graph, cfg.routing_table_input_file);
    const auto import_end = std::chrono::steady_clock::now();

    std::cerr << "Routing table import: " << std::chrono::duration<double>(import_end - import_start).count() << " s\n";

    if (!imported) {
        std::cerr << "Error loading routing table: " << imported.error().message << '\n';
        return 1;
    }

    // 3. Verify the imported routing scheme.
    const auto verify_start = std::chrono::steady_clock::now();

    auto verified = RoutingSchemeVerifier::verify(graph, *imported.value());

    const auto verify_end = std::chrono::steady_clock::now();

    std::cerr << "Routing table verification: " << std::chrono::duration<double>(verify_end - verify_start).count() << " s\n";

    if (!verified) {
        std::cerr << "Error validating routing scheme: " << verified.error().message << '\n';
        return 1;
    }

    WorkspaceConfiguration modifications;

    modifications.capacity_upgrades = cfg.capacity_upgrades;
    modifications.demand_scale = cfg.demand_scale;

    for (const auto& item : command.commodity_scales) {
        modifications.commodity_overrides.push_back({item.source, item.target, item.factor});
    }

    auto workspace = WorkspacePreparer::prepare(graph, cfg, modifications);

    if (!workspace) {
        std::cerr << "Workspace preparation failed: " << getError(workspace).error().message << '\n';
        return 1;
    }

    // 4. Reconstruct the routing result metadata.
    IRoutingResult result{};

    result.status = ResultStatus::OK;
    result.type = cfg.solvers.front();
    result.solver_name = getSolverTypeName(result.type);
    result.graph_name = std::filesystem::path(cfg.filename).stem().string();

    result.nodes = graph.getNumNodes();
    result.edges = graph.getNumUndirectedEdges();

    result.total_runtime_microseconds = 0.0;
    result.solve_runtime_microseconds = 0.0;
    result.oblivious_ratio = -1.0;

    result.scheme = std::move(imported.value());

    // 5. Generate the baseline independently of requested growth.
    Config baselineConfig = cfg;
    baselineConfig.demand_scale = 1.0;

    auto generated = ScenarioEngine::computeTotalDemandVolume(graph, baselineConfig);

    if (!generated) {
        std::cerr << "Traffic demand generation error: " << getError(generated).error().message << '\n';
        return 1;
    }

    // 6. Configure global and per-commodity traffic growth.
    TrafficGrowth scenario;

    auto configured = scenario.setGlobalScale(cfg.demand_scale);

    if (!configured) {
        std::cerr << "Traffic configuration error: " << getError(configured).error().message << '\n';
        return 1;
    }

    for (const auto& override : command.commodity_scales) {
        auto added = scenario.growPerCommodity(override.source, override.target, override.factor);

        if (!added) {
            std::cerr << "Invalid commodity override (" << override.source << ", " << override.target << "): " << getError(added).error().message << '\n';

            return 1;
        }
    }

    // 7. Evaluate the modified demand against fixed routing.
    auto analysis = ScenarioEngine::evaluateTrafficGrowth(graph, result, workspace->original_demands, workspace->effective_demands);

    if (!analysis) {
        std::cerr << "Traffic analysis error: " << getError(analysis).error().message << '\n';
        return 1;
    }

    // 8. Preserve the original requested transformations.
    analysis->requested_global_scale = cfg.demand_scale;

    analysis->commodity_overrides.clear();
    analysis->commodity_overrides.reserve(command.commodity_scales.size());

    for (const auto& item : command.commodity_scales) {
        analysis->commodity_overrides.push_back({item.source, item.target, item.factor});
    }

    // Canonicalize the report order for reproducibility.
    std::sort(analysis->commodity_overrides.begin(), analysis->commodity_overrides.end(), [](const auto& a, const auto& b) { return std::tie(a.source, a.target) < std::tie(b.source, b.target); });

    // 9. Serialize the traffic analysis.
    auto written = TrafficGrowthResultJsonWriter::write(analysis.value(), cfg.output_filename);

    if (!written) {
        std::cerr << "Traffic report error: " << getError(written).error().message << '\n';
        return 1;
    }

    std::cout << "Traffic analysis report written to: " << cfg.output_filename << '\n';
    std::cout << "Baseline max congestion: " << analysis->baseline.max_congestion << '\n';
    std::cout << "Modified max congestion: " << analysis->modified.max_congestion << '\n';

    return 0;
}
