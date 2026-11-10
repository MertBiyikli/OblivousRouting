//
// Created by Mert on 02.11.26.
//

#include "application/resilience_application.h"
#include "analysis/resilience/resilience_analyzer.h"
#include "analysis/resilience/resilience_result_json.h"
#include "analysis/scenario_engine.h"
#include "analysis/workspace/workspace_configuration.h"
#include "analysis/workspace/workspace_preparer.h"
#include "io/parse_argument_io.h"

#include <iostream>

int ResilienceApplication::run(const Command& command) {
    if (command.action != Action::AnalyzeN1) {
        std::cerr << "Invalid command for ResilienceApplication.\n";
        return 1;
    }

    auto cfg = command.config;

    // 1. Load the original network.
    auto loadedGraph = load_graph_optimized(cfg);

    if (!loadedGraph) {
        std::cerr << "Error loading N-1 graph: " << loadedGraph.error().message << '\n';
        return 1;
    }

    optimized::Graph<EdgeData> graph = *loadedGraph.value();

    // 2. Construct the effective workspace from the original graph.
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

    // 3. Compute the baseline on the effective workspace.
    auto prepared = ScenarioEngine::prepare(graph, cfg, command.execution_mode, workspace->effective_demands);

    if (!prepared) {
        std::cerr << "N-1 preparation failed: " << getError(prepared).error().message << '\n';
        return 1;
    }

    // 4. Evaluate all physical-link failures independently.
    auto analysis = analyzeN1(graph, prepared.value());

    if (optimized::GraphFingerprint::compute(graph) != workspace->effective_fingerprint) {
        std::cerr << "N-1 failed to restore workspace graph\n";
        return 1;
    }

    if (!analysis) {
        std::cerr << "N-1 execution failed: " << getError(analysis).error().message << '\n';
        return 1;
    }


    // 5. Persist the resilience report.
    auto written = ResilienceResultJsonWriter::write(analysis.value(), cfg.output_filename);

    if (!written) {
        std::cerr << "N-1 analysis error: " << getError(written).error().message << '\n';
        return 1;
    }

    std::cout << "N-1 report written to: " << cfg.output_filename << '\n';

    return 0;
}
