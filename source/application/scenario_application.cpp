//
// Created by Mert on 02.11.26.
//

#include "application/scenario_application.h"

#include "analysis/scenario_engine.h"
#include "analysis/scenario_result_json.h"
#include "analysis/workspace/workspace_executer.h"
#include "analysis/workspace/workspace_preparer.h"
#include "io/parse_argument_io.h"

#include <cstdint>
#include <iostream>
#include <utility>

int ScenarioApplication::run(const Command& cmd) {
    if (cmd.action != Action::AnalyzeScenario) {
        std::cerr << "Invalid command for ScenarioApplication.\n";
        return 1;
    }

    if (!cmd.scenario_request) {
        std::cerr << "Missing scenario analysis request.\n";
        return 1;
    }

    const auto& request = *cmd.scenario_request;

    // Preserve the existing configuration handling.
    Config cfg = cmd.config;

    // 1. Load the original network.
    auto loaded = load_graph_optimized(cfg);

    if (!loaded) {
        std::cerr << "Graph loading failed: " << getError(loaded).error().message << '\n';
        return 1;
    }

    optimized::Graph<EdgeData> graph = *loaded.value();

    // 2. Construct the effective workspace exactly once.
    auto workspace = WorkspacePreparer::prepare(graph, cfg, request.workspace);

    if (!workspace) {
        std::cerr << "Workspace preparation failed: " << getError(workspace).error().message << '\n';
        return 1;
    }

    // 3. Prepare the baseline using effective demands.
    auto prepared = ScenarioEngine::prepare(graph, cfg, request.execution_mode, workspace->effective_demands);

    if (!prepared) {
        std::cerr << "Analysis preparation failed: " << getError(prepared).error().message << '\n';
        return 1;
    }

    // 4. Execute one requested analysis.
    Result<ScenarioResult> result;
    switch (request.kind) {

    case AnalysisKind::LinkFailure: {
        if (!request.failed_edge_id) {
            std::cerr << "Missing link failure target.\n";
            return 1;
        }

        result = WorkspaceExecutor::analyzeLinkFailure(graph, prepared.value(), *request.failed_edge_id);
        break;
    }

    case AnalysisKind::CapacityComparison: {
        if (request.capacity_upgrade_events.empty()) {
            std::cerr << "Missing capacity upgrade events.\n";
            return 1;
        }

        result = WorkspaceExecutor::capacityUpgrade(graph, prepared.value(), request.capacity_upgrade_events);
        break;
    }

    default:
        std::cerr << "Unsupported scenario analysis kind.\n";
        return 1;
    }

    // 5. Verify restoration on success and failure.
    const auto workspace_fingerprint = workspace->effective_fingerprint;

    if (optimized::GraphFingerprint::compute(graph) != workspace_fingerprint) {
        std::cerr << "Scenario did not restore workspace graph.\n";
        return 1;
    }

    if (!result) {
        std::cerr << "Scenario execution failed: " << getError(result).error().message << '\n';
        return 1;
    }

    // 6. Attach provenance without modifying graph or demands.
    WorkspaceProvenance provenance;

    provenance.original_graph_fingerprint = workspace->original_fingerprint;
    provenance.effective_graph_fingerprint = workspace->effective_fingerprint;
    provenance.original_demand_volume = workspace->original_demand_volume;
    provenance.effective_demand_volume = workspace->effective_demand_volume;
    provenance.modifications = request.workspace;

    result->workspace_provenance = std::move(provenance);

    // 7. Serialize completed or disconnected outcomes.
    auto written = ScenarioResultJsonWriter::write(result.value(), cfg.output_filename);

    if (!written) {
        std::cerr << "Scenario report failed: " << getError(written).error().message << '\n';
        return 1;
    }

    std::cout << "Scenario report written to: " << cfg.output_filename << '\n';

    return 0;
}
