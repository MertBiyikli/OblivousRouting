//
// Created by Mert on 30.10.26.
//

#ifndef E_ROUTING_ANALYSIS_REQUEST_H
#define E_ROUTING_ANALYSIS_REQUEST_H
#include "analysis/scenario_results.h"
#include "analysis/workspace/workspace_configuration.h"

#include <vector>

enum class AnalysisKind {
    LinkFailure,
    CapacityComparison
};

struct CapacityUpgradeEvent {
    int edge_id = -1;
    double factor = 1.0;
};

// Using the class as a "stack" for keeping track of the changes made for the current analysis request
struct ScenarioAnalysisRequest {
    AnalysisKind kind = AnalysisKind::LinkFailure;
    WorkspaceConfiguration workspace;
    // Present only for analyses requiring one failure target.
    std::optional<int> failed_edge_id;
    ExecutionMode execution_mode = ExecutionMode::FAST;
    std::vector<CapacityUpgradeEvent> capacity_upgrade_events;
};

#endif // E_ROUTING_ANALYSIS_REQUEST_H
