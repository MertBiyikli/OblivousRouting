//
// Created by Mert on 30.10.26.
//

#ifndef E_ROUTING_WORKSPACE_PREPARER_H
#define E_ROUTING_WORKSPACE_PREPARER_H
#include "analysis/workspace/workspace_configuration.h"
#include "analysis/scenario_engine.h"
#include "analysis/capacity_upgrade/capacity_upgrade.h"
#include "analysis/traffic_change/traffic_growth.h"

#include <cstdint>

struct PreparedWorkspace {
    demands original_demands;
    demands effective_demands;

    double original_demand_volume = 0.0;
    double effective_demand_volume = 0.0;

    std::uint64_t original_fingerprint = 0;
    std::uint64_t effective_fingerprint = 0;
};

class WorkspacePreparer {
public:
    // On success, working_graph contains the effective capacities.
    // On failure, working_graph must remain unchanged.
    //
    // Always call with a fresh graph copy, not a shared registered graph.
    static Result<PreparedWorkspace> prepare(
        optimized::Graph<EdgeData>& working_graph,
        const Config& baseline_config,
        const WorkspaceConfiguration& modifications
    );
};

#endif // E_ROUTING_WORKSPACE_PREPARER_H
