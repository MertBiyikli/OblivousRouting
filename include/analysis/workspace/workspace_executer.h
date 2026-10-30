//
// Created by Mert on 30.10.26.
//

#ifndef E_ROUTING_WORKSPACE_EXECUTER_H
#define E_ROUTING_WORKSPACE_EXECUTER_H

#include "analysis/capacity_upgrade/capacity_upgrade.h"
#include "analysis/failures/link_failures.h"
#include "analysis/scenario_engine.h"

class WorkspaceExecutor {
  public:
    static Result<ScenarioResult> linkFailure(optimized::Graph<EdgeData>& graph, const PreparedAnalysis& prepared, int edge_id) {
        LinkFailure failure;

        auto configured = failure.linkFailure(edge_id);
        if (!configured)
            return getError(configured);

        auto verified = failure.verify(graph);
        if (!verified)
            return getError(verified);

        return ScenarioEngine::execute(graph, failure, prepared);
    }

    static Result<ScenarioResult> capacityUpgrade(optimized::Graph<EdgeData>& graph, const PreparedAnalysis& prepared, int edge_id, double factor) {
        CapacityUpgrade upgrade;

        auto configured = upgrade.upgrade(edge_id, factor);
        if (!configured)
            return getError(configured);

        auto verified = upgrade.verify(graph);
        if (!verified)
            return getError(verified);

        return ScenarioEngine::execute(graph, upgrade, prepared);
    }
};

#endif // E_ROUTING_WORKSPACE_EXECUTER_H
