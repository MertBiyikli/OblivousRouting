//
// Created by Mert on 30.10.26.
//

#include "analysis/workspace/workspace_preparer.h"

#include "analysis/traffic_change/demand_vector.h"
#include "core/errors.h"

#include <utility>

Result<PreparedWorkspace> WorkspacePreparer::prepare(optimized::Graph<EdgeData>& working_graph, const Config& baseline_config, const WorkspaceConfiguration& modifications) {
    PreparedWorkspace prepared;

    prepared.original_fingerprint = optimized::GraphFingerprint::compute(working_graph);

    // 1. Generate demands from the ORIGINAL graph.
    Config demand_config = baseline_config;
    demand_config.demand_scale = 1.0;

    auto generated = ScenarioEngine::computeTotalDemandVolume(working_graph, demand_config);

    if (!generated)
        return getError(generated);

    auto normalized = DemandVector::fromDemands(generated->first, working_graph.getNumNodes());

    if (!normalized)
        return getError(normalized);

    prepared.original_demand_volume = normalized->totalVolume();

    auto original = normalized->toDemands();

    if (!original)
        return getError(original);

    prepared.original_demands = std::move(original.value());

    // 2. Build the requested traffic transformation.
    TrafficGrowth traffic;

    auto scaled = traffic.setGlobalScale(modifications.demand_scale);

    if (!scaled)
        return getError(scaled);

    for (const auto& item : modifications.commodity_overrides) {
        auto registered = traffic.growPerCommodity(item.source, item.target, item.factor);

        if (!registered)
            return getError(registered);
    }

    auto traffic_valid = traffic.verify(working_graph);

    if (!traffic_valid)
        return getError(traffic_valid);

    auto transformed = traffic.transform(normalized.value());

    if (!transformed)
        return getError(transformed);

    prepared.effective_demand_volume = transformed->totalVolume();

    auto effective = transformed->toDemands();

    if (!effective)
        return getError(effective);

    prepared.effective_demands = std::move(effective.value());

    // 3. Validate all capacity upgrades before changing the graph.
    CapacityUpgrade capacity;

    Config capacity_config = baseline_config;
    capacity_config.capacity_upgrades = modifications.capacity_upgrades;

    auto configured = capacity.configure(capacity_config);

    if (!configured)
        return getError(configured);

    auto capacity_valid = capacity.verify(working_graph);

    if (!capacity_valid)
        return getError(capacity_valid);

    // 4. Apply topology/capacity modifications transactionally.
    const int checkpoint = working_graph.checkpoint();

    auto applied = capacity.apply(working_graph);

    if (!applied) {
        auto restored = working_graph.rollback(checkpoint);

        if (!restored)
            return getError(restored);

        return getError(applied);
    }

    // 5. Publish the effective graph fingerprint.
    prepared.effective_fingerprint = optimized::GraphFingerprint::compute(working_graph);

    return prepared;
}
