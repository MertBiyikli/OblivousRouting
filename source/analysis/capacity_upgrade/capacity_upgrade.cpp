//
// Created by Mert on 28.10.26.
//

#include "analysis/capacity_upgrade/capacity_upgrade.h"

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

Result<void> CapacityUpgrade::configure(const Config& cfg) {
    CapacityUpgrade pending;

    for (const auto& [edge_id, factor] : cfg.capacity_upgrades) {
        auto result = pending.upgrade(edge_id, factor);

        if (!result)
            return getError(result);
    }

    upgrades = std::move(pending.upgrades);
    return {};
}

Result<void> CapacityUpgrade::upgrade(int edge_id, double scale) {
    if (edge_id < 0 || !std::isfinite(scale) || scale < 1.0) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid capacity upgrade.");
    }

    if (!upgrades.emplace(edge_id, scale).second) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Duplicate capacity upgrade.");
    }

    return {};
}

Result<void> CapacityUpgrade::verify(const optimized::Graph<EdgeData>& graph) const {
    std::unordered_set<int> canonical_ids;

    const int max_id = graph.globalEdgeCount() * 2;

    for (const auto& [edge_id, scale] : upgrades) {
        if (edge_id < 0 || edge_id >= max_id || !graph.edgeAlive(edge_id)) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid or inactive capacity upgrade edge.");
        }

        const int canonical = std::min(edge_id, graph.reverse(edge_id).id);

        if (!canonical_ids.insert(canonical).second) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Duplicate physical capacity upgrade.");
        }

        const double capacity = graph.edgeData(edge_id).capacity;

        const double updated = capacity * scale;

        if (!std::isfinite(scale) || scale < 1.0 || !std::isfinite(updated) || updated <= 0.0) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Capacity upgrade produces invalid capacity.");
        }
    }

    return {};
}

Result<void> CapacityUpgrade::apply(optimized::Graph<EdgeData>& graph) const {
    for (const auto& [edge_id, scale] : upgrades) {
        auto data = graph.edgeData(edge_id);

        data.capacity *= scale;

        // Physical edge data is shared by both arcs.
        graph.setEdgeData(edge_id, data);
    }

    return {};
}