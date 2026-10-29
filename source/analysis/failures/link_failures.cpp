//
// Created by Mert on 28.10.26.
//

#include "analysis/failures/link_failures.h"

Result<void> LinkFailure::linkFailure(int edge_id) {
    if (edge_id < 0) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Edge ID must be non-negative.");
    }

    if (!failed_links.insert(edge_id).second) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Duplicate link failure.");
    }

    return {};
}

Result<void> LinkFailure::verify(const optimized::Graph<EdgeData>& graph) const {
    std::unordered_set<int> canonical_ids;

    const int max_id = graph.globalEdgeCount() * 2;

    for (const int edge_id : failed_links) {
        if (edge_id < 0 || edge_id >= max_id) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid link failure edge ID.");
        }

        if (!graph.edgeAlive(edge_id)) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Cannot fail an inactive link.");
        }

        const int canonical = std::min(edge_id, graph.reverse(edge_id).id);

        if (!canonical_ids.insert(canonical).second) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Duplicate physical link failure.");
        }
    }

    return {};
}

Result<void> LinkFailure::configure(const Config& cfg) {
    LinkFailure pending;

    for (const int edge_id : cfg.failed_links) {
        auto result = pending.linkFailure(edge_id);

        if (!result)
            return getError(result);
    }

    failed_links = std::move(pending.failed_links);
    return {};
}

Result<void> LinkFailure::apply(optimized::Graph<EdgeData>& graph) const {
    for (const int edge_id : failed_links) {
        if (!graph.removeEdge(edge_id)) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Failed to remove requested link.");
        }
    }

    return {};
}