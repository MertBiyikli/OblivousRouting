//
// Created by Mert on 28.10.26.
//

#include "analysis/traffic_change/traffic_growth.h"
/*
Result<void> TrafficGrowth::setGlobalScale(double scale) {
    if (scale < 0 || !std::isfinite(scale)) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Global scale must be finite and non-negative");
    }

    globalScale = scale;
    return {};
}

Result<void> TrafficGrowth::growPerCommodity(int source, int target, double scale) {
    if (!validNode(source) || !validNode(target)) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid commodity endpoint");
    }

    if (source == target) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Commodity endpoints must differ");
    }

    if ((scale < 0 || !std::isfinite(scale))) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Commodity scale must be finite and positive");
    }

    if (!commodities.insert({{source, target}, scale}).second) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Commodity scale already configured");
    }
    return {};
}

Result<void> TrafficGrowth::verify() {
    if ((globalScale < 0 || !std::isfinite(globalScale))) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid global traffic scale");
    }

    // Detect duplicate ordered source-target pairs.
    std::set<std::pair<int, int>> seen;

    for (const auto& [source_target, scale] : commodities) {
        const auto& [source, target] = source_target;
        if (!validNode(source) || !validNode(target)) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid commodity endpoint");
        }

        if (source == target || !(scale > 0 && std::isfinite(scale))) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid commodity scale");
        }

        if (!seen.emplace(source, target).second) {
            return makeErrorMessage(ErrorCode::LogicError, "Duplicate commodity scale");
        }
    }

    return {};
}*/