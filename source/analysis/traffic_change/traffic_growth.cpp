//
// Created by Mert on 28.10.26.
//

#include "analysis/traffic_change/traffic_growth.h"

Result<void> TrafficGrowth::setGlobalScale(double scale) {
    if (scale <= 0 || !std::isfinite(scale)) {
        return makeErrorMessage(ErrorCode::InvalidDemand, "Global scale must be finite and non-negative");
    }

    globalScale = scale;
    return {};
}

Result<void> TrafficGrowth::growPerCommodity(int source, int target, double scale) {
    if (source < 0 || target < 0 || source == target) {
        return makeErrorMessage(ErrorCode::InvalidDemand, "Invalid commodity endpoints");
    }

    if (!std::isfinite(scale) || scale <= 0.0) {
        return makeErrorMessage(ErrorCode::InvalidDemand, "Commodity scale must be finite and positive");
    }

    commodities[{source, target}] = scale;
    return {};
}

Result<void> TrafficGrowth::configure(const Config& cfg) {
    // Validate before changing the current configuration.
    if (!std::isfinite(cfg.demand_scale) || cfg.demand_scale <= 0.0) {
        return makeErrorMessage(ErrorCode::InvalidDemand, "Invalid configured demand scale");
    }

    globalScale = cfg.demand_scale;

    // Configuration replaces prior commodity overrides.
    // Explicit overrides can be added afterwards.
    commodities.clear();

    return {};
}

Result<void> TrafficGrowth::verify(const optimized::Graph<EdgeData>& graph) const {
    if (graph.getNumNodes() <= 0) {
        return makeErrorMessage(ErrorCode::InvalidGraph, "Traffic analysis requires a nonempty graph");
    }

    if (!std::isfinite(globalScale) || globalScale <= 0.0) {
        return makeErrorMessage(ErrorCode::InvalidDemand, "Invalid global demand scale");
    }

    for (const auto& [pair, factor] : commodities) {
        if (pair.first < 0 || pair.second < 0 || pair.first >= graph.getNumNodes() || pair.second >= graph.getNumNodes() || pair.first == pair.second || !std::isfinite(factor) || factor <= 0.0) {
            return makeErrorMessage(ErrorCode::InvalidDemand, "Invalid commodity configuration");
        }
    }

    return {};
}


Result<DemandVector> TrafficGrowth::transform(const DemandVector& baseline) const {
    auto scaled = baseline.scaled(globalScale);

    if (!scaled) {
        return getError(scaled);
    }

    DemandVector result = std::move(scaled.value());

    // Iterate over ordered baseline entries rather than the
    // unordered configuration map to ensure deterministic behavior.
    for (const auto& entry : baseline.entries()) {
        const auto it = commodities.find({entry.source, entry.target});

        if (it == commodities.end()) {
            continue;
        }

        auto modified = result.scaledCommodity(entry.source, entry.target, it->second);

        if (!modified) {
            return getError(modified);
        }

        result = std::move(modified.value());
    }

    // Reject overrides for commodities not present in the baseline.
    // They should never be silently ignored.
    for (const auto& [pair, factor] : commodities) {
        (void)factor;

        bool found = false;

        for (const auto& entry : baseline.entries()) {
            if (entry.source == pair.first && entry.target == pair.second) {
                found = true;
                break;
            }
        }

        if (!found) {
            return makeErrorMessage(ErrorCode::InvalidDemand, "Commodity override not found in baseline demands");
        }
    }

    return result;
}