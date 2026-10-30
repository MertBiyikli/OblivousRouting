//
// Created by Mert on 22.10.26.
//

#include "analysis/fixed_routing_scenario.h"
#include <cmath>
#include <unordered_set>

#include "analysis/fixed_routing_scenario.h"

#include <algorithm>
#include <cmath>
#include <unordered_set>

Result<NetworkData> transformScenario(const NetworkData& baseline, double demandScale, const std::unordered_map<int, double>& capacityFactors) {
    if (!std::isfinite(demandScale) || demandScale <= 0.0) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Demand scale must be finite and positive.");
    }

    if (baseline.links.empty()) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Baseline links must not be empty.");
    }

    for (const auto& [edgeId, factor] : capacityFactors) {
        if (edgeId < 0 || !std::isfinite(factor) || factor < 1.0) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid capacity upgrade.");
        }
    }

    NetworkData result = baseline;
    result.max_congestion = 0.0;
    result.average_congestion = 0.0;

    double utilizationSum = 0.0;
    std::unordered_set<int> observed;

    for (auto& link : result.links) {
        if (!observed.insert(link.edge_id).second) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Duplicate physical edge ID.");
        }

        if (!std::isfinite(link.load) || !std::isfinite(link.capacity) || link.load < 0.0 || link.capacity <= 0.0) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid baseline link metrics.");
        }

        const auto it = capacityFactors.find(link.edge_id);
        const double factor = it == capacityFactors.end() ? 1.0 : it->second;

        link.load *= demandScale;
        link.capacity *= factor;
        link.utilization = link.load / link.capacity;

        if (!std::isfinite(link.load) || !std::isfinite(link.capacity) || !std::isfinite(link.utilization) || link.capacity <= 0.0) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Scenario transformation overflow.");
        }

        result.max_congestion = std::max(result.max_congestion, link.utilization);

        utilizationSum += link.utilization;
    }

    for (const auto& [edgeId, factor] : capacityFactors) {
        if (!observed.contains(edgeId)) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Capacity upgrade references unknown edge.");
        }
    }

    result.average_congestion = utilizationSum / static_cast<double>(result.links.size());

    result.total_demand_volume = baseline.total_demand_volume * demandScale;

    if (!std::isfinite(result.average_congestion) || !std::isfinite(result.total_demand_volume)) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Scenario aggregate overflow.");
    }

    return result;
}
