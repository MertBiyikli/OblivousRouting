//
// Created by Mert on 29.10.26.
//

#ifndef E_ROUTING_NETWORK_DATA_H
#define E_ROUTING_NETWORK_DATA_H
#include "core/errors.h"
#include <optional>
#include <cmath>
#include <vector>
#include <unordered_map>

struct LinkData {
    int edge_id = -1;
    double capacity = 0.0;
    double load = 0.0;
    double utilization = 0.0;
    int source = -1;
    int target = -1;
};

struct NetworkData {
    double max_congestion = 0.0;
    double average_congestion = 0.0;
    double total_demand_volume = 0.0;
    std::vector<LinkData> links;
};

// Derived capacity-planning metrics for one fixed demand matrix.
// All improvement values are fractions (0.15 = 15%); negative means worse.
// Undefined relative improvements (zero denominator) remain nullopt.
struct ScenarioPerformanceMetrics {
    double baseline_congestion = 0.0;
    double fixed_routing_congestion = 0.0;
    double reoptimized_congestion = 0.0;

    std::optional<double> capacity_only_improvement;
    std::optional<double> total_improvement;
    std::optional<double> reoptimization_gain;

    // Sum of absolute changes of PHYSICAL-link loads, in demand-volume units.
    // nullopt when per-link data is missing.
    std::optional<double> total_absolute_load_redistribution;
    std::optional<double> maximum_absolute_link_load_change;
    std::optional<std::size_t> changed_physical_links;
};

inline Result<ScenarioPerformanceMetrics> deriveScenarioPerformanceMetrics(const NetworkData& baseline, const NetworkData& fixed, const NetworkData& reoptimized) {
    auto valid = [](const NetworkData& network) {
        return std::isfinite(network.max_congestion) && network.max_congestion >= 0.0 && std::isfinite(network.total_demand_volume) && network.total_demand_volume >= 0.0;
    };

    if (!valid(baseline) || !valid(fixed) || !valid(reoptimized)) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Non-finite or negative scenario metrics.");
    }

    const double volume_scale = std::max({1.0, baseline.total_demand_volume, fixed.total_demand_volume, reoptimized.total_demand_volume});

    if (std::abs(baseline.total_demand_volume - fixed.total_demand_volume) > 1e-9 * volume_scale || std::abs(baseline.total_demand_volume - reoptimized.total_demand_volume) > 1e-9 * volume_scale) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Scenario metrics use different demand volumes.");
    }

    ScenarioPerformanceMetrics metrics;
    metrics.baseline_congestion = baseline.max_congestion;
    metrics.fixed_routing_congestion = fixed.max_congestion;
    metrics.reoptimized_congestion = reoptimized.max_congestion;

    if (baseline.max_congestion > 0.0) {
        metrics.capacity_only_improvement = 1.0 - fixed.max_congestion / baseline.max_congestion;

        metrics.total_improvement = 1.0 - reoptimized.max_congestion / baseline.max_congestion;
    }

    if (fixed.max_congestion > 0.0) {
        metrics.reoptimization_gain = 1.0 - reoptimized.max_congestion / fixed.max_congestion;
    }

    if (!baseline.links.empty() && !reoptimized.links.empty()) {
        std::unordered_map<int, double> loads;

        for (const auto& link : baseline.links) {
            if (link.edge_id < 0 || !std::isfinite(link.load) || link.load < 0.0 || !loads.emplace(link.edge_id, link.load).second) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid or duplicate baseline physical link.");
            }
        }

        if (loads.size() != reoptimized.links.size()) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Physical link sets differ across scenario evaluations.");
        }

        double total = 0.0;
        double maximum = 0.0;
        std::size_t changed = 0;

        for (const auto& link : reoptimized.links) {
            const auto it = loads.find(link.edge_id);

            if (it == loads.end() || !std::isfinite(link.load) || link.load < 0.0) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "Missing or invalid modified physical link.");
            }

            const double difference = std::abs(link.load - it->second);

            total += difference;
            maximum = std::max(maximum, difference);

            const double tolerance = 1e-9 * std::max({1.0, link.load, it->second});

            if (difference > tolerance) {
                ++changed;
            }

            loads.erase(it);
        }

        if (!loads.empty()) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Duplicate modified physical link.");
        }

        metrics.total_absolute_load_redistribution = total;
        metrics.maximum_absolute_link_load_change = maximum;
        metrics.changed_physical_links = changed;
    }

    return metrics;
}

#endif // E_ROUTING_NETWORK_DATA_H
