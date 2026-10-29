//
// Created by Mert on 26.10.26.
//

#ifndef E_ROUTING_SCENARIO_RESULTS_H
#define E_ROUTING_SCENARIO_RESULTS_H

#include "core/types.h"
#include "routing/routing_result.h"
#include "routing/routing_table.h"
#include <optional>
#include <vector>

// Keep this shared between the analysis entry point and result model.
enum class ExecutionMode { FAST, STANDARD, RESEARCH };
enum class ScenarioKind { GENERIC, CAPACITY_UPGRADE, LINK_FAILURE, TRAFFIC_GROWTH };

struct LinkData {
    int edge_id = -1; // canonical physical edge ID
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

struct ComparisonData {
    double congestion_before = 0.0;
    double congestion_after = 0.0;
    double absolute_reduction = 0.0;                         // positive indicates improvement
    std::optional<double> relative_reduction = std::nullopt; // undefined for zero baseline
};

// Call only after both Data have been evaluated against the same demands.
inline ComparisonData compareNetworkData(const NetworkData& before, const NetworkData& after) {
    ComparisonData comparison;
    comparison.congestion_before = before.max_congestion;
    comparison.congestion_after = after.max_congestion;
    comparison.absolute_reduction = before.max_congestion - after.max_congestion;
    if (before.max_congestion > 0.0)
        comparison.relative_reduction = comparison.absolute_reduction / before.max_congestion;
    return comparison;
}

struct CapacityUpgradeDetail {
    int edge_id = -1;
    double original_capacity = 0.0;
    double upgraded_capacity = 0.0;
    double capacity_factor = 1.0;
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
            if (difference > tolerance) ++changed;
            loads.erase(it); // Detect duplicates on modified side.
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
struct SolverAnalysisResult {
    std::vector<std::string> solver_algorithms;

    // Graphs must outlive their routing schemes.
    std::shared_ptr<optimized::Graph<EdgeData>> baseline_graph;
    std::shared_ptr<optimized::Graph<EdgeData>> modified_graph;

    std::optional<IRoutingResult> baseline_routing;
    std::optional<IRoutingResult> modified_routing;

    // Solver execution times.
    double baseline_solver_seconds = 0.0;
    double modified_solver_seconds = 0.0;

    // Demand evaluation times.
    double baseline_evaluation_seconds = 0.0;
    double modified_evaluation_seconds = 0.0;

    // Original network with original routing.
    std::optional<NetworkData> baseline_Data;

    // Upgraded network with recomputed routing.
    std::optional<NetworkData> modified_Data;

    // Upgraded network with unchanged original routing.
    std::optional<NetworkData> fixed_routing_Data;

    // Baseline versus reoptimized routing.
    std::optional<ComparisonData> comparison;

    // Baseline versus fixed routing.
    std::optional<ComparisonData> fixed_routing_comparison;

    // Present only after all three states were evaluated.
    std::optional<ScenarioPerformanceMetrics> performance_metrics;
};

struct ScenarioResult {
    ScenarioResult() = default;

    ScenarioResult(const ScenarioResult&) = delete;
    ScenarioResult& operator=(const ScenarioResult&) = delete;

    ScenarioResult(ScenarioResult&&) noexcept = default;
    ScenarioResult& operator=(ScenarioResult&&) noexcept = default;

    ScenarioKind kind = ScenarioKind::GENERIC;
    ExecutionMode execution_mode = ExecutionMode::FAST;
    std::uint64_t original_fingerprint = 0;
    std::uint64_t modified_fingerprint = 0;
    std::vector<CapacityUpgradeDetail> capacity_upgrades;
    std::vector<int> link_failures;
    std::vector<SolverAnalysisResult> solver_results;
};


#endif // E_ROUTING_SCENARIO_RESULTS_H
