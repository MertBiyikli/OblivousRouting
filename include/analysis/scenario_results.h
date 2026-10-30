//
// Created by Mert on 26.10.26.
//

#ifndef E_ROUTING_SCENARIO_RESULTS_H
#define E_ROUTING_SCENARIO_RESULTS_H

#include "core/types.h"
#include "routing/routing_result.h"
#include "network_data.h"
#include "routing/routing_table.h"
#include <optional>
#include <vector>
#include "analysis/network_data.h"
#include "routing/routing_result.h"
#include "core/types.h"

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

// Keep this shared between the analysis entry point and result model.
enum class ExecutionMode { FAST, STANDARD, RESEARCH };
enum class ScenarioKind { GENERIC, CAPACITY_UPGRADE, LINK_FAILURE, TRAFFIC_GROWTH };


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



struct SolverAnalysisResult {
    std::vector<std::string> solver_algorithms;

    // Graphs must outlive the associated routing schemes.
    std::shared_ptr<optimized::Graph<EdgeData>> baseline_graph;
    std::shared_ptr<optimized::Graph<EdgeData>> modified_graph;

    std::optional<IRoutingResult> baseline_routing;
    std::optional<IRoutingResult> modified_routing;

    double baseline_solver_seconds = 0.0;
    double modified_solver_seconds = 0.0;

    double baseline_evaluation_seconds = 0.0;
    double modified_evaluation_seconds = 0.0;

    std::optional<NetworkData> baseline_Data;
    std::optional<NetworkData> modified_Data;
    std::optional<NetworkData> fixed_routing_Data;

    std::optional<ComparisonData> comparison;
    std::optional<ComparisonData> fixed_routing_comparison;

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
