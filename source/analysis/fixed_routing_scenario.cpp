//
// Created by Mert on 22.10.26.
//

#include "analysis/fixed_routing_scenario.h"

Result<RoutingVisualizationResult> transformScenario(const RoutingVisualizationResult& baseline, double demandScale, const std::unordered_map<int, double>& capacityFactors) {
    if (demandScale <= 0 || !std::isfinite(demandScale)) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Demand scale must be finite and non-negative.");
    }

    if (baseline.edges.empty()) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Baseline edges must not be empty.");
    }

    for (const auto& [edgeId, factor] : capacityFactors) {
        if (edgeId < 0 || factor < 1 || !std::isfinite(factor)) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid edge upgrade for edge " + std::to_string(edgeId) + " with factor " + std::to_string(factor) + ".");
        }
    }

    RoutingVisualizationResult result = baseline;

    // A transformed result has no newly computed failure evidence.
    result.link_failures.clear();
    result.failure_summary = {};
    result.failure_analysis_runtime_microseconds = 0.0;
    result.demand_evaluation_runtime_microseconds = 0.0;

    result.summary.maximum_congestion = 0.0;
    result.summary.average_utilization = 0.0;
    result.summary.overloaded_edges = 0;

    double utilizationSum = 0.0;
    std::unordered_set<int> observedEdges;

    for (auto& edge : result.edges) {
        if (!observedEdges.insert(edge.id).second) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Duplicate baseline edge ID");
        }

        if (!std::isfinite(edge.load) || !std::isfinite(edge.capacity) || edge.load < 0.0 || edge.capacity <= 0.0) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid baseline link metrics");
        }

        const auto it = capacityFactors.find(edge.id);

        const double capacityFactor = it == capacityFactors.end() ? 1.0 : it->second;

        const double newLoad = edge.load * demandScale;
        const double newCapacity = edge.capacity * capacityFactor;
        const double newUtilization = newLoad / newCapacity;

        if (!std::isfinite(newLoad) || !std::isfinite(newCapacity) || !std::isfinite(newUtilization) || newCapacity <= 0.0) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Scenario transformation overflow");
        }

        edge.load = newLoad;
        edge.capacity = newCapacity;
        edge.utilization = newUtilization;

        result.summary.maximum_congestion = std::max(result.summary.maximum_congestion, newUtilization);

        utilizationSum += newUtilization;

        if (newUtilization > 1.0) {
            ++result.summary.overloaded_edges;
        }
    }

    for (const auto& [edgeId, factor] : capacityFactors) {
        if (!observedEdges.contains(edgeId)) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Capacity upgrade references unknown edge");
        }
    }

    result.summary.average_utilization = utilizationSum / static_cast<double>(result.edges.size());

    result.summary.total_routed_demand = baseline.summary.total_routed_demand * demandScale;

    if (!std::isfinite(result.summary.average_utilization) || !std::isfinite(result.summary.total_routed_demand)) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Scenario aggregate overflow");
    }

    return result;
}