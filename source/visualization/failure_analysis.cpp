//
// Created by Mert Biyikli on 24.07.26.
//

#include "visualization/failure_analysis.h"
#include "routing/storage/linear_routing_table.h"
#include "utils/my_math.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

constexpr double flow_epsilon = 1e-12;

/*
 * RoutingScheme::getFlow() may encode the orientation of an undirected edge
 * through its sign.
 *
 * For failure analysis, either orientation of the physical link is lost.
 * Therefore, we sum the absolute routed fractions represented by the failed
 * edge and its anti-edge carefully.
 */
double failedLinkFraction(const optimized::Graph<EdgeData>& graph, const RoutingScheme& scheme, const int edge, const int anti_edge, const int source, const int target) {
    const double forward = scheme.getFlow(edge, source, target);

    if (!std::isfinite(forward)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    /*
     * In several of your routing schemes, the canonical edge stores signed
     * undirected flow and the anti-edge is not independently populated.
     *
     * If that is the universal convention in RoutingScheme, abs(forward)
     * alone is sufficient.
     *
     * The anti-edge fallback supports schemes that explicitly store both
     * directed orientations.
     */
    if (std::abs(forward) > flow_epsilon) {
        return std::abs(forward);
    }

    if (anti_edge == INVALID_EDGE_ID) {
        return 0.0;
    }

    const double reverse = scheme.getFlow(anti_edge, source, target);

    if (!std::isfinite(reverse)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return std::abs(reverse);
}

Result<LinkFailureSummary> LinkFailureAnalyzer::analyze(const optimized::Graph<EdgeData>& graph, const RoutingScheme& scheme, const demands& demand_map, const int failed_edge, double total_demand) {
    if (failed_edge < 0 || failed_edge >= graph.getNumDirectedEdges()) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "LinkFailureAnalyzer: invalid failed edge ID.");
    }

    const int anti_edge = graph.reverse(failed_edge).id;

    if (anti_edge == INVALID_EDGE_ID) {
        return makeErrorMessage(ErrorCode::InvalidGraph, "LinkFailureAnalyzer: failed edge has no anti-edge.");
    }

    const auto [edge_source, edge_target] = graph.getEdgeEndpoints(failed_edge);

    LinkFailureSummary result;

    result.failed_edge_id = failed_edge;
    result.anti_edge_id = anti_edge;
    result.source = edge_source;
    result.target = edge_target;
    result.capacity = graph.edgeData(failed_edge).capacity;
    result.total_demand = total_demand;

    /*
     * Only aggregate failure statistics are retained.
     *
     * Previously, every affected commodity was appended to
     * result.commodities for every physical edge. On medium-sized graphs,
     * this produced visualization JSON files hundreds of megabytes large.
     *
     * The browser currently uses only:
     *
     * - affected_commodities
     * - affected_demand_volume
     * - lost_traffic
     * - total_demand
     * - affected_demand_fraction
     * - lost_traffic_fraction
     */
    const auto* linear = dynamic_cast<const LinearRoutingScheme*>(&scheme);
    std::vector<double> basis;

    if (linear != nullptr) {
        linear->getEdgeBasis(failed_edge, basis);
    }

    for (std::size_t index = 0; index < demand_map.size(); ++index) {
        const auto& [source, target] = demand_map.getDemandPair(index);
        const double demand = demand_map.getDemandValue(index);

        if (!std::isfinite(demand) || demand < 0.0) {
            return makeErrorMessage(ErrorCode::InvalidDemand, "LinkFailureAnalyzer: demand contains an invalid value.");
        }

        if (demand <= flow_epsilon) {
            continue;
        }

        double failed_fraction = 0.0;

        if (linear != nullptr) {
            const double signed_flow = basis[source] - basis[target];

            // Reproduce LinearRoutingScheme::getFlow() semantics.
            const double routed_flow = std::abs(signed_flow) < EPS ? 0.0 : signed_flow;
            failed_fraction = std::abs(routed_flow);
        } else {
            failed_fraction = failedLinkFraction(graph, scheme, failed_edge, anti_edge, source, target);
        }

        if (!std::isfinite(failed_fraction)) {
            return makeErrorMessage(ErrorCode::InvalidRouting, "LinkFailureAnalyzer: routing scheme returned non-finite flow.");
        }

        if (failed_fraction <= flow_epsilon) {
            continue;
        }

        const double clamped_fraction = std::clamp(failed_fraction, 0.0, 1.0);
        const double lost_traffic = demand * clamped_fraction;

        ++result.affected_commodities;
        result.affected_demand_volume += demand;
        result.lost_traffic += lost_traffic;
    }

    if (result.total_demand > flow_epsilon) {
        result.affected_demand_fraction = result.affected_demand_volume / result.total_demand;

        result.lost_traffic_fraction = result.lost_traffic / result.total_demand;
    }

    return result;
}