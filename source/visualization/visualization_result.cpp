//
// Created by Mert Biyikli on 24.07.26.
//

#include "visualization/visualization_result.h"

#include "routing/routing_result.h"

Result<RoutingVisualizationResult>
RoutingAnalyzer::analyze(const optimized::Graph<EdgeData>& graph,const RoutingScheme& scheme,const demands& demand_map,
    const IRoutingResult& _result,
    std::string demand_model
) {
    RoutingVisualizationResult result;
    result.graph_name = _result.graph_name;


    result.solver_name = _result.solver_name;

    result.demand_model = std::move(demand_model);

    /*
     * ================================================================
     * Normal demand evaluation
     * ================================================================
     */

    const auto evaluation_start =
        std::chrono::steady_clock::now();

    result.nodes.reserve(
        static_cast<std::size_t>(
            graph.getNumNodes()
        )
    );

    for (
        int vertex = 0;
        vertex < graph.getNumNodes();
        ++vertex
    ) {
        result.nodes.push_back({
            .id = vertex,
            .label = std::to_string(vertex)
        });
    }

    std::vector<double> utilization_by_edge(
        static_cast<std::size_t>(
            graph.getNumDirectedEdges()
        ),
        0.0
    );

    scheme.routeDemands(
        utilization_by_edge,
        demand_map
    );

    result.edges.reserve(
        static_cast<std::size_t>(
            graph.getNumUndirectedEdges()
        )
    );

    double utilization_sum = 0.0;

    for (
        int edge = 0;
        edge < graph.getNumDirectedEdges();
        ++edge
    ) {
        const auto& [source, target] =
            graph.getEdgeEndpoints(edge);

        if (source > target) {
            continue;
        }

        const double capacity =
            graph.edgeData(edge).capacity;

        if (capacity <= 0.0) {
            return makeErrorMessage(
                ErrorCode::RuntimeError,
                "Edge capacity must be positive"
            );
        }

        const double utilization =
            utilization_by_edge.at(
                static_cast<std::size_t>(
                    edge
                )
            );

        const double load =
            utilization * capacity;

        result.edges.push_back({
            .id = edge,
            .source = source,
            .target = target,
            .capacity = capacity,
            .load = load,
            .utilization = utilization,
            .enabled = true
        });

        utilization_sum +=
            utilization;

        result.summary.maximum_congestion =
            std::max(
                result.summary.maximum_congestion,
                utilization
            );

        if (utilization > 1.0) {
            ++result.summary.overloaded_edges;
        }
    }

    if (!result.edges.empty()) {
        result.summary.average_utilization =
            utilization_sum /
            static_cast<double>(
                result.edges.size()
            );
    }

    for (
        std::size_t index = 0;
        index < demand_map.size();
        ++index
    ) {
        result.summary.total_routed_demand +=
            demand_map.getDemandValue(
                index
            );
    }

    const auto evaluation_end =
        std::chrono::steady_clock::now();

    result.demand_evaluation_runtime_microseconds =
        std::chrono::duration<
            double,
            std::micro
        >(
            evaluation_end
            - evaluation_start
        ).count();

    /*
     * ================================================================
     * Layer-1 static link-failure analysis
     * ================================================================
     */

    const auto failure_start =
        std::chrono::steady_clock::now();

    auto failure_analysis =
        analyzeSingleLinkFailures(
            graph,
            scheme,
            demand_map,
            result
        );

    const auto failure_end =
        std::chrono::steady_clock::now();

    result.failure_analysis_runtime_microseconds =
        std::chrono::duration<
            double,
            std::micro
        >(
            failure_end
            - failure_start
        ).count();

    if (!failure_analysis) {
        return getError(
            failure_analysis
        );
    }

    return result;
}
Result<void>
RoutingAnalyzer::analyzeSingleLinkFailures(
    const optimized::Graph<EdgeData>& graph,
    const RoutingScheme& scheme,
    const demands& demand_map,
    RoutingVisualizationResult& result
) {
    result.link_failures.clear();

    result.link_failures.reserve(
        static_cast<std::size_t>(
            graph.getNumUndirectedEdges()
        )
    );

    for (int e = 0;e < graph.getNumDirectedEdges();++e) {
        const int anti = graph.reverse(e).id;

        if (anti == INVALID_EDGE_ID) {
            return makeErrorMessage(ErrorCode::InvalidGraph,
                "RoutingAnalyzer: edge has no anti-edge.");
        }

        /*
         * Process each physical undirected link exactly once.
         */
        if (e > anti) {
            continue;
        }

        auto failure =LinkFailureAnalyzer::analyze(graph,scheme,demand_map,e);

        if (!failure) {
            return getError(failure);
        }

        result.link_failures.push_back(std::move(failure.value()));
    }

        result.failure_summary.tested_links =
    result.link_failures.size();

    if (result.link_failures.empty()) {
        return {};
    }

    std::vector<double> lost_fractions;

    lost_fractions.reserve(result.link_failures.size());

    double lost_fraction_sum = 0.0;

    for (const auto& failure :result.link_failures)
    {
        const double lost_fraction = failure.lost_traffic_fraction;

        lost_fractions.push_back(lost_fraction);

        lost_fraction_sum += lost_fraction;

        if (lost_fraction > 0.0) {
            ++result.failure_summary.traffic_carrying_links;
        }

        if (lost_fraction >= 0.10) {
            ++result.failure_summary.critical_links_10_percent;
        }

        if (lost_fraction >= 0.25) {
            ++result.failure_summary.critical_links_25_percent;
        }

        if (lost_fraction >= 0.50) {
            ++result.failure_summary.critical_links_50_percent;
        }

        if (failure.affected_demand_fraction>result.failure_summary.maximum_affected_demand_fraction) {
            result.failure_summary
                .maximum_affected_demand_fraction =
                failure.affected_demand_fraction;
        }

        if (lost_fraction>result.failure_summary.maximum_lost_traffic_fraction) {
            result.failure_summary.maximum_lost_traffic_fraction = lost_fraction;
            result.failure_summary.most_critical_edge_id = failure.failed_edge_id;
            result.failure_summary.most_critical_source = failure.source;
            result.failure_summary.most_critical_target = failure.target;
        }
    }

    result.failure_summary.average_lost_traffic_fraction =lost_fraction_sum /static_cast<double>(result.link_failures.size());

    std::sort(lost_fractions.begin(),lost_fractions.end());

    const std::size_t count = lost_fractions.size();

    if (count % 2 == 1) {
        result.failure_summary.median_lost_traffic_fraction =lost_fractions[count / 2];
    } else {
        result.failure_summary.median_lost_traffic_fraction =0.5 *(lost_fractions[count / 2 - 1]+lost_fractions[count / 2]);
    }
    return {};
}