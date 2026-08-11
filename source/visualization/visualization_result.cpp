//
// Created by Mert Biyikli on 24.07.26.
//

#include "visualization/visualization_result.h"

Result<RoutingVisualizationResult> RoutingAnalyzer::analyze(const optimized::Graph<EdgeData> &graph, const RoutingScheme &scheme, const demands &demand_map, std::string graph_name, std::string solver_name, std::string demand_model) {
    RoutingVisualizationResult result;

    result.graph_name = std::move(graph_name);
    result.solver_name = std::move(solver_name);
    result.demand_model = std::move(demand_model);

    /*
     * Export all graph vertices.
     */
    result.nodes.reserve(
        static_cast<std::size_t>(graph.getNumNodes())
    );

    for (int vertex = 0; vertex < graph.getNumNodes(); ++vertex) {
        result.nodes.push_back({
            .id = vertex,
            .label = std::to_string(vertex)
        });
    }

    /*
     * routeDemands() fills the vector with per-edge utilization.
     *
     * In the current routing implementation, opposite directed edges are
     * accumulated into the canonical representation of the undirected edge.
     */
    std::vector<double> utilization_by_edge(
        static_cast<std::size_t>(graph.getNumDirectedEdges()),
        0.0
    );

    scheme.routeDemands(
        utilization_by_edge,
        demand_map
    );

    result.edges.reserve(
        static_cast<std::size_t>(graph.getNumUndirectedEdges())
    );

    double utilization_sum = 0.0;

    /*
     * The graph stores two directed representations for each undirected edge.
     * We export only one canonical representation.
     *
     * This assumes the graph contains both (u,v) and (v,u), and vertex IDs
     * provide a stable ordering.
     */
    for (int edge = 0; edge < graph.getNumDirectedEdges(); ++edge) {
        const auto& [source, target] = graph.getEdgeEndpoints(edge);

        if (source > target) {
            continue;
        }

        const double capacity = graph.edgeData(edge).capacity;

        if (capacity <= 0.0) {
            return makeErrorMessage(ErrorCode::RuntimeError, "Edge capacity must be positive");
        }

        const double utilization =utilization_by_edge.at(static_cast<std::size_t>(edge));

        const double load = utilization * capacity;

        result.edges.push_back({
            .id = edge,
            .source = source,
            .target = target,
            .capacity = capacity,
            .load = load,
            .utilization = utilization,
            .enabled = true
        });

        utilization_sum += utilization;

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
        result.summary.average_utilization = utilization_sum /static_cast<double>(result.edges.size());
    }

    /*
     * Sum the original demand values.
     */
    for (std::size_t index = 0;index < demand_map.size();++index) {
        result.summary.total_routed_demand +=
            demand_map.getDemandValue(index);
    }

    auto failure_analysis =
    analyzeSingleLinkFailures(
        graph,
        scheme,
        demand_map,
        result
    );

    if (!failure_analysis) {
        return getError(failure_analysis);
    }

    return result;
}

Result<void> RoutingAnalyzer::analyzeSingleLinkFailures(const optimized::Graph<EdgeData>& graph,const RoutingScheme& scheme,const demands& demand_map,RoutingVisualizationResult& result) {
    result.link_failures.reserve(
        static_cast<std::size_t>(
            graph.getNumUndirectedEdges()
        )
    );

    for (
        int edge = 0;
        edge < graph.getNumDirectedEdges();
        ++edge
    ) {
        const auto& rev_edge = graph.reverse(graph.getEdge(graph.getEdgeEndpoints(edge).first, 0));
        int anti = rev_edge.id;

        if (anti == INVALID_EDGE_ID) {
            return makeErrorMessage(
                ErrorCode::InvalidGraph,
                "RoutingAnalyzer: edge has no anti-edge."
            );
        }

        /*
         * Process each physical link once.
         */
        if (edge > anti) {
            continue;
        }

        auto failure =
            LinkFailureAnalyzer::analyze(
                graph,
                scheme,
                demand_map,
                edge
            );

        if (!failure) {
            return getError(failure);
        }

        result.link_failures.push_back(
            std::move(failure.value())
        );
    }

    return {};
}