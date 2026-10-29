//
// Created by Mert on 29.10.26.
//

#ifndef E_ROUTING_VISUALIZATION_RESULT_H
#define E_ROUTING_VISUALIZATION_RESULT_H

#include <cstddef>
#include <string>
#include <vector>


struct LinkFailureRecoveryResult {
    int failed_edge_id = -1;

    int source = -1;
    int target = -1;

    bool graph_disconnected = false;

    double total_demand = 0.0;
    double unroutable_demand = 0.0;
    double unroutable_demand_fraction = 0.0;

    double baseline_congestion = -1.0;
    double post_failure_congestion = -1.0;
    double congestion_increase_factor = -1.0;

    double recomputation_runtime_microseconds = -1.0;

    bool recomputation_attempted = false;
    bool recomputation_succeeded = false;
};

struct LinkFailureRecoverySummary {
    std::size_t tested_links = 0;
    std::size_t disconnected_failures = 0;
    std::size_t successful_recomputations = 0;
    std::size_t failed_recomputations = 0;

    double maximum_unroutable_demand_fraction = 0.0;
    double average_unroutable_demand_fraction = 0.0;

    double maximum_post_failure_congestion = -1.0;
    double average_post_failure_congestion = -1.0;

    double maximum_congestion_increase_factor = -1.0;
    double average_congestion_increase_factor = -1.0;

    double average_recomputation_runtime_microseconds = -1.0;
    double maximum_recomputation_runtime_microseconds = -1.0;

    int worst_failed_edge_id = -1;
    int worst_failed_source = -1;
    int worst_failed_target = -1;

    int worst_disconnect_edge_id = -1;
    int worst_disconnect_source = -1;
    int worst_disconnect_target = -1;

    double worst_disconnect_unroutable_demand_fraction = -1.0;

    int worst_congestion_edge_id = -1;
    int worst_congestion_source = -1;
    int worst_congestion_target = -1;

    double worst_congestion_baseline = -1.0;
    double worst_congestion_post_failure = -1.0;
    double worst_congestion_increase_factor = -1.0;

    int slowest_recovery_edge_id = -1;
    int slowest_recovery_source = -1;
    int slowest_recovery_target = -1;

    double slowest_recovery_runtime_microseconds = -1.0;
};

struct FailureRecoveryAnalysis {
    LinkFailureRecoverySummary summary;
    std::vector<LinkFailureRecoveryResult> failures;
};


struct AffectedCommodity {
    int source = -1;
    int target = -1;

    double demand = 0.0;

    /*
     * Fraction of the unit routing for this commodity that crosses the
     * failed undirected link.
     */
    double failed_fraction = 0.0;

    /*
     * demand × failed_fraction
     */
    double lost_traffic = 0.0;
};

struct LinkFailureSummary {
    int failed_edge_id = -1;
    int anti_edge_id = -1;

    int source = -1;
    int target = -1;

    double capacity = 0.0;

    std::size_t affected_commodities = 0;

    double affected_demand_volume = 0.0;
    double lost_traffic = 0.0;
    double total_demand = 0.0;

    double affected_demand_fraction = 0.0;
    double lost_traffic_fraction = 0.0;

    std::vector<AffectedCommodity> commodities;
};
struct VisualizationNode {
    int id = -1;
    std::string label;
};

struct VisualizationEdge {
    int id = -1;
    int source = -1;
    int target = -1;

    double capacity = 0.0;
    double load = 0.0;
    double utilization = 0.0;

    bool enabled = true;
};


struct RoutingAnalysisSummary {
    double maximum_congestion = 0.0;
    double average_utilization = 0.0;
    double total_routed_demand = 0.0;

    std::size_t overloaded_edges = 0;
};


struct LinkFailureAnalysisSummary {
    std::size_t tested_links = 0;

    /*
     * Link whose failure affects the largest fraction
     * of total demand.
     */
    int most_critical_edge_id = -1;
    int most_critical_source = -1;
    int most_critical_target = -1;

    /*
     * Maximum lost fraction across all single-link
     * failure scenarios.
     */
    double maximum_lost_traffic_fraction = 0.0;

    /*
     * Average over all tested physical links.
     */
    double average_lost_traffic_fraction = 0.0;

    /*
     * Median over all tested physical links.
     */
    double median_lost_traffic_fraction = 0.0;

    /*
     * Maximum fraction of total demand volume whose
     * commodities touch a single failed link.
     */
    double maximum_affected_demand_fraction = 0.0;

    /*
     * Number of links whose failure removes traffic
     * from at least one commodity.
     */
    std::size_t traffic_carrying_links = 0;

    /*
     * Number of links whose failure would lose at
     * least 10%, 25%, or 50% of total demand.
     */
    std::size_t critical_links_10_percent = 0;
    std::size_t critical_links_25_percent = 0;
    std::size_t critical_links_50_percent = 0;
};

struct RoutingVisualizationResult {
    std::string graph_name;
    std::string solver_name;
    std::string demand_model;

    std::vector<VisualizationNode> nodes;
    std::vector<VisualizationEdge> edges;

    RoutingAnalysisSummary summary;
    LinkFailureAnalysisSummary failure_summary;

    /*
     * Runtime for normal demand evaluation:
     * route demand + compute utilization/congestion.
     *
     * Excludes failure analysis.
     */
    double demand_evaluation_runtime_microseconds = -1.0;

    /*
     * Runtime for static N-1 link-failure analysis only.
     */
    double failure_analysis_runtime_microseconds = -1.0;

    std::vector<LinkFailureSummary> link_failures;
};



#endif // E_ROUTING_VISUALIZATION_RESULT_H
