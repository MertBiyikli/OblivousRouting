//
// Created by Mert on 08.10.26.
//

#ifndef E_ROUTING_FAILURE_RECOVERY_TYPES_H
#define E_ROUTING_FAILURE_RECOVERY_TYPES_H


#include <cstddef>
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

#endif // E_ROUTING_FAILURE_RECOVERY_TYPES_H
