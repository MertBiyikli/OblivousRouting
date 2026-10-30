//
// Created by Mert on 29.10.26.
//

#ifndef E_ROUTING_LINK_FAILURE_RECOVERY_RESULT_H
#define E_ROUTING_LINK_FAILURE_RECOVERY_RESULT_H

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

#endif // E_ROUTING_LINK_FAILURE_RECOVERY_RESULT_H
