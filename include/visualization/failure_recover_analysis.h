//
// Created by Mert on 07.10.26.
//

#ifndef E_ROUTING_FAILURE_RECOVER_ANALYSIS_H
#define E_ROUTING_FAILURE_RECOVER_ANALYSIS_H

#include <cstddef>
#include <vector>

#include "core/errors.h"
#include "core/types.h"
#include "data_structures/graph/graph.h"
#include "routing/routing_table.h"
#include "utils/demands.h"
#include "algorithms/semi_oblivious/semi_oblivious_solver.h"
#include "algorithms/semi_oblivious/semi_routing_engine.h"
#include "algorithms/semi_oblivious/postprocessing/or_tools_optimizer.h"


/*
 * Layer-2 N-1 recovery result for one physical link failure.
 *
 * Layer 1 asks:
 *
 *     "How much of the current routing depends on this link?"
 *
 * Layer 2 asks:
 *
 *     "What happens if the link is actually removed and the
 *      routing algorithm is allowed to recompute?"
 */
struct LinkFailureRecoveryResult {
    int failed_edge_id = -1;

    int source = -1;
    int target = -1;

    /*
     * True if removing the physical link disconnects the topology.
     */
    bool graph_disconnected = false;

    /*
     * Demand that cannot communicate after the topology failure.
     */
    double total_demand = 0.0;

    double unroutable_demand = 0.0;

    double unroutable_demand_fraction = 0.0;

    /*
     * Congestion before and after failure.
     */
    double baseline_congestion = -1.0;

    double post_failure_congestion = -1.0;

    double congestion_increase_factor = -1.0;

    /*
     * Fresh solver construction + solve time on the failed topology.
     *
     * Demand evaluation is deliberately not included here.
     */
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

    /*
 * ------------------------------------------------------------
 * Overall worst failure
 * ------------------------------------------------------------
 *
 * Preserved for backwards compatibility.
 *
 * The ordering is primarily based on unroutable demand and
 * secondarily on congestion degradation.
 */
    int worst_failed_edge_id = -1;
    int worst_failed_source = -1;
    int worst_failed_target = -1;


    /*
     * ------------------------------------------------------------
     * Worst disconnecting failure
     * ------------------------------------------------------------
     *
     * Among failures that make demand unreachable, this is the
     * physical link whose removal produces the largest fraction
     * of unroutable demand.
     */
    int worst_disconnect_edge_id = -1;
    int worst_disconnect_source = -1;
    int worst_disconnect_target = -1;

    double worst_disconnect_unroutable_demand_fraction = -1.0;


    /*
     * ------------------------------------------------------------
     * Worst survivable congestion failure
     * ------------------------------------------------------------
     *
     * The topology remains connected and the solver successfully
     * recomputes a routing scheme, but this failure produces the
     * largest absolute post-failure congestion.
     */
    int worst_congestion_edge_id = -1;
    int worst_congestion_source = -1;
    int worst_congestion_target = -1;

    double worst_congestion_baseline = -1.0;
    double worst_congestion_post_failure = -1.0;
    double worst_congestion_increase_factor = -1.0;


    /*
     * ------------------------------------------------------------
     * Slowest successful recovery
     * ------------------------------------------------------------
     */
    int slowest_recovery_edge_id = -1;
    int slowest_recovery_source = -1;
    int slowest_recovery_target = -1;

    double slowest_recovery_runtime_microseconds = -1.0;


};


struct FailureRecoveryAnalysis {
    LinkFailureRecoverySummary summary;

    std::vector<LinkFailureRecoveryResult> failures;
};


class FailureRecoveryAnalyzer {
public:

    /*
     * Performs complete N-1 physical-link recovery analysis.
     *
     * Currently intended for normal oblivious solvers:
     *
     *   ELECTRICAL_SKETCHING
     *   RAECKE_CKR_FLAT
     *   RAECKE_FRT_FLAT
     *   EXPANDER_HIERARCHY
     *   ...
     *
     * Semi-oblivious solvers will be added after this path is
     * validated because their recovery lifecycle also includes
     * candidate-path preprocessing + demand-specific optimization.
     */
    static Result<FailureRecoveryAnalysis> analyze(optimized::Graph<EdgeData>& graph, SolverType solver_type, const demands& demand_map, const DemandModelType& demand_model, double baseline_congestion);

private:
    static Result<std::unique_ptr<RoutingScheme>> recomputeSemiOblivious(
        optimized::Graph<EdgeData>& graph,
        SolverType solver_type,
        const demands& demand_map,
        DemandModelType demand_model
    );

    static bool isSemiObliviousSolver(SolverType solver_type);
};

#endif //E_ROUTING_FAILURE_RECOVER_ANALYSIS_H
