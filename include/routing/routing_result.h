//
// Created by Mert Biyikli on 22.06.26.
//

#ifndef OBLIVIOUSROUTING_ROUTING_RESULT_H
#define OBLIVIOUSROUTING_ROUTING_RESULT_H

#include <unordered_map>
#include <vector>
#include <memory>
#include "io/demand_io.h"
#include "routing_table.h"
#include "core/types.h"
#include "utils/my_math.h"
#include "visualization/failure_analysis.h"


struct DemandEvaluationResult {
    DemandModelType demand_type{};

    double congestion = -1.0;


    /*
     * Demand evaluation only:
     *
     * routing the demand through the already constructed
     * routing scheme and computing congestion.
     *
     * Does NOT include failure analysis.
     */
    double runtime_microseconds = -1.0;

    /*
     * Static N-1 link-failure exposure analysis only.
     */
    double failure_analysis_runtime_microseconds = -1.0;

    /*
     * Static single-link failure exposure metrics.
     */
    std::size_t failure_tested_links = 0;

    int failure_most_critical_edge_id = -1;
    int failure_most_critical_source = -1;
    int failure_most_critical_target = -1;

    double failure_maximum_lost_traffic_fraction = 0.0;
    double failure_average_lost_traffic_fraction = 0.0;
    double failure_median_lost_traffic_fraction = 0.0;

    double failure_maximum_affected_demand_fraction = 0.0;

    std::size_t failure_traffic_carrying_links = 0;

    std::size_t failure_critical_links_10_percent = 0;
    std::size_t failure_critical_links_25_percent = 0;
    std::size_t failure_critical_links_50_percent = 0;

    /*
 * ------------------------------------------------------------
 * Layer-2 failure recovery metrics
 * ------------------------------------------------------------
 */

    bool failure_recovery_available = false;

    std::size_t recovery_tested_links = 0;

    std::size_t recovery_disconnected_failures = 0;

    std::size_t recovery_successful_recomputations = 0;

    std::size_t recovery_failed_recomputations = 0;

    double recovery_maximum_unroutable_demand_fraction = -1.0;

    double recovery_average_unroutable_demand_fraction = -1.0;

    double recovery_maximum_post_failure_congestion = -1.0;

    double recovery_average_post_failure_congestion = -1.0;

    double recovery_maximum_congestion_increase_factor = -1.0;

    double recovery_average_congestion_increase_factor = -1.0;

    double recovery_average_recomputation_runtime_microseconds = -1.0;

    double recovery_maximum_recomputation_runtime_microseconds = -1.0;

    int recovery_worst_failed_edge_id = -1;
    int recovery_worst_failed_source = -1;
    int recovery_worst_failed_target = -1;

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

    /*
 * Worst disconnecting failure.
 */
    int recovery_worst_disconnect_edge_id = -1;
    int recovery_worst_disconnect_source = -1;
    int recovery_worst_disconnect_target = -1;

    double recovery_worst_disconnect_unroutable_demand_fraction = -1.0;


    /*
     * Worst survivable congestion failure.
     */
    int recovery_worst_congestion_edge_id = -1;
    int recovery_worst_congestion_source = -1;
    int recovery_worst_congestion_target = -1;

    double recovery_worst_congestion_baseline = -1.0;
    double recovery_worst_congestion_post_failure = -1.0;
    double recovery_worst_congestion_increase_factor = -1.0;


    /*
     * Slowest successful recovery.
     */
    int recovery_slowest_recovery_edge_id = -1;
    int recovery_slowest_recovery_source = -1;
    int recovery_slowest_recovery_target = -1;

    double recovery_slowest_recovery_runtime_microseconds = -1.0;
};

class MWUMetrics{
public:
    MWUMetrics() {
        iteration_count = 0;
        solve_time = 0;
        transformation_time = 0;
        mwu_weight_update_time = 0;
    }

    virtual ~MWUMetrics() = default;


    std::vector<double> oracle_running_times;
    int iteration_count;
    double solve_time;
    double transformation_time;
    double mwu_weight_update_time;
    double load_computation_time{};

    [[nodiscard]] bool empty() const {
        return iteration_count == 0;
    }

    [[nodiscard]] double averageOracleTime() const {
        double sum = 0.0;
        for (double t : oracle_running_times) {
            sum += t;
        }
        return ((sum > EPS) ? sum / static_cast<double>(oracle_running_times.size()) : -1.0);
    }

    [[nodiscard]] int getIterationCount() const {
        return iteration_count;
    }

};

struct ExpanderMetrics {
    double hierarchy_runtime_microseconds = 0.0;
    double tree_runtime_microseconds = 0.0;
    double basis_flow_runtime_microseconds = 0.0;

    // Hierarchy structure.
    std::size_t hierarchy_levels = 0;
    std::size_t hierarchy_clusters = 0;
    std::vector<std::size_t> clusters_per_level;

    std::size_t max_cluster_vertices = 0;
    double average_cluster_vertices = 0.0;

    // Tree sparsifier structure.
    std::size_t tree_nodes = 0;
    std::size_t tree_edges = 0;
    int tree_depth = 0;

    // Routing construction.
    std::size_t basis_flows = 0;
    std::size_t total_electrical_solves = 0;

    // Numerical diagnostics over all root-target basis flows.
    double max_basis_embedding_congestion = 0.0;
    double max_conservation_error = 0.0;

    double preprocessingRuntime() const noexcept {
        return hierarchy_runtime_microseconds +
               tree_runtime_microseconds;
    }

    double averageElectricalSolves() const noexcept {
        if (basis_flows == 0) {
            return 0.0;
        }

        return static_cast<double>(total_electrical_solves) /
               static_cast<double>(basis_flows);
    }

    bool empty() const noexcept {
        return hierarchy_levels == 0 &&
               hierarchy_clusters == 0 &&
               tree_nodes == 0 &&
               basis_flows == 0;
    }
};



struct IRoutingResult {

    ResultStatus status = ResultStatus::OK; // By default
    std::string solver_name;
    std::string routing_base;
    std::string graph_name;
    int nodes, edges;
    SolverType type;
    double oblivious_ratio = -1;


    // Runtime
    double total_runtime_microseconds = -1;
    double preprocessing_runtime_microseconds = -1;
    double solve_runtime_microseconds = -1;

    // Optional: only filled for normal solvers or last scheme if needed
    std::unique_ptr<RoutingScheme> scheme;
    int candidate_paths = -1;
    double average_paths_per_pair = -1.0;

    MWUMetrics mwu_metrics;
    ExpanderMetrics expander_metrics;

    std::vector<DemandEvaluationResult> demand_evaluations;
    // Semi-oblivious: one demand-specific scheme per demand model
    std::unordered_map<std::string, std::unique_ptr<RoutingScheme>> demand_schemes;

    /*
     * One visualization result per evaluated demand model.
     */
    std::vector<RoutingVisualizationResult> visualization_results;

};

/*
 * The result of a routing experiment, which may include multiple solvers.
 */
struct RoutingExperimentResult{
    std::string graph_name;
    std::string graph_path;
    int nodes = 0, edges = 0;

    std::vector<IRoutingResult> solver_results;

};

struct SemiObliviousRoutingResult :public IRoutingResult{

    DemandModelType demand_type{};
    double congestion = -1.0;
    std::string path_selection_strategy;

    /*
    * Visualization for the exact demand-specific routing solution.
    */
    RoutingVisualizationResult visualization;
};




inline void printTimeStats(const MWUMetrics& metrics) {
    std::cout << "Solve time: " << metrics.solve_time << " micro seconds\n";
    std::cout << "Transformation time: " << metrics.transformation_time << " micro seconds\n";
    std::cout << "MWU iterations: " << metrics.iteration_count << "\n";
    std::cout << "MWU load computation: " << metrics.load_computation_time << " micro seconds\n";
    std::cout << "Average oracle time: " << metrics.averageOracleTime() << " micro seconds\n";
    std::cout << "Total MWU weight update time: " << metrics.mwu_weight_update_time << " micro seconds\n";
}


inline void printSemiObliviousResult(
    const SemiObliviousRoutingResult& r
) {
    std::cout << "Routing base: " << r.path_selection_strategy << std::endl;
    std::cout << "Demand [" << demandModelName(r.demand_type) << "]\n";
    std::cout << "  Congestion: " << r.congestion << '\n';
    std::cout << "  Runtime: " << r.total_runtime_microseconds << " us\n";
    std::cout << "  Candidate paths: " << r.candidate_paths << '\n';
    std::cout << "  Avg paths/pair: " << r.average_paths_per_pair << '\n';
}



#endif //OBLIVIOUSROUTING_ROUTING_RESULT_H