//
// Created by Mert on 07.10.26.
//
#include "../include/visualization/failure_recover_analysis.h"
#include "algorithms/oblivious/mwu/electrical_mwu.h"
#include "algorithms/oblivious/mwu/oracle/tree/fast_ckr/fast_ckr.h"
#include "algorithms/oblivious/mwu/tree_mwu.h"
#include "algorithms/semi_oblivious/expander_hierarchy/tree_sparsifier_solver.h"
#include "algorithms/semi_oblivious/postprocessing/or_tools_optimizer.h"
#include "algorithms/semi_oblivious/semi_oblivious_solver.h"
#include "algorithms/semi_oblivious/semi_routing_engine.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <span>
#include <vector>

#include "io/solver_io.h"
#include "utils/time_tracking.h"

constexpr double recovery_epsilon = 1e-12;

bool FailureRecoveryAnalyzer::isSemiObliviousSolver(const SolverType solver_type) {
    return solver_type == SolverType::SEMI_ELECTRICAL || solver_type == SolverType::SEMI_TREE || solver_type == SolverType::SEMI_EXPANDER_HIERARCHY;
}

bool reachable(const optimized::Graph<EdgeData>& graph, const int source, const int target) {
    if (source == target) {
        return true;
    }

    if (!graph.alive(source) || !graph.alive(target)) {
        return false;
    }

    std::vector<std::uint8_t> visited(static_cast<std::size_t>(graph.getNumNodes()), 0);

    /*
     * graph.size() is the number of active vertices.
     *
     * At Layer 2 we currently never remove vertices, so the
     * original IDs remain 0 .. n-1.
     */
    visited.resize(static_cast<std::size_t>(graph.getNumNodes()), 0);

    std::queue<int> queue;

    visited[static_cast<std::size_t>(source)] = 1;

    queue.push(source);

    while (!queue.empty()) {
        const int u = queue.front();

        queue.pop();

        for (const auto& edge : graph.edgesOf(u)) {
            const int v = edge.tail;

            if (visited[static_cast<std::size_t>(v)] != 0) {
                continue;
            }

            if (v == target) {
                return true;
            }

            visited[static_cast<std::size_t>(v)] = 1;

            queue.push(v);
        }
    }

    return false;
}

/*
 * Does the entire active topology remain connected?
 */
bool graphConnected(const optimized::Graph<EdgeData>& graph) {
    if (graph.getNumNodes() <= 1) {
        return true;
    }

    const int start = *graph.begin();

    std::vector<std::uint8_t> visited(static_cast<std::size_t>(graph.getNumNodes()), 0);

    std::queue<int> queue;

    visited[static_cast<std::size_t>(start)] = 1;

    queue.push(start);

    std::size_t reached = 1;

    while (!queue.empty()) {
        const int u = queue.front();

        queue.pop();

        for (const auto& edge : graph.edgesOf(u)) {
            const int v = edge.tail;

            if (visited[static_cast<std::size_t>(v)] != 0) {
                continue;
            }

            visited[static_cast<std::size_t>(v)] = 1;

            ++reached;

            queue.push(v);
        }
    }

    return reached == static_cast<std::size_t>(graph.getNumNodes());
}

/*
 * --------------------------------------------------------------------------
 * Demand reachability
 * --------------------------------------------------------------------------
 */

void calculateUnroutableDemand(const optimized::Graph<EdgeData>& graph, const demands& demand_map, LinkFailureRecoveryResult& result) {
    for (std::size_t index = 0; index < demand_map.size(); ++index) {
        const auto [source, target] = demand_map.getDemandPair(index);

        const double demand = demand_map.getDemandValue(index);

        if (!std::isfinite(demand) || demand < 0.0) {
            continue;
        }

        result.total_demand += demand;

        if (demand <= recovery_epsilon) {
            continue;
        }

        if (!reachable(graph, source, target)) {
            result.unroutable_demand += demand;
        }
    }

    if (result.total_demand > recovery_epsilon) {
        result.unroutable_demand_fraction = result.unroutable_demand / result.total_demand;
    }
}

/*
 * --------------------------------------------------------------------------
 * Congestion
 * --------------------------------------------------------------------------
 */

double evaluateCongestion(const optimized::Graph<EdgeData>& graph, const RoutingScheme& scheme, const demands& demand_map) {
    /*
     * The reindexed failure graph has contiguous directed edge IDs.
     */
    std::vector<double> utilization(static_cast<std::size_t>(graph.getNumDirectedEdges()), 0.0);

    scheme.routeDemands(utilization, demand_map);

    double maximum_congestion = 0.0;

    for (const double value : utilization) {
        if (!std::isfinite(value)) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        maximum_congestion = std::max(maximum_congestion, value);
    }

    return maximum_congestion;
}

/*
 * --------------------------------------------------------------------------
 * Failure graph
 * --------------------------------------------------------------------------
 *
 * Graph::removeEdge() preserves stable IDs. Most routing algorithms were
 * originally written assuming a compact graph.
 *
 * Therefore, after removing the failed physical link we create an
 * independent reindexed graph containing all active vertices.
 *
 * Because Layer 2 removes only links and no vertices, supplying vertices
 * in natural order preserves the original vertex IDs while reindexing
 * edges contiguously.
 */

Result<optimized::Graph<EdgeData>::ReindexedSubgraph> makeFailureGraph(const optimized::Graph<EdgeData>& graph) {
    std::vector<int> vertices;

    vertices.reserve(static_cast<std::size_t>(graph.getNumNodes()));

    for (const int vertex : graph) {
        vertices.push_back(vertex);
    }

    return graph.reindexedSubgraph(vertices.begin(), vertices.end());
}

/*
 * --------------------------------------------------------------------------
 * Worst-failure ordering
 * --------------------------------------------------------------------------
 *
 * Primary severity:
 *
 *     unroutable demand fraction
 *
 * Secondary severity:
 *
 *     congestion increase factor
 */

bool worseFailure(const LinkFailureRecoveryResult& candidate, const LinkFailureRecoveryResult& current) {
    if (candidate.unroutable_demand_fraction > current.unroutable_demand_fraction + recovery_epsilon) {
        return true;
    }

    if (candidate.unroutable_demand_fraction + recovery_epsilon < current.unroutable_demand_fraction) {
        return false;
    }

    if (candidate.congestion_increase_factor >= 0.0 && (current.congestion_increase_factor < 0.0 || candidate.congestion_increase_factor > current.congestion_increase_factor)) {
        return true;
    }

    return false;
}

Result<FailureRecoveryAnalysis> FailureRecoveryAnalyzer::analyze(optimized::Graph<EdgeData>& graph, const SolverType solver_type, const demands& demand_map, const DemandModelType& demand_model, double baseline_congestion) {
    FailureRecoveryAnalysis analysis;

    double unroutable_fraction_sum = 0.0;

    double post_failure_congestion_sum = 0.0;

    double congestion_increase_sum = 0.0;

    double recomputation_runtime_sum = 0.0;

    std::size_t post_failure_congestion_count = 0;

    std::size_t congestion_increase_count = 0;

    std::size_t recomputation_runtime_count = 0;

    bool have_worst_failure = false;

    LinkFailureRecoveryResult worst_failure;

    bool have_worst_disconnect = false;
    LinkFailureRecoveryResult worst_disconnect;

    bool have_worst_congestion = false;
    LinkFailureRecoveryResult worst_congestion;

    bool have_slowest_recovery = false;
    LinkFailureRecoveryResult slowest_recovery;

    /*
     * globalEdgeCount() is intentionally used here.
     *
     * Physical edge IDs remain stable even while edges are temporarily
     * removed.
     */
    const int physical_edges = graph.globalEdgeCount();

    analysis.failures.reserve(static_cast<std::size_t>(physical_edges));

    for (int undirected_edge = 0; undirected_edge < physical_edges; ++undirected_edge) {
        /*
         * Canonical directed representation of the physical edge.
         */
        const int edge_id = 2 * undirected_edge;

        if (!graph.edgeAlive(edge_id)) {
            continue;
        }

        const auto [source, target] = graph.getEdgeEndpoints(edge_id);

        LinkFailureRecoveryResult failure;

        failure.failed_edge_id = edge_id;

        failure.source = source;

        failure.target = target;

        failure.baseline_congestion = baseline_congestion;

        /*
         * --------------------------------------------------------------
         * Remove physical link
         * --------------------------------------------------------------
         */

        const int checkpoint = graph.checkpoint();

        if (!graph.removeEdge(edge_id)) {
            auto rollback = graph.rollback(checkpoint);

            if (!rollback) {
                return getError(rollback);
            }

            continue;
        }

        ++analysis.summary.tested_links;

        /*
         * --------------------------------------------------------------
         * Topological failure impact
         * --------------------------------------------------------------
         */

        failure.graph_disconnected = !graphConnected(graph);

        if (failure.graph_disconnected) {
            ++analysis.summary.disconnected_failures;
        }

        calculateUnroutableDemand(graph, demand_map, failure);

        /*
         * --------------------------------------------------------------
         * Worst disconnecting physical link
         * --------------------------------------------------------------
         */
        if (failure.graph_disconnected && (!have_worst_disconnect || failure.unroutable_demand_fraction > worst_disconnect.unroutable_demand_fraction)) {
            worst_disconnect = failure;

            have_worst_disconnect = true;
        }

        unroutable_fraction_sum += failure.unroutable_demand_fraction;

        analysis.summary.maximum_unroutable_demand_fraction = std::max(analysis.summary.maximum_unroutable_demand_fraction,

                                                                       failure.unroutable_demand_fraction);

        /*
         * --------------------------------------------------------------
         * Recovery / recomputation
         * --------------------------------------------------------------
         *
         * For this first Layer-2 version, recomputation is attempted only
         * when the entire graph remains connected.
         *
         * This prevents an all-pairs oblivious solver from being asked to
         * construct routes between disconnected components.
         */

        if (!failure.graph_disconnected) {
            failure.recomputation_attempted = true;

            auto failed_graph_result = makeFailureGraph(graph);

            if (failed_graph_result) {
                auto failed_graph_container = std::move(failed_graph_result.value());

                auto& failed_graph = failed_graph_container.graph;

                const auto recovery_start = timeNow();

                try {
                    Result<std::unique_ptr<RoutingScheme>> scheme_result = makeErrorMessage(ErrorCode::InvalidSolver, "FailureRecoveryAnalyzer: solver was not initialized.");

                    if (isSemiObliviousSolver(solver_type)) {
                        /*
                         * Semi-oblivious recovery:
                         *
                         * 1. rebuild candidate paths on the failed topology
                         * 2. optimize path loads for the same demand matrix
                         */
                        scheme_result = recomputeSemiOblivious(failed_graph, solver_type, demand_map, demand_model );

                    } else {

                        /*
                         * Standard oblivious recovery.
                         */
                        auto solver_option = makeSolver(solver_type, failed_graph);

                        if (solver_option && *solver_option) {
                            auto& solver = *solver_option;

                            scheme_result = solver->solve();
                        }
                    }

                    if (scheme_result && scheme_result.value()) {
                        failure.recomputation_runtime_microseconds = duration(timeNow() - recovery_start);

                        auto scheme = std::move(scheme_result.value());

                        failure.post_failure_congestion = evaluateCongestion(failed_graph, *scheme, demand_map);

                        if (std::isfinite(failure.post_failure_congestion)) {
                            failure.recomputation_succeeded = true;

                            ++analysis.summary.successful_recomputations;

                            if (baseline_congestion > recovery_epsilon) {
                                failure.congestion_increase_factor = failure.post_failure_congestion / baseline_congestion;
                            }
                        }
                    }
                } catch (...) {
                    /*
                     * A solver failure on one N-1 topology is itself a
                     * useful recovery result.
                     *
                     * Do not abort the complete failure sweep.
                     */
                }
            }

            if (!failure.recomputation_succeeded) {
                ++analysis.summary.failed_recomputations;
            }
        }

        /*
         * --------------------------------------------------------------
         * Aggregate successful recovery metrics
         * --------------------------------------------------------------
         */

        if (failure.recomputation_succeeded) {
            /*
             * --------------------------------------------------------------
             * Worst survivable post-failure congestion
             * --------------------------------------------------------------
             *
             * Only successful recomputations are considered here.
             */
            if (failure.post_failure_congestion >= 0.0 && (!have_worst_congestion || failure.post_failure_congestion > worst_congestion.post_failure_congestion)) {
                worst_congestion = failure;

                have_worst_congestion = true;
            }

            /*
             * --------------------------------------------------------------
             * Slowest successful recomputation
             * --------------------------------------------------------------
             */
            if (failure.recomputation_runtime_microseconds >= 0.0 && (!have_slowest_recovery || failure.recomputation_runtime_microseconds > slowest_recovery.recomputation_runtime_microseconds)) {
                slowest_recovery = failure;

                have_slowest_recovery = true;
            }
            /*
             * --------------------------------------------------------------
             * Aggregate post-failure congestion
             * --------------------------------------------------------------
             */
            if (failure.post_failure_congestion >= 0.0) {
                post_failure_congestion_sum += failure.post_failure_congestion;

                ++post_failure_congestion_count;

                analysis.summary.maximum_post_failure_congestion = std::max(analysis.summary.maximum_post_failure_congestion,

                                                                            failure.post_failure_congestion);
            }

            if (failure.congestion_increase_factor >= 0.0) {
                congestion_increase_sum += failure.congestion_increase_factor;

                ++congestion_increase_count;

                analysis.summary.maximum_congestion_increase_factor = std::max(analysis.summary.maximum_congestion_increase_factor,

                                                                               failure.congestion_increase_factor);
            }

            if (failure.recomputation_runtime_microseconds >= 0.0) {
                recomputation_runtime_sum += failure.recomputation_runtime_microseconds;

                ++recomputation_runtime_count;

                analysis.summary.maximum_recomputation_runtime_microseconds = std::max(analysis.summary.maximum_recomputation_runtime_microseconds,

                                                                                       failure.recomputation_runtime_microseconds);
            }
        }

        /*
         * --------------------------------------------------------------
         * Track worst physical-link failure
         * --------------------------------------------------------------
         */

        if (!have_worst_failure || worseFailure(failure, worst_failure)) {
            worst_failure = failure;

            have_worst_failure = true;
        }

        analysis.failures.push_back(failure);

        /*
         * --------------------------------------------------------------
         * Restore original topology
         * --------------------------------------------------------------
         */

        auto rollback = graph.rollback(checkpoint);

        if (!rollback) {
            return getError(rollback);
        }
    }

    /*
     * ------------------------------------------------------------------
     * Final aggregate statistics
     * ------------------------------------------------------------------
     */

    if (analysis.summary.tested_links > 0) {
        analysis.summary.average_unroutable_demand_fraction = unroutable_fraction_sum / static_cast<double>(analysis.summary.tested_links);
    }

    if (post_failure_congestion_count > 0) {
        analysis.summary.average_post_failure_congestion = post_failure_congestion_sum / static_cast<double>(post_failure_congestion_count);
    }

    if (congestion_increase_count > 0) {
        analysis.summary.average_congestion_increase_factor = congestion_increase_sum / static_cast<double>(congestion_increase_count);
    }

    if (recomputation_runtime_count > 0) {
        analysis.summary.average_recomputation_runtime_microseconds = recomputation_runtime_sum / static_cast<double>(recomputation_runtime_count);
    }

    if (have_worst_failure) {
        analysis.summary.worst_failed_edge_id = worst_failure.failed_edge_id;

        analysis.summary.worst_failed_source = worst_failure.source;

        analysis.summary.worst_failed_target = worst_failure.target;
    }
    /*
     * --------------------------------------------------------------
     * Worst disconnect
     * --------------------------------------------------------------
     */
    if (have_worst_disconnect) {
        analysis.summary.worst_disconnect_edge_id = worst_disconnect.failed_edge_id;

        analysis.summary.worst_disconnect_source = worst_disconnect.source;

        analysis.summary.worst_disconnect_target = worst_disconnect.target;

        analysis.summary.worst_disconnect_unroutable_demand_fraction = worst_disconnect.unroutable_demand_fraction;
    }

    /*
     * --------------------------------------------------------------
     * Worst survivable congestion
     * --------------------------------------------------------------
     */
    if (have_worst_congestion) {
        analysis.summary.worst_congestion_edge_id = worst_congestion.failed_edge_id;

        analysis.summary.worst_congestion_source = worst_congestion.source;

        analysis.summary.worst_congestion_target = worst_congestion.target;

        analysis.summary.worst_congestion_baseline = worst_congestion.baseline_congestion;

        analysis.summary.worst_congestion_post_failure = worst_congestion.post_failure_congestion;

        analysis.summary.worst_congestion_increase_factor = worst_congestion.congestion_increase_factor;
    }

    /*
     * --------------------------------------------------------------
     * Slowest successful recovery
     * --------------------------------------------------------------
     */
    if (have_slowest_recovery) {
        analysis.summary.slowest_recovery_edge_id = slowest_recovery.failed_edge_id;

        analysis.summary.slowest_recovery_source = slowest_recovery.source;

        analysis.summary.slowest_recovery_target = slowest_recovery.target;

        analysis.summary.slowest_recovery_runtime_microseconds = slowest_recovery.recomputation_runtime_microseconds;
    }
    return analysis;
}

Result<std::unique_ptr<RoutingScheme>> FailureRecoveryAnalyzer::recomputeSemiOblivious(optimized::Graph<EdgeData>& graph, const SolverType solver_type, const demands& demand_map,
                                                                                       const DemandModelType demand_type) {
    std::shared_ptr<SemiSolverRoutingEngine> routing_engine;

    switch (solver_type) {

    case SolverType::SEMI_ELECTRICAL:
        routing_engine = std::make_shared<SemiSolverRoutingEngine>(std::make_shared<ElectricalMWU>(graph, 0, true));
        break;

    case SolverType::SEMI_TREE:
        routing_engine = std::make_shared<SemiSolverRoutingEngine>(std::make_shared<TreeMWU<FlatHST>>(graph, 0, std::make_unique<FastCKR<FlatHST>>(graph)));
        break;

    case SolverType::SEMI_EXPANDER_HIERARCHY:
        routing_engine = std::make_shared<SemiSolverRoutingEngine>(std::make_shared<ElectrifiedExpanderHierarchySolver>(graph, 0));
        break;

    default:
        return makeErrorMessage(ErrorCode::InvalidSolver, "FailureRecoveryAnalyzer: requested semi recomputation for non-semi solver.");
    }

    auto optimizer = std::make_shared<OrToolsSemiObliviousLoadOptimizer>();

    SemiObliviousRoutingSolver solver(graph, routing_engine, optimizer);

    auto pre = solver.preprocess();

    if (!pre) {
        return getError(pre);
    }

    auto routed = solver.route(demand_map, demand_type);

    if (!routed) {
        return getError(routed);
    }

    auto result = std::move(routed.value());

    if (!result.scheme) {
        return makeErrorMessage(ErrorCode::InvalidRouting, "FailureRecoveryAnalyzer: semi-oblivious recovery returned no routing scheme.");
    }

    return std::move(result.scheme);
}