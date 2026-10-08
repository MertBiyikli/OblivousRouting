//
// Created by Mert Biyikli on 24.06.26.
//

#include "algorithms/semi_oblivious/semi_routing_runner.h"
#include "algorithms/oblivious/oblivious_routing_runner.h"
#include "core/errors.h"

static Result<std::shared_ptr<SemiSolverRoutingEngine>> makeSemiRoutingEngine(SolverType type, optimized::Graph<EdgeData>& graph) {
    switch (type) {
    case SolverType::SEMI_ELECTRICAL:
        return std::make_shared<SemiSolverRoutingEngine>(std::make_shared<ElectricalMWU>(graph, 0, true));

    case SolverType::SEMI_TREE:
        return std::make_shared<SemiSolverRoutingEngine>(std::make_shared<TreeMWU<FlatHST>>(graph, 0, std::make_unique<FastCKR<FlatHST>>(graph)));

    case SolverType::SEMI_EXPANDER_HIERARCHY:
        return std::make_shared<SemiSolverRoutingEngine>(std::make_shared<ElectrifiedExpanderHierarchySolver>(graph, 0));

    default:
        return makeErrorMessage(ErrorCode::InvalidSolver, "Requested semi-oblivious routing engine for non-semi solver");
    }
}

Result<IRoutingResult> SemiObliviousSolverRunner::run(optimized::Graph<EdgeData>& graph, const Config& cfg, SolverType type) const {
    SemiObliviousRoutingResult semiResult;
    semiResult.type = type;
    semiResult.graph_name = cfg.filename;
    semiResult.nodes = graph.getNumNodes();
    semiResult.edges = graph.getNumUndirectedEdges();
    semiResult.solver_name = getSolverName(type);

    if (!cfg.evaluate_demand_models || cfg.demand_models.empty()) {
        semiResult.status = ResultStatus::ERROR_MISSING_DEMAND_MODELS;
        return makeErrorMessage(ErrorCode::InvalidDemand, "Semi-oblivious solver requires demand models.");
    }

    auto engine_factory = makeSemiRoutingEngine(type, graph);
    auto optimizer_factory = std::make_shared<OrToolsSemiObliviousLoadOptimizer>();

    if (!optimizer_factory || !engine_factory) {
        return getError(engine_factory);
    }

    std::shared_ptr<SemiSolverRoutingEngine> routingEngine = engine_factory.value();
    std::shared_ptr<OrToolsSemiObliviousLoadOptimizer> optimizer = optimizer_factory;

    SemiObliviousRoutingSolver solver(graph, routingEngine, optimizer);

    const auto preprocessStart = timeNow();

    auto pre = solver.preprocess();
    if (!pre) {
        return getError(pre);
    }
    semiResult.preprocessing_runtime_microseconds = duration(timeNow() - preprocessStart);

    auto pairs = generateAllDemandPairs(graph);

    for (const auto& demandType : cfg.demand_models) {
        auto model = makeDemandModel(demandType);

        auto dmap = model->generate(graph, pairs);

        if (!dmap) {
            return getError(dmap);
        }

        /*
         * ------------------------------------------------------------
         * Demand-specific semi-oblivious solve
         * ------------------------------------------------------------
         *
         * SemiObliviousRoutingSolver::route() already measures the
         * actual load-optimization runtime separately and stores it in
         *
         *     result.total_runtime_microseconds
         *
         * RoutingAnalyzer then separately records:
         *
         *     demand_evaluation_runtime_microseconds
         *     failure_analysis_runtime_microseconds
         *
         * Therefore we do NOT time the complete route() call here.
         */

        SemiObliviousRoutingResult result;

        auto routed = solver.route(dmap.value(), demandType);

        if (!routed) {
            return getError(routed);
        }

        result = std::move(routed.value());

        /*
         * Demand-specific optimization runtime only.
         */
        semiResult.solve_runtime_microseconds += result.total_runtime_microseconds;

        semiResult.routing_base = result.path_selection_strategy;

        semiResult.candidate_paths = result.candidate_paths;

        semiResult.average_paths_per_pair = result.average_paths_per_pair;

        /*
         * ------------------------------------------------------------
         * Pull analysis metrics from the exact demand-specific scheme.
         * ------------------------------------------------------------
         */

        const auto& analysis = result.visualization;

        const auto& failure_summary = analysis.failure_summary;

        FailureRecoveryAnalysis recovery_analysis;

        bool recovery_available = false;

        if (cfg.failure_recovery) {
            auto recovery = FailureRecoveryAnalyzer::analyze(graph, type, dmap.value(), demandType, analysis.summary.maximum_congestion);

            if (!recovery) {
                return getError(recovery);
            }

            recovery_analysis = std::move(recovery.value());

            recovery_available = true;
        }
        semiResult.demand_evaluations.emplace_back(DemandEvaluationResult{
            .demand_type = demandType,

            .congestion = analysis.summary.maximum_congestion,

            /*
             * Baseline demand evaluation only.
             *
             * Does not include semi-oblivious optimization.
             * Does not include failure analysis.
             */
            .runtime_microseconds = analysis.demand_evaluation_runtime_microseconds,

            /*
             * ------------------------------------------------------------
             * Layer 1: static single-link exposure analysis
             * ------------------------------------------------------------
             */
            .failure_analysis_runtime_microseconds = analysis.failure_analysis_runtime_microseconds,

            .failure_tested_links = failure_summary.tested_links,

            .failure_most_critical_edge_id = failure_summary.most_critical_edge_id,

            .failure_most_critical_source = failure_summary.most_critical_source,

            .failure_most_critical_target = failure_summary.most_critical_target,

            .failure_maximum_lost_traffic_fraction = failure_summary.maximum_lost_traffic_fraction,

            .failure_average_lost_traffic_fraction = failure_summary.average_lost_traffic_fraction,

            .failure_median_lost_traffic_fraction = failure_summary.median_lost_traffic_fraction,

            .failure_maximum_affected_demand_fraction = failure_summary.maximum_affected_demand_fraction,

            .failure_traffic_carrying_links = failure_summary.traffic_carrying_links,

            .failure_critical_links_10_percent = failure_summary.critical_links_10_percent,

            .failure_critical_links_25_percent = failure_summary.critical_links_25_percent,

            .failure_critical_links_50_percent = failure_summary.critical_links_50_percent,

            /*
             * ------------------------------------------------------------
             * Layer 2: actual N-1 recovery analysis
             * ------------------------------------------------------------
             */
            .failure_recovery_available = recovery_available,

            .recovery_tested_links = recovery_available ? recovery_analysis.summary.tested_links : 0,

            .recovery_disconnected_failures = recovery_available ? recovery_analysis.summary.disconnected_failures : 0,

            .recovery_successful_recomputations = recovery_available ? recovery_analysis.summary.successful_recomputations : 0,

            .recovery_failed_recomputations = recovery_available ? recovery_analysis.summary.failed_recomputations : 0,

            .recovery_maximum_unroutable_demand_fraction = recovery_available ? recovery_analysis.summary.maximum_unroutable_demand_fraction : -1.0,

            .recovery_average_unroutable_demand_fraction = recovery_available ? recovery_analysis.summary.average_unroutable_demand_fraction : -1.0,

            .recovery_maximum_post_failure_congestion = recovery_available ? recovery_analysis.summary.maximum_post_failure_congestion : -1.0,

            .recovery_average_post_failure_congestion = recovery_available ? recovery_analysis.summary.average_post_failure_congestion : -1.0,

            .recovery_maximum_congestion_increase_factor = recovery_available ? recovery_analysis.summary.maximum_congestion_increase_factor : -1.0,

            .recovery_average_congestion_increase_factor = recovery_available ? recovery_analysis.summary.average_congestion_increase_factor : -1.0,

            .recovery_average_recomputation_runtime_microseconds = recovery_available ? recovery_analysis.summary.average_recomputation_runtime_microseconds : -1.0,

            .recovery_maximum_recomputation_runtime_microseconds = recovery_available ? recovery_analysis.summary.maximum_recomputation_runtime_microseconds : -1.0,

            /*
             * Overall worst failure.
             */
            .recovery_worst_failed_edge_id = recovery_available ? recovery_analysis.summary.worst_failed_edge_id : -1,

            .recovery_worst_failed_source = recovery_available ? recovery_analysis.summary.worst_failed_source : -1,

            .recovery_worst_failed_target = recovery_available ? recovery_analysis.summary.worst_failed_target : -1,

            /*
             * Worst disconnecting failure.
             */
            .recovery_worst_disconnect_edge_id = recovery_available ? recovery_analysis.summary.worst_disconnect_edge_id : -1,

            .recovery_worst_disconnect_source = recovery_available ? recovery_analysis.summary.worst_disconnect_source : -1,

            .recovery_worst_disconnect_target = recovery_available ? recovery_analysis.summary.worst_disconnect_target : -1,

            .recovery_worst_disconnect_unroutable_demand_fraction = recovery_available ? recovery_analysis.summary.worst_disconnect_unroutable_demand_fraction : -1.0,

            /*
             * Worst survivable congestion failure.
             */
            .recovery_worst_congestion_edge_id = recovery_available ? recovery_analysis.summary.worst_congestion_edge_id : -1,

            .recovery_worst_congestion_source = recovery_available ? recovery_analysis.summary.worst_congestion_source : -1,

            .recovery_worst_congestion_target = recovery_available ? recovery_analysis.summary.worst_congestion_target : -1,

            .recovery_worst_congestion_baseline = recovery_available ? recovery_analysis.summary.worst_congestion_baseline : -1.0,

            .recovery_worst_congestion_post_failure = recovery_available ? recovery_analysis.summary.worst_congestion_post_failure : -1.0,

            .recovery_worst_congestion_increase_factor = recovery_available ? recovery_analysis.summary.worst_congestion_increase_factor : -1.0,

            /*
             * Slowest successful recovery.
             */
            .recovery_slowest_recovery_edge_id = recovery_available ? recovery_analysis.summary.slowest_recovery_edge_id : -1,

            .recovery_slowest_recovery_source = recovery_available ? recovery_analysis.summary.slowest_recovery_source : -1,

            .recovery_slowest_recovery_target = recovery_available ? recovery_analysis.summary.slowest_recovery_target : -1,

            .recovery_slowest_recovery_runtime_microseconds = recovery_available ? recovery_analysis.summary.slowest_recovery_runtime_microseconds : -1.0,
        });

        /*
         * Preserve the exact demand-specific routing scheme.
         */
        semiResult.demand_schemes[demandModelName(demandType)] = std::move(result.scheme);

        /*
         * Preserve visualization / analysis result.
         */
        RoutingVisualizationResult visualization = std::move(result.visualization);

        visualization.graph_name = cfg.filename;

        visualization.solver_name = getSolverName(type);

        semiResult.visualization_results.push_back(std::move(visualization));
    }

    semiResult.total_runtime_microseconds = semiResult.preprocessing_runtime_microseconds + semiResult.solve_runtime_microseconds;
    semiResult.oblivious_ratio = -1;

    ObliviousSolverRunner::appendMetricsIfAvailable(routingEngine->solver_, semiResult);
    // ObliviousSolverRunner::appendObjectiveIfAvailable(routingEngine, *semiResult.scheme, semiResult);

    semiResult.status = ResultStatus::OK;
    return semiResult;
}