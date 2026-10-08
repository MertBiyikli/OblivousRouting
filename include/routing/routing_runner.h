//
// Created by Mert Biyikli on 24.06.26.
//

#ifndef OBLIVIOUSROUTING_ROUTING_RUNNER_H
#define OBLIVIOUSROUTING_ROUTING_RUNNER_H

#include "core/errors.h"
#include "routing_engine.h"
#include "visualization/failure_recover_analysis.h"
#include "visualization/visualization_result.h"

class IRoutingExperimentRunner {
  public:
    virtual ~IRoutingExperimentRunner() = default;

    virtual Result<IRoutingResult> run(optimized::Graph<EdgeData>& graph, const Config& cfg, SolverType type) const = 0;
};

class DemandEvaluator {
  public:
    static Result<void> evaluate(optimized::Graph<EdgeData>& graph, const std::unique_ptr<RoutingScheme>& scheme, const Config& cfg, IRoutingResult& result) {
        if (!cfg.evaluate_demand_models) {
            return makeErrorMessage(ErrorCode::InvalidDemand, "Evaluating demand model is set off.");
        }

        if (!scheme) {
            result.status = ResultStatus::ERROR_INVALID_ROUTING_SCHEME;

            return makeErrorMessage(ErrorCode::InvalidRouting, "Routing scheme is invalid, when evaluating congestion.");
        }

        auto pairs = generateAllDemandPairs(graph);

        for (auto demandType : cfg.demand_models) {
            auto model = makeDemandModel(demandType);

            auto dmap = model->generate(graph, pairs);

            if (!dmap) {
                return getError(dmap);
            }

            /*
             * RoutingAnalyzer now measures:
             *
             * 1. normal demand evaluation runtime
             * 2. failure-analysis runtime
             *
             * independently.
             */
            auto visualization = RoutingAnalyzer::analyze(graph, *scheme, dmap.value(), result, demandModelName(demandType));

            if (!visualization) {
                return getError(visualization);
            }

            const auto& analysis = visualization.value();

            const auto& failure_summary = analysis.failure_summary;

            /*
             * ------------------------------------------------------------
             * Optional Layer-2 N-1 recovery analysis
             * ------------------------------------------------------------
             */

            FailureRecoveryAnalysis recovery_analysis;

            bool recovery_available = false;

            if (cfg.failure_recovery) {
                auto recovery = FailureRecoveryAnalyzer::analyze(graph, result.type, dmap.value(),  demandType, analysis.summary.maximum_congestion);

                if (!recovery) {
                    return getError(recovery);
                }

                recovery_analysis = std::move(recovery.value());

                recovery_available = true;
            }

            result.demand_evaluations.emplace_back(DemandEvaluationResult{
                .demand_type = demandType,

                .congestion = analysis.summary.maximum_congestion,

                /*
                 * Demand routing + congestion calculation only.
                 */
                .runtime_microseconds = analysis.demand_evaluation_runtime_microseconds,

                /*
                 * Static N-1 failure analysis only.
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


            if (recovery_available) {
                result.demand_evaluations.back().recovery_events = std::move(recovery_analysis.failures);
            }


            result.visualization_results.push_back(std::move(visualization.value()));
        }

        return {};
    }
};

#endif // OBLIVIOUSROUTING_ROUTING_RUNNER_H
