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
    static Result<void> evaluate(optimized::Graph<EdgeData>& graph, const std::unique_ptr<RoutingScheme>& scheme, const Config& cfg, IRoutingResult& result, bool run_failure_analysis = false, double demand_scale = 1.0) {
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


            if (!std::isfinite(demand_scale) || demand_scale <= 0.0) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "Demand scale must be finite and positive");
            }

            demands scaledDemands;

            for (std::size_t i = 0; i < dmap->size(); ++i) {
                const auto [source, target] = dmap->getDemandPair(i);
                const double value = dmap->getDemandValue(i) * demand_scale;

                if (!std::isfinite(value)) {
                    return makeErrorMessage(ErrorCode::InvalidDemand, "Scaled demand overflow");
                }

                scaledDemands.addDemand(source, target, value);
            }

            dmap = std::move(scaledDemands);


            /*
             * RoutingAnalyzer now measures:
             *
             * 1. normal demand evaluation runtime
             * 2. failure-analysis runtime
             *
             * independently.
             */
            auto visualization = RoutingAnalyzer::analyze(graph, *scheme, dmap.value(), result, demandModelName(demandType), run_failure_analysis);

            if (!visualization) {
                return getError(visualization);
            }

            const auto& analysis = visualization.value();

            FailureRecoveryAnalysis recovery_analysis;
            bool recovery_available = false;

            // Layer 2: Optional N-1 recovery analysis.
            if (cfg.failure_recovery) {
                auto recovery = FailureRecoveryAnalyzer::analyze(graph, result.type, dmap.value(), demandType, analysis.summary.maximum_congestion);

                if (!recovery) {
                    return getError(recovery);
                }

                recovery_analysis = std::move(recovery.value());
                recovery_available = true;
            }

            // Always record the demand evaluation.
            DemandEvaluationResult evaluation{};

            evaluation.demand_type = demandType;
            evaluation.congestion = analysis.summary.maximum_congestion;
            evaluation.runtime_microseconds = analysis.demand_evaluation_runtime_microseconds;

            // Only populate Layer-1 failure metrics when requested.
            if (cfg.failure_analysis) {
                const auto& failure_summary = analysis.failure_summary;

                evaluation.failure_analysis_runtime_microseconds = analysis.failure_analysis_runtime_microseconds;
                evaluation.failure_tested_links = failure_summary.tested_links;
                evaluation.failure_most_critical_edge_id = failure_summary.most_critical_edge_id;
                evaluation.failure_most_critical_source = failure_summary.most_critical_source;
                evaluation.failure_most_critical_target = failure_summary.most_critical_target;
                evaluation.failure_maximum_lost_traffic_fraction = failure_summary.maximum_lost_traffic_fraction;
                evaluation.failure_average_lost_traffic_fraction = failure_summary.average_lost_traffic_fraction;
                evaluation.failure_median_lost_traffic_fraction = failure_summary.median_lost_traffic_fraction;
                evaluation.failure_maximum_affected_demand_fraction = failure_summary.maximum_affected_demand_fraction;
                evaluation.failure_traffic_carrying_links = failure_summary.traffic_carrying_links;
                evaluation.failure_critical_links_10_percent = failure_summary.critical_links_10_percent;
                evaluation.failure_critical_links_25_percent = failure_summary.critical_links_25_percent;
                evaluation.failure_critical_links_50_percent = failure_summary.critical_links_50_percent;
            }

            // Preserve the existing recovery assignments from your original
            // initializer, using recovery_available and recovery_analysis.
            evaluation.failure_recovery_available = recovery_available;

            if (recovery_available) {
                evaluation.recovery_events = std::move(recovery_analysis.failures);
            }

            result.demand_evaluations.push_back(std::move(evaluation));
            result.visualization_results.push_back(std::move(visualization.value()));
        }

        return {};
    }
};

#endif // OBLIVIOUSROUTING_ROUTING_RUNNER_H
