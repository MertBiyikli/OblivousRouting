//
// Created by Mert Biyikli on 24.06.26.
//

#ifndef OBLIVIOUSROUTING_ROUTING_RUNNER_H
#define OBLIVIOUSROUTING_ROUTING_RUNNER_H

#include "analysis/scenario_engine.h"
#include "core/errors.h"
#include "io/parse_argument_io.h"
#include "routing_result.h"

#include "analysis/visualization_result.h"

class IRoutingExperimentRunner {
  public:
    virtual ~IRoutingExperimentRunner() = default;

    virtual Result<IRoutingResult> run(optimized::Graph<EdgeData>& graph, const Config& cfg, SolverType type) const = 0;
};

class DemandEvaluator {
  public:
    static Result<void> evaluate(optimized::Graph<EdgeData>& graph, const std::unique_ptr<RoutingScheme>& scheme, const Config& cfg, IRoutingResult& result, bool run_failure_analysis = false,
                                 double demand_scale = 1.0) {
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

            const auto evaluation_start = std::chrono::steady_clock::now();

            auto evaluated = DemandEvaluator::evaluate(graph, result, dmap.value());

            if (!evaluated) {
                return getError(evaluated);
            }

            const auto evaluation_end = std::chrono::steady_clock::now();

            DemandEvaluationResult summary{};

            summary.demand_type = demandType;
            summary.congestion = evaluated->max_congestion;

            summary.runtime_microseconds = std::chrono::duration<double, std::micro>(evaluation_end - evaluation_start).count();



            result.demand_evaluations.push_back(std::move(summary));
        }

        return {};
    }

    static Result<NetworkData> evaluate(const optimized::Graph<EdgeData>& graph, const IRoutingResult& routing, const demands& demand) {
        double total_volume = 0.0;

        for (std::size_t i = 0; i < demand.size(); ++i) {
            const double value = demand.getDemandValue(i);

            if (!std::isfinite(value) || value < 0.0) {
                return makeErrorMessage(ErrorCode::InvalidDemand, "Demand must be finite and nonnegative.");
            }

            total_volume += value;

            if (!std::isfinite(total_volume)) {
                return makeErrorMessage(ErrorCode::InvalidDemand, "Total demand volume overflow.");
            }
        }

        return evaluateRouting(graph, routing, demand, total_volume);
    }
};

#endif // OBLIVIOUSROUTING_ROUTING_RUNNER_H
