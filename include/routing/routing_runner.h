//
// Created by Mert Biyikli on 24.06.26.
//

#ifndef OBLIVIOUSROUTING_ROUTING_RUNNER_H
#define OBLIVIOUSROUTING_ROUTING_RUNNER_H

#include "routing_engine.h"
#include "core/errors.h"



class IRoutingExperimentRunner {
public:
    virtual ~IRoutingExperimentRunner() = default;

    virtual Result<IRoutingResult> run(
        IGraph& graph,
        const Config& cfg,
        SolverType type
    ) const = 0;
};


class DemandEvaluator {
public:
    static Result<void> evaluate(
        IGraph& graph,
        const std::unique_ptr<RoutingScheme>& scheme,
        const Config& cfg,
        IRoutingResult& result
    ) {
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

            const auto t0 = timeNow();
            double congestion =
                computeRoutingSchemeCongestion(graph, scheme, dmap.value());
            double time = duration(timeNow() - t0);

            result.demand_evaluations.emplace_back(
                DemandEvaluationResult{
                .demand_type = demandType,
                .congestion = congestion,
                .runtime_microseconds = time
            });
        }
        return {};
    }
};



#endif //OBLIVIOUSROUTING_ROUTING_RUNNER_H