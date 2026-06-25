//
// Created by Mert Biyikli on 24.06.26.
//

#ifndef OBLIVIOUSROUTING_ROUTING_RUNNER_H
#define OBLIVIOUSROUTING_ROUTING_RUNNER_H

#include "routing_engine.h"



class IRoutingExperimentRunner {
public:
    virtual ~IRoutingExperimentRunner() = default;

    virtual IRoutingResult run(
        IGraph& graph,
        const Config& cfg,
        SolverType type
    ) const = 0;
};


class DemandEvaluator {
public:
    static void evaluate(
        IGraph& graph,
        const std::unique_ptr<RoutingScheme>& scheme,
        const Config& cfg,
        IRoutingResult& result
    ) {
        if (!cfg.evaluate_demand_models) {
            return;
        }

        if (!scheme) {
            result.status = ResultStatus::ERROR_INVALID_ROUTING_SCHEME;
            return;
        }

        auto pairs = generateAllDemandPairs(graph);

        for (auto demandType : cfg.demand_models) {
            auto model = makeDemandModel(demandType);
            demands dmap = model->generate(graph, pairs);

            const auto t0 = timeNow();
            double congestion =
                computeRoutingSchemeCongestion(graph, scheme, dmap);
            double time = duration(timeNow() - t0);
            result.demand_evaluations.emplace_back(
                DemandEvaluationResult{
                .demand_type = demandType,
                .congestion = congestion,
                .runtime_microseconds = time
            });
        }
    }
};



#endif //OBLIVIOUSROUTING_ROUTING_RUNNER_H