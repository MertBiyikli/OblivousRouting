//
// Created by Mert Biyikli on 24.06.26.
//

#ifndef OBLIVIOUSROUTING_OBLIVIOUS_ROUTING_RUNNER_H
#define OBLIVIOUSROUTING_OBLIVIOUS_ROUTING_RUNNER_H

#include "routing/routing_runner.h"


class ObliviousSolverRunner final : public IRoutingExperimentRunner {
public:
    IRoutingResult run(IGraph& graph,const Config& cfg,SolverType type) const override;

private:
    static void appendMetricsIfAvailable(const std::unique_ptr<ISolver>& solver,IRoutingResult& result) {
        if (const auto* mwu = dynamic_cast<const MWUFramework*>(solver.get())) {
            result.mwu_metrics = mwu->getMetrics();
        }
    }

    static void appendObjectiveIfAvailable(const std::unique_ptr<ISolver>& solver,const RoutingScheme& scheme,IRoutingResult& result) {
        if (const auto* linearScheme =
                dynamic_cast<const LinearRoutingScheme*>(&scheme)) {
            result.oblivious_ratio = linearScheme->computeObliviousRatio();
                }

        if (const auto* lp = dynamic_cast<const LP*>(solver.get())) {
            if (lp->alpha) {
                result.oblivious_ratio = lp->alpha->solution_value();
            }
        }
    }
};

#endif //OBLIVIOUSROUTING_OBLIVIOUS_ROUTING_RUNNER_H