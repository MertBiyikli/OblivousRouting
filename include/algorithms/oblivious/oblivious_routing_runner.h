//
// Created by Mert Biyikli on 24.06.26.
//

#ifndef OBLIVIOUSROUTING_OBLIVIOUS_ROUTING_RUNNER_H
#define OBLIVIOUSROUTING_OBLIVIOUS_ROUTING_RUNNER_H

#include "routing/routing_runner.h"


class ObliviousSolverRunner final : public IRoutingExperimentRunner {
public:
    Result<IRoutingResult> run(IGraph& graph,const Config& cfg,SolverType type) const override;

    template<typename SolverPtr>
    static void appendMetricsIfAvailable(const SolverPtr& solver,IRoutingResult& result) {
        if (const auto* mwu = dynamic_cast<const MWUFramework*>(solver.get())) {
            result.mwu_metrics = mwu->getMetrics();
        }
        if (const auto* expander =
            dynamic_cast<const ElectrifiedExpanderHierarchySolver*>(
                solver.get()
            )) {
            result.expander_metrics =
                expander->getMetrics();

            /*
             * Give the generic runtime fields meaningful EEH values.
             */
            result.preprocessing_runtime_microseconds =
                result.expander_metrics.preprocessingRuntime();

            result.solve_runtime_microseconds =
                result.expander_metrics
                    .basis_flow_runtime_microseconds;
            }
    }

    template<typename SolverPtr>
    static void appendObjectiveIfAvailable(const SolverPtr& solver,const RoutingScheme& scheme,IRoutingResult& result) {
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