//
// Created by Mert Biyikli on 24.06.26.
//

#ifndef OBLIVIOUSROUTING_SEMI_ROUTING_RUNNER_H
#define OBLIVIOUSROUTING_SEMI_ROUTING_RUNNER_H

#include "routing/routing_runner.h"

static Result<std::shared_ptr<SemiSolverRoutingEngine>> makeSemiRoutingEngine(SolverType type,IGraph& graph);

class SemiObliviousSolverRunner final : public IRoutingExperimentRunner {
public:
    Result<IRoutingResult> run(IGraph& graph,const Config& cfg,SolverType type) const override ;
};
#endif //OBLIVIOUSROUTING_SEMI_ROUTING_RUNNER_H