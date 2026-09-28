//
// Created by Mert Biyikli on 24.06.26.
//

#ifndef OBLIVIOUSROUTING_SEMI_ROUTING_RUNNER_H
#define OBLIVIOUSROUTING_SEMI_ROUTING_RUNNER_H

#include "routing/routing_runner.h"

static Result<std::shared_ptr<SemiSolverRoutingEngine>> makeSemiRoutingEngine(SolverType type, optimized::Graph<EdgeData>& graph);

class SemiObliviousSolverRunner final : public IRoutingExperimentRunner {
public:
    Result<IRoutingResult> run(optimized::Graph<EdgeData>& graph,const Config& cfg,SolverType type) const override ;
};
#endif //OBLIVIOUSROUTING_SEMI_ROUTING_RUNNER_H