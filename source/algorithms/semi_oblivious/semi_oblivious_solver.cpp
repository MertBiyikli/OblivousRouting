//
// Created by Mert Biyikli on 24.06.26.
//

#include "algorithms/semi_oblivious/semi_oblivious_solver.h"

#include "utils/time_tracking.h"


std::unique_ptr<RoutingScheme> SemiObliviousRoutingSolver::solve() {
    current_result = SemiObliviousRoutingResult();
    if (!demand_ || !demandType_) {
        throw std::logic_error(
            "SemiObliviousRoutingSolver::setDemand must be called before solve"
        );
    }

    if (!candidateScheme_) {
        candidateScheme_ = routingEngine_->preprocess(graph);
    }

    current_result = route(*demand_, *demandType_);

    return std::move(current_result.scheme); // NOTE: after that the scheme in the current result is empty/null
}


void SemiObliviousRoutingSolver::setDemand(const demands& demand, DemandModelType demandType) {
    demand_ = demand;
    demandType_ = demandType;
}

CandidateRoutingScheme SemiObliviousRoutingSolver::preprocess() {
    candidateScheme_ = routingEngine_->preprocess(graph);
    return *candidateScheme_;
}


SemiObliviousRoutingResult SemiObliviousRoutingSolver::route(const demands& demand,DemandModelType demandType) {
    setDemand(demand, demandType);

    if (!candidateScheme_) {
        candidateScheme_ = routingEngine_->preprocess(graph);
    }

    SemiObliviousRoutingResult res;
    res.demand_type = demandType;
    res.path_selection_strategy = routingEngine_->getSolverBase();

    res.candidate_paths = candidateScheme_->numPaths();
    res.average_paths_per_pair =
        candidateScheme_->averagePathsPerPair(graph.getNumNodes());

    const auto start = timeNow();

    auto optResult = loadOptimizer_->optimize(
        graph,
        *candidateScheme_,
        demand
    );

    res.total_runtime_microseconds = duration(timeNow() - start);
    res.scheme = std::move(optResult.scheme);

    std::vector<double> cong;
    res.scheme->routeDemands(cong, demand);

    res.congestion = 0.0;
    for (double c : cong) {
        res.congestion = std::max(res.congestion, c);
    }

    return res;
}