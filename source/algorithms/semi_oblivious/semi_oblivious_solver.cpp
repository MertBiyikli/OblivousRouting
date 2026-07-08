//
// Created by Mert Biyikli on 24.06.26.
//

#include "algorithms/semi_oblivious/semi_oblivious_solver.h"

#include "utils/time_tracking.h"


Result<std::unique_ptr<RoutingScheme>> SemiObliviousRoutingSolver::solve() {
    current_result = SemiObliviousRoutingResult();
    if (!demand_ || !demandType_) {
        return makeErrorMessage(ErrorCode::LogicError, "Semi oblivious :demand must be set before solve.");
    }

    if (!candidateScheme_) {;
        if (auto pre = this->preprocess()) {
            candidateScheme_ = pre.value();
        }else {
            return getError(pre);
        }
    }

    if (auto routed = route(*demand_, *demandType_)) {
        current_result = std::move(routed.value());
    }else {
        getError(routed);
    }

    return std::move(current_result.scheme); // NOTE: after that the scheme in the current result is empty/null
}


void SemiObliviousRoutingSolver::setDemand(const demands& demand,  const DemandModelType& demandType) {
    demand_ = demand;
    demandType_ = demandType;
}

Result<CandidateRoutingScheme> SemiObliviousRoutingSolver::preprocess() {
    if (auto candidateScheme_res = routingEngine_->preprocess(graph)) {
        candidateScheme_ = candidateScheme_res.value();
        return *(candidateScheme_);
    }else {
        return getError(candidateScheme_res);
    }
}


Result<SemiObliviousRoutingResult> SemiObliviousRoutingSolver::route(const demands& demand, const DemandModelType& demandType) {
    if (!candidateScheme_) {
        return makeErrorMessage(ErrorCode::LogicError, "Semi oblivious: Paths must be precomputed before routing any demand.");
    }

    if (!demand.size()) {
        return makeErrorMessage(ErrorCode::InvalidDemand, "Semi oblivious: Demand must be set before routing.");
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

    if (!optResult) {
        return getError(optResult);
    }

    res.total_runtime_microseconds = duration(timeNow() - start);
    if (!optResult.value().scheme) {
        return makeErrorMessage(ErrorCode::InvalidRouting, "Semi oblivious: failed to compute demand-specific routing scheme.");
    }
    res.scheme = std::move(optResult.value().scheme);

    std::vector<double> cong;
    res.scheme->routeDemands(cong, demand);

    res.congestion = 0.0;
    for (double c : cong) {
        res.congestion = std::max(res.congestion, c);
    }

    return res;
}