//
// Created by Mert Biyikli on 12.06.26.
//

#ifndef OBLIVIOUSROUTING_SEMI_ROUTING_ENGINE_H
#define OBLIVIOUSROUTING_SEMI_ROUTING_ENGINE_H

#include "../../data_structures/graph/graph.h"
#include "algorithms/oblivious/oblivious_solver.h"
#include "preprocessing/candidate_routing_scheme.h"
#include "core/errors.h"

class SemiSolverRoutingEngine  {
public:
    virtual ~SemiSolverRoutingEngine() = default;

    // TODO: in the future maybe make it also applicable for general oblivious routing as well
    explicit SemiSolverRoutingEngine(std::shared_ptr<ILinearObliviousSolverBase> solver)
        : solver_(std::move(solver)) {}

    Result<CandidateRoutingScheme> preprocess(const optimized::Graph<EdgeData>& graph);
    Result<void> extractPath(const optimized::Graph<EdgeData>& g,const RoutingScheme& scheme,int s,int t,CandidateRoutingScheme& out) const;
    virtual const std::string getSolverBase() const;
    bool dfsDecompose(const optimized::Graph<EdgeData>& g,int u,int t,std::vector<double>& residual,std::vector<int>& currentPath,std::vector<bool>& visited) const;


    std::shared_ptr<ILinearObliviousSolverBase> solver_;
};

#endif //OBLIVIOUSROUTING_SEMI_ROUTING_ENGINE_H