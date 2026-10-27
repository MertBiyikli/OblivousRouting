//
// Created by Mert Biyikli on 09.06.26.
//


#pragma once

#include "../core/errors.h"
#include "algorithms/oblivious/oblivious_routing_runner.h"
#include "algorithms/semi_oblivious/semi_routing_runner.h"
#include "io/parse_argument_io.h"
#include "routing_result.h"

#include <optional>

class RoutingEngine {
public:

    Result<RoutingExperimentResult> entry(Config& cfg);


    Result<IRoutingResult> solve(optimized::Graph<EdgeData>& graph,const Config& cfg,SolverType type);

    Config cfg;
    std::unique_ptr<optimized::Graph<EdgeData>> m_graph;
    Result<optimized::Graph<EdgeData>*> getGraph() {
        if (!m_graph) {
            return makeErrorMessage(ErrorCode::RuntimeError, "RoutingEngine::getGraph: graph not set");
        }else {
            return m_graph.get();
        }
    }


    Config getConfig() {return cfg;}
};


inline bool isSemiObliviousSolver(SolverType type) {
    return type == SolverType::SEMI_ELECTRICAL ||
           type == SolverType::SEMI_TREE ||
           type == SolverType::SEMI_EXPANDER_HIERARCHY;
}

inline std::unique_ptr<IRoutingExperimentRunner> makeRunner(SolverType type) {
    if (isSemiObliviousSolver(type)) {
        return std::make_unique<SemiObliviousSolverRunner>();
    }
    return std::make_unique<ObliviousSolverRunner>();
}