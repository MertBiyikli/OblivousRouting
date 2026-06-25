//
// Created by Mert Biyikli on 23.06.26.
//

//
// Created by Mert Biyikli on 23.06.26.
//
#include "routing/routing_engine.h"
#include "routing/routing_runner.h"
#include "algorithms/oblivious/oblivious_routing_runner.h"
#include "algorithms/semi_oblivious/semi_routing_runner.h"

#include "algorithms/oblivious/mwu/electrical_mwu.h"
#include "algorithms/oblivious/mwu/tree_mwu.h"
#include "algorithms/semi_oblivious/semi_oblivious_solver.h"


bool isSemiObliviousSolver(SolverType type) {
    return type == SolverType::SEMI_ELECTRICAL ||
           type == SolverType::SEMI_TREE;
}

std::unique_ptr<IRoutingExperimentRunner> makeRunner(SolverType type) {
    if (isSemiObliviousSolver(type)) {
        return std::make_unique<SemiObliviousSolverRunner>();
    }

    return std::make_unique<ObliviousSolverRunner>();
}

std::optional<IRoutingResult> RoutingEngine::solve(
    IGraph& graph,
    const Config& cfg,
    SolverType type
) {
    try {
        auto runner = makeRunner(type);
        return runner->run(graph, cfg, type);
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] RoutingEngine failed: " << e.what() << '\n';

        IRoutingResult result;
        result.type = type;
        result.solver_name = getSolverName(type);
        result.status = ResultStatus::ERROR_EXCEPTION;
        return result;
    }
}




