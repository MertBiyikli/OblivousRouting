//
// Created by Mert Biyikli on 24.06.26.
//

#include "algorithms/semi_oblivious/semi_routing_runner.h"

static std::shared_ptr<SemiSolverRoutingEngine> makeSemiRoutingEngine(SolverType type,IGraph& graph) {
    switch (type) {
        case SolverType::SEMI_ELECTRICAL:
            return std::make_shared<SemiSolverRoutingEngine>(
                std::make_shared<ElectricalMWU>(graph, 0, true)
            );

        case SolverType::SEMI_TREE:
            return std::make_shared<SemiSolverRoutingEngine>(
                std::make_shared<TreeMWU<FlatHST>>(
                    graph,
                    0,
                    std::make_unique<FastCKR<FlatHST>>(graph)
                )
            );

        default:
            throw std::invalid_argument(
                "Requested semi-oblivious routing engine for non-semi solver"
            );
    }
}

IRoutingResult SemiObliviousSolverRunner::run(IGraph& graph,const Config& cfg,SolverType type) const {
    IRoutingResult result;
    result.type = type;
    result.graph_name = cfg.filename;
    result.solver_name = getSolverName(type);

    if (!cfg.evaluate_demand_models || cfg.demand_models.empty()) {
        std::cerr << "[ERROR] Semi-oblivious solver requires demand models.\n";
        result.status = ResultStatus::ERROR_MISSING_DEMAND_MODELS;
        return result;
    }

    auto routingEngine = makeSemiRoutingEngine(type, graph);
    auto optimizer =
        std::make_shared<OrToolsSemiObliviousLoadOptimizer>();

    SemiObliviousRoutingSolver solver(
        graph,
        routingEngine,
        optimizer
    );

    const auto preprocessStart = timeNow();
    solver.preprocess();
    result.preprocessing_runtime_microseconds =
        duration(timeNow() - preprocessStart);

    auto pairs = generateAllDemandPairs(graph);

    for (const auto demandType : cfg.demand_models) {
        auto model = makeDemandModel(demandType);
        demands dmap = model->generate(graph, pairs);

        const auto t0 = timeNow();
        auto semiResult = solver.route(dmap, demandType);
        const double runtime = duration(timeNow() - t0);

        result.solve_runtime_microseconds += runtime;

        result.routing_base = semiResult.path_selection_strategy;
        result.candidate_paths = semiResult.candidate_paths;
        result.average_paths_per_pair =
            semiResult.average_paths_per_pair;

        result.demand_evaluations.push_back({
            .demand_type = demandType,
            .congestion = semiResult.congestion,
            .runtime_microseconds = runtime
        });

        result.demand_schemes[demandModelName(demandType)] =
            std::move(semiResult.scheme);
    }

    result.total_runtime_microseconds =
        result.preprocessing_runtime_microseconds +
        result.solve_runtime_microseconds;

    result.status = ResultStatus::OK;
    return result;
}