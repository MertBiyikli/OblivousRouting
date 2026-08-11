//
// Created by Mert Biyikli on 24.06.26.
//

#include "algorithms/semi_oblivious/semi_routing_runner.h"
#include "algorithms/oblivious/oblivious_routing_runner.h"
#include "core/errors.h"

static Result<std::shared_ptr<SemiSolverRoutingEngine>> makeSemiRoutingEngine(SolverType type, optimized::Graph<EdgeData>& graph) {
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

        case SolverType::SEMI_EXPANDER_HIERARCHY:
            return std::make_shared<SemiSolverRoutingEngine>(
                std::make_shared<ElectrifiedExpanderHierarchySolver>(graph, 0)
            );

        default:
            return makeErrorMessage(ErrorCode::InvalidSolver, "Requested semi-oblivious routing engine for non-semi solver");
    }
}

Result<IRoutingResult> SemiObliviousSolverRunner::run(optimized::Graph<EdgeData>& graph,const Config& cfg,SolverType type) const {
    SemiObliviousRoutingResult semiResult;
    semiResult.type = type;
    semiResult.graph_name = cfg.filename;
    semiResult.nodes = graph.getNumNodes();
    semiResult.edges = graph.getNumUndirectedEdges();
    semiResult.solver_name = getSolverName(type);

    if (!cfg.evaluate_demand_models || cfg.demand_models.empty()) {
        semiResult.status = ResultStatus::ERROR_MISSING_DEMAND_MODELS;
        return makeErrorMessage(ErrorCode::InvalidDemand, "Semi-oblivious solver requires demand models.");
    }

    auto engine_factory = makeSemiRoutingEngine(type, graph);
    auto optimizer_factory = std::make_shared<OrToolsSemiObliviousLoadOptimizer>();

    if (!optimizer_factory || !engine_factory) {
        return getError(engine_factory);
    }

    std::shared_ptr<SemiSolverRoutingEngine> routingEngine = engine_factory.value();
    std::shared_ptr<OrToolsSemiObliviousLoadOptimizer> optimizer = optimizer_factory;

    SemiObliviousRoutingSolver solver(
        graph,
        routingEngine,
        optimizer
    );

    const auto preprocessStart = timeNow();

    auto pre = solver.preprocess();
    if (!pre) {
        return getError(pre);
    }
    semiResult.preprocessing_runtime_microseconds = duration(timeNow() - preprocessStart);


    auto pairs = generateAllDemandPairs(graph);

    for (const auto& demandType : cfg.demand_models) {
        auto model = makeDemandModel(demandType);
        auto dmap = model->generate(graph, pairs);
        if (!dmap) {
            return getError(dmap);
        }

        const auto t0 = timeNow();

        SemiObliviousRoutingResult result;
        if (auto res = solver.route(dmap.value(), demandType)) {
            result = std::move(res.value());
        }else {
            return getError(res);
        }

        const double runtime = duration(timeNow() - t0);

        semiResult.solve_runtime_microseconds += runtime;

        semiResult.routing_base = result.path_selection_strategy;
        semiResult.candidate_paths = result.candidate_paths;
        semiResult.average_paths_per_pair = result.average_paths_per_pair;

        semiResult.demand_evaluations.push_back({
            .demand_type = demandType,
            .congestion = result.congestion,
            .runtime_microseconds = runtime
        });

        semiResult.demand_schemes[demandModelName(demandType)] =
            std::move(result.scheme);

        RoutingVisualizationResult visualization = result.visualization;

        visualization.graph_name =
            cfg.filename;

        visualization.solver_name =
            getSolverName(type);

        semiResult.visualization_results.push_back(
            std::move(visualization)
        );
    }

    semiResult.total_runtime_microseconds =
        semiResult.preprocessing_runtime_microseconds +
        semiResult.solve_runtime_microseconds;
    semiResult.oblivious_ratio = -1;

    ObliviousSolverRunner::appendMetricsIfAvailable(routingEngine->solver_, semiResult);
    //ObliviousSolverRunner::appendObjectiveIfAvailable(routingEngine, *semiResult.scheme, semiResult);

    semiResult.status = ResultStatus::OK;
    return semiResult;
}