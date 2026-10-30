//
// Created by Mert Biyikli on 24.06.26.
//

#include "algorithms/semi_oblivious/semi_routing_runner.h"
#include "algorithms/oblivious/oblivious_routing_runner.h"
#include "analysis/failures/link_exposure_result_adapter.h"
#include "core/errors.h"

static Result<std::shared_ptr<SemiSolverRoutingEngine>> makeSemiRoutingEngine(SolverType type, optimized::Graph<EdgeData>& graph) {
    switch (type) {
    case SolverType::SEMI_ELECTRICAL:
        return std::make_shared<SemiSolverRoutingEngine>(std::make_shared<ElectricalMWU>(graph, 0, true));

    case SolverType::SEMI_TREE:
        return std::make_shared<SemiSolverRoutingEngine>(std::make_shared<TreeMWU<FlatHST>>(graph, 0, std::make_unique<FastCKR<FlatHST>>(graph)));

    case SolverType::SEMI_EXPANDER_HIERARCHY:
        return std::make_shared<SemiSolverRoutingEngine>(std::make_shared<ElectrifiedExpanderHierarchySolver>(graph, 0));

    default:
        return makeErrorMessage(ErrorCode::InvalidSolver, "Requested semi-oblivious routing engine for non-semi solver");
    }
}

Result<IRoutingResult> SemiObliviousSolverRunner::run(optimized::Graph<EdgeData>& graph, const Config& cfg, SolverType type) const {
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

    SemiObliviousRoutingSolver solver(graph, routingEngine, optimizer);

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


        SemiObliviousRoutingResult result;

        auto routed = solver.route(dmap.value(), demandType);

        if (!routed) {
            return getError(routed);
        }

        result = std::move(routed.value());

        /*
         * Demand-specific optimization runtime only.
         */
        semiResult.solve_runtime_microseconds += result.total_runtime_microseconds;

        semiResult.routing_base = result.path_selection_strategy;

        semiResult.candidate_paths = result.candidate_paths;

        semiResult.average_paths_per_pair = result.average_paths_per_pair;

        /*
         * ------------------------------------------------------------
         * Pull analysis metrics from the exact demand-specific scheme.
         * ------------------------------------------------------------
         */

        // Evaluate the exact demand-specific scheme.
        // Demand evaluation does not perform optimization or N-1.
        const auto evaluation_start = std::chrono::steady_clock::now();

        auto evaluated = DemandEvaluator::evaluate(graph, result, dmap.value());

        if (!evaluated) {
            return getError(evaluated);
        }

        const auto evaluation_end = std::chrono::steady_clock::now();

        // Preserve normal routing-evaluation metrics.
        DemandEvaluationResult evaluation{};
        evaluation.demand_type = demandType;
        evaluation.congestion = evaluated->max_congestion;
        evaluation.runtime_microseconds = std::chrono::duration<double, std::micro>(evaluation_end - evaluation_start).count();
        evaluation.network_data = std::move(evaluated.value());
        if (cfg.failure_analysis) {
            const auto failure_start = std::chrono::steady_clock::now();

            auto exposure = LinkExposureAnalyzer::analyze(
                graph,
                *result.scheme,
                dmap.value()
            );

            if (!exposure)
                return getError(exposure);

            const auto failure_end = std::chrono::steady_clock::now();

            const double runtime =
                std::chrono::duration<double, std::micro>(
                    failure_end - failure_start
                ).count();

            populateLinkExposureMetrics(
                evaluation,
                exposure.value(),
                runtime
            );
        }
        semiResult.demand_evaluations.push_back(std::move(evaluation));

        // Preserve the exact demand-specific routing scheme.
        semiResult.demand_schemes[demandModelName(demandType)] = std::move(result.scheme);

    }

    semiResult.total_runtime_microseconds = semiResult.preprocessing_runtime_microseconds + semiResult.solve_runtime_microseconds;
    semiResult.oblivious_ratio = -1;

    ObliviousSolverRunner::appendMetricsIfAvailable(routingEngine->solver_, semiResult);
    // ObliviousSolverRunner::appendObjectiveIfAvailable(routingEngine, *semiResult.scheme, semiResult);

    semiResult.status = ResultStatus::OK;
    return semiResult;
}