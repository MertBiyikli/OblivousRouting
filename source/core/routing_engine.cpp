//
// Created by Mert Biyikli on 23.06.26.
//

#include "routing/routing_engine.h"

std::optional<RoutingRunResult> RoutingEngine::solve(
        IGraph& graph,
        const Config& cfg,
        const SolverType& type) {


        if (isSemiObliviousSolver(type)) {
            return solveSemiOblivious(graph, cfg, type);
        }


        RoutingRunResult result;
        result.type = type;
        result.solver_name = getSolverName(type);

        auto solver_opt = makeSolver(type, graph);
        if (!solver_opt) {
            std::cerr << "[ERROR] Failed to create solver of type " << static_cast<int>(type) << "\n";
            result.status = ResultStatus::ERROR_INVALID_SOLVER;
        }
        auto& solver = *solver_opt;
        auto t0 = timeNow();
        result.scheme = solver->solve();
        auto t1 = timeNow();
        result.total_runtime_microseconds = duration(t1-t0);


        //result.scheme->printRoutingTable();


        // Print time statistics if available
        if (auto mwu = dynamic_cast<MWUFramework*>(solver.get())) {
            //mwu->printTimeStats();
            result.mwu_iterations = mwu->getIterationCount();
        }else {
            result.mwu_iterations = -1;
        }

        // Compute oblivious ratio for linear schemes
        if (auto linear_scheme = dynamic_cast<LinearRoutingScheme*>(result.scheme.get())) {
            result.oblivious_ratio = linear_scheme->computeObliviousRatio();
        }
        if ( auto lp = dynamic_cast<LP*>(solver.get())) {
            result.oblivious_ratio = lp->alpha->solution_value();
        }

        // Evaluate demand models if provided
        if (cfg.evaluate_demand_models) {
            auto pairs = generateAllDemandPairs(graph);
            for (const auto& type : cfg.demand_models) {

                auto model = makeDemandModel(type);
                demands dmap = model->generate(graph, pairs);

                if (!result.scheme) {
                    std::cerr << "[ERROR]: Solver returned null routing scheme\n";
                    return std::nullopt;
                }
                double scheme_congestion = computeRoutingSchemeCongestion(graph, result.scheme, dmap);
                result.demand_evaluations.push_back({.demand_type = type, .congestion = scheme_congestion});
                //printStatsForDemandModel(model_name, {offline_cong, scheme_cong});
            }
        }
        result.status = ResultStatus::OK;
        return result;
    }


    bool RoutingEngine::isSemiObliviousSolver(SolverType type) {
        return type == SolverType::SEMI_ELECTRICAL ||
               type == SolverType::SEMI_TREE;
    }



    std::optional<RoutingRunResult> RoutingEngine::solveSemiOblivious(
    IGraph& graph,
    const Config& cfg,
    SolverType type
) {
        RoutingRunResult result;
        result.type = type;
        result.solver_name = getSolverName(type);

        result.status = ResultStatus::OK;
        result.oblivious_ratio = -1.0;
        result.mwu_iterations = -1;
        result.total_runtime_microseconds = 0.0;

        auto routingEngine = makeSemiRoutingEngine(type, graph);
        auto optimizer = std::make_shared<OrToolsSemiObliviousLoadOptimizer>();

        SemiObliviousRoutingSolver solver(
            graph,
            routingEngine,
            optimizer
        );

        const auto preprocessStart = timeNow();
        auto candidateScheme = solver.preprocess();
        result.preprocessing_runtime_microseconds = duration(timeNow() - preprocessStart);

        if (!cfg.evaluate_demand_models) {
            std::cerr << "[ERROR] Semi-oblivious solver requires demand models.\n";
            result.status = ResultStatus::ERROR_INVALID_SOLVER;
            return std::nullopt;
        }

        auto pairs = generateAllDemandPairs(graph);

        result.solve_runtime_microseconds = 0;
        for (const auto& demandType : cfg.demand_models) {
            auto model = makeDemandModel(demandType);
            demands dmap = model->generate(graph, pairs);

            auto t0 = timeNow();
            auto semiResult = solver.route(dmap, demandType);
            result.solve_runtime_microseconds += duration(timeNow() - t0);
            // store the scheme for each demand

            result.demand_schemes[demandModelName(demandType)] = std::move(semiResult.scheme);

            printSemiObliviousResult(semiResult);

            result.demand_evaluations.push_back({
                .demand_type = demandType,
                .congestion = semiResult.congestion,
                .runtime_microseconds = semiResult.runtime_microseconds
                    });

        }
        result.total_runtime_microseconds =
            result.preprocessing_runtime_microseconds +
            result.solve_runtime_microseconds;

        return result;
    }