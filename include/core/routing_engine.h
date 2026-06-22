//
// Created by Mert Biyikli on 09.06.26.
//

#ifndef OBLIVIOUSROUTING_ROUTING_ENGINE_H
#define OBLIVIOUSROUTING_ROUTING_ENGINE_H

#include "routing_result.h"
#include "utils.h"
#include <optional>

class RoutingEngine
{
public:
    std::optional<RoutingRunResult> solve(
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

private:
    static bool isSemiObliviousSolver(SolverType type) {
        return type == SolverType::SEMI_ELECTRICAL ||
               type == SolverType::SEMI_TREE;
    }

    static std::shared_ptr<IRoutingEngine> makeSemiRoutingEngine(
        SolverType type,
        IGraph& graph
    ) {
        switch (type) {
            case SolverType::SEMI_ELECTRICAL:
                return std::make_shared<ExistingSolverRoutingEngine>(std::make_shared<ElectricalMWU>(graph, 0, true));

            case SolverType::SEMI_TREE:
                return std::make_shared<ExistingSolverRoutingEngine>(std::make_shared<TreeMWU<FlatHST>>(graph,0, std::make_unique<FastCKR<FlatHST>>(graph)));

            default:
                throw std::invalid_argument(
                    "Requested semi-oblivious routing engine for non-semi solver"
                );
        }
    }

    std::optional<RoutingRunResult> solveSemiOblivious(
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
            routingEngine,
            optimizer
        );

        const auto preprocessStart = timeNow();
        auto candidateScheme = solver.preprocess(graph);
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
};



#endif //OBLIVIOUSROUTING_ROUTING_ENGINE_H