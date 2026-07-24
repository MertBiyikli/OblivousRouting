//
// Created by Mert Biyikli on 23.06.26.
//

//
// Created by Mert Biyikli on 23.06.26.
//
#include "routing/routing_engine.h"
#include "routing/routing_runner.h"
#include "routing/routing_validation.h"
#include "algorithms/oblivious/oblivious_routing_runner.h"
#include "algorithms/semi_oblivious/semi_routing_runner.h"
#include "core/errors.h"
#include "io/result_io.h"


bool isSemiObliviousSolver(SolverType type) {
    return type == SolverType::SEMI_ELECTRICAL ||
           type == SolverType::SEMI_TREE ||
           type == SolverType::SEMI_EXPANDER_HIERARCHY;
}

std::unique_ptr<IRoutingExperimentRunner> makeRunner(SolverType type) {
    if (isSemiObliviousSolver(type)) {
        return std::make_unique<SemiObliviousSolverRunner>();
    }
    return std::make_unique<ObliviousSolverRunner>();
}


Result<void> RoutingEngine::entry(int argc, char **argv) {
    auto cfg = parse_parameter(argc, argv);
    if (!cfg) {
        return getError(cfg);
    }

    auto graph = load_graph(cfg.value(), argc, argv);
    if (!graph) {
        return getError(graph);
    }
/*
    auto offline = offlineOptimal(graph.value(), cfg.value());
    if (!offline) {
        return getError(offline);
    }
    */

    RoutingEngine engine;

    for (SolverType type : cfg.value().solvers ) {
        auto result = engine.solve(*graph.value(), cfg.value(), type);

        if (!result) {
            return getError(result);
        }

        std::string out = (cfg->output_filename.empty() ? "result/run_" +getSolverName(type)+".json" : cfg->output_filename);

        auto output = RoutingResultWriter::write(result.value(), out, cfg.value().output_format);
        if (!output) {
            return getError(output);
        }

    }

    m_graph = std::move(graph.value());
    this->cfg = cfg.value();
    return {};
}



Result<IRoutingResult> RoutingEngine::solve(IGraph& graph,const Config& cfg,SolverType type) {

    RoutingValidation validator;

    auto input = validator.input(graph, cfg);
    if (!input) {
        return getError(input);
    }

    auto runner = makeRunner(type);

    auto result = runner->run(graph, cfg, type);
    if (!result) {
        return getError(result);
    }

    auto output = validator.output(result.value());
    if (!output) {
        return getError(output);
    }

    return result;
}

