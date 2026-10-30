//
// Created by Mert Biyikli on 23.06.26.
//
#include "core/errors.h"
#include "io/result_io.h"
#include "routing/routing_engine.h"
#include "routing/routing_runner.h"
#include "routing/routing_validation.h"


std::string safeFileComponent(std::string value) {
    for (char &character: value) {
        const auto byte =
                static_cast<unsigned char>(character);

        if (
            !std::isalnum(byte) &&
            character != '-' &&
            character != '_'
        ) {
            character = '_';
        }
    }

    return value;
}




Result<RoutingExperimentResult> RoutingEngine::entry(Config &cfg) {

    auto graph = load_graph_optimized(cfg);
    if (!graph) {
        return getError(graph);
    }

    RoutingEngine engine;

    RoutingExperimentResult exp_result;
    exp_result.graph_path = cfg.filename;

    for (SolverType type: cfg.solvers) {
        auto result = engine.solve(*graph.value(), cfg, type);

        if (!result) {
            return getError(result);
        }

        if (exp_result.graph_name.empty()) {
            exp_result.graph_name = std::filesystem::path(result->graph_name).stem().string();
            exp_result.nodes = result->nodes;
            exp_result.edges = result->edges;
        }


        exp_result.solver_results.push_back(std::move(result.value()));
    }


    auto output = RoutingResultWriter::write(exp_result, cfg);
    if (!output) {
        return getError(output);
    }

    m_graph = std::move(graph.value());
    this->cfg = cfg;
    return exp_result;
}


Result<IRoutingResult> RoutingEngine::solve(optimized::Graph<EdgeData> &graph, const Config &cfg, SolverType type) {
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
