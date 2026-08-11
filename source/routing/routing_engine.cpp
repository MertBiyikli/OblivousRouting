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
#include "visualization/json_exporter.h"


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

    auto graph = load_graph_optimized(cfg.value(), argc, argv);
    if (!graph) {
        return getError(graph);
    }

    RoutingEngine engine;

    for (SolverType type: cfg.value().solvers) {
        auto result = engine.solve(*graph.value(), cfg.value(), type);

        if (!result) {
            return getError(result);
        }

        std::string out = (cfg->output_filename.empty() ? "result/run_" + getSolverName(type) + ".json" : cfg->output_filename);

        auto output = RoutingResultWriter::write(result.value(), out, cfg.value().output_format);
        if (!output) {
            return getError(output);
        }

        if (!cfg->visualization_output_directory.empty()) {
            for (const auto &visualization: result->visualization_results) {
                const std::string filename = safeFileComponent(visualization.solver_name) + "__" + safeFileComponent(visualization.demand_model) + ".json";

                const auto visualization_path = std::filesystem::path(cfg->visualization_output_directory) / filename;

                auto visualization_output = RoutingVisualizationJsonExporter::write(visualization, visualization_path);

                if (!visualization_output) {
                    return getError(visualization_output);
                }
            }
        }
    }

    m_graph = std::move(graph.value());
    this->cfg = cfg.value();
    return {};
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
