//
// Created by Mert on 02.11.26.
//

#include "application/routing_application.h"

#include "io/parse_argument_io.h"
#include "io/result_io.h"

#include "routing/demand_preparation.h"
#include "routing/routing_engine.h"
#include "routing/routing_runner.h"

#include "routing/storage/routing_scheme_verifier.h"
#include "routing/storage/routing_table_exporter.h"
#include "routing/storage/routing_table_importer.h"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <utility>

int RoutingApplication::run(const Command& command) {
    if (command.action != Action::Run) {
        std::cerr << "Invalid command for RoutingApplication.\n";
        return 1;
    }

    auto cfg = command.config;

    // ------------------------------------------------------------
    // 1. Execute using an imported routing scheme.
    // ------------------------------------------------------------
    if (!cfg.routing_table_input_file.empty()) {
        auto loadedGraph = load_graph_optimized(cfg);

        if (!loadedGraph) {
            std::cerr << "Error loading graph: " << loadedGraph.error().message << '\n';
            return 1;
        }

        optimized::Graph<EdgeData> graph = *loadedGraph.value();

        const auto import_start = std::chrono::steady_clock::now();

        auto imported = RoutingTableImporter::load(graph, cfg.routing_table_input_file);

        const auto import_end = std::chrono::steady_clock::now();

        std::cerr << "Routing table import: " << std::chrono::duration<double>(import_end - import_start).count() << " s\n";

        if (!imported) {
            std::cerr << "Error loading routing table: " << imported.error().message << '\n';
            return 1;
        }

        const auto verify_start = std::chrono::steady_clock::now();

        auto verified = RoutingSchemeVerifier::verify(graph, *imported.value());

        const auto verify_end = std::chrono::steady_clock::now();

        std::cerr << "Routing table verification: " << std::chrono::duration<double>(verify_end - verify_start).count() << " s\n";

        if (!verified) {
            std::cerr << "Error validating routing scheme: " << verified.error().message << '\n';
            return 1;
        }

        IRoutingResult result{};

        result.status = ResultStatus::OK;
        result.type = cfg.solvers.front();
        result.solver_name = getSolverTypeName(result.type);
        result.graph_name = std::filesystem::path(cfg.filename).stem().string();

        result.nodes = graph.getNumNodes();
        result.edges = graph.getNumUndirectedEdges();

        result.total_runtime_microseconds = 0.0;
        result.solve_runtime_microseconds = 0.0;
        result.oblivious_ratio = -1.0;

        result.scheme = std::move(imported.value());

        // Preserve optional demand-model evaluations.
        if (cfg.evaluate_demand_models) {
            for (const auto demandType : cfg.demand_models) {
                auto prepared = prepareDemands(graph, demandType, cfg.demand_scale);

                if (!prepared) {
                    std::cerr << "Demand preparation error: " << prepared.error().message << '\n';
                    return 1;
                }

                const auto start = std::chrono::steady_clock::now();

                auto evaluated = DemandEvaluator::evaluate(graph, result, prepared.value());

                if (!evaluated) {
                    std::cerr << "Demand evaluation error: " << evaluated.error().message << '\n';
                    return 1;
                }

                const auto end = std::chrono::steady_clock::now();

                DemandEvaluationResult summary{};

                summary.demand_type = demandType;
                summary.congestion = evaluated->max_congestion;

                summary.runtime_microseconds = std::chrono::duration<double, std::micro>(end - start).count();

                result.demand_evaluations.push_back(std::move(summary));
            }
        }

        RoutingExperimentResult experiment{};

        experiment.graph_name = result.graph_name;
        experiment.graph_path = cfg.filename;
        experiment.nodes = result.nodes;
        experiment.edges = result.edges;

        experiment.solver_results.push_back(std::move(result));

        auto written = RoutingResultWriter::write(experiment, cfg);

        if (!written) {
            std::cerr << "Error writing results: " << written.error().message << '\n';
            return 1;
        }

        std::cout << "Routing scheme loaded and verified. "
                  << "Optimization solver was not executed.\n";

        return 0;
    }

    // ------------------------------------------------------------
    // 2. Compute a new routing scheme.
    // ------------------------------------------------------------
    RoutingEngine engine;

    auto res = engine.entry(cfg);

    if (!res) {
        std::cerr << "Error: " << getError(res).error().message << '\n';
        return 1;
    }

    if (!engine.m_graph) {
        std::cerr << "Error: RoutingEngine returned no graph.\n";
        return 1;
    }

    optimized::Graph<EdgeData>& graph = *engine.m_graph;

    // ------------------------------------------------------------
    // 3. Optionally persist newly computed routing schemes.
    // ------------------------------------------------------------
    if (!cfg.routing_table_output_directory.empty()) {
        auto exported = RoutingTableExporter::save(res.value(), graph, cfg.routing_table_output_directory);

        if (!exported) {
            std::cerr << "Routing table export error: " << getError(exported).error().message << '\n';
            return 1;
        }

        std::cout << "Routing tables exported to: " << cfg.routing_table_output_directory << '\n';
    }
    
    // Persist routing experiment results.
    auto written = RoutingResultWriter::write(res.value(), cfg);

    if (!written) {
        std::cerr << "Error writing routing results: "
                  << getError(written).error().message << '\n';
        return 1;
    }

    return 0;
}
