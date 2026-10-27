#include "../include/io/cli.h"
#include "../include/routing/routing_engine.h"

#include "analysis/scenario_engine.h"
#include "analysis/scenario_result_json.h"
#include "io/result_io.h"
#include "routing/routing_runner.h"
#include "routing/storage/routing_scheme_verifier.h"
#include "routing/storage/routing_table_adapter.h"
#include "routing/storage/routing_table_exporter.h"
#include "routing/storage/routing_table_importer.h"

int main(int argc, char** argv) {

    // Parse command line arguments
    auto cmd = parse(argc, argv);
    if (!cmd) {
        std::cerr << "Error: " << getError(cmd).error().message << std::endl;
        return 1;
    }

    switch (cmd->action) {

    case Action::Help:
        std::cout << help(argv[0]);
        return 0;

    case Action::ListSolvers:
        std::cout << solverList();
        return 0;

    case Action::Run:
        break;
    }

    // Running engine
    auto cfg = cmd.value().config;

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

        if (cfg.evaluate_demand_models) {

            auto evaluation = DemandEvaluator::evaluate(graph, result.scheme, cfg, result, cfg.failure_analysis, cfg.demand_scale);

            if (!evaluation) {
                std::cerr << "Error evaluating imported scheme: " << evaluation.error().message << '\n';
                return 1;
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

        std::cout << "Routing scheme loaded and verified. Optimization solver was not executed.\n";
    }

    RoutingEngine engine;
    auto res = engine.entry(cfg);

    if (!res) {
        std::cerr << "Error: " << getError(res).error().message << '\n';
        return 1;
    }

    // Obtain the graph when exporting routing tables
    // or executing a scenario.
    if (!cfg.routing_table_output_directory.empty() || !cfg.capacity_upgrades.empty()) {

        auto graphResult = engine.getGraph();

        if (!graphResult) {
            std::cerr << "Error: " << getError(graphResult).error().message << '\n';
            return 1;
        }

        auto& graph = *graphResult.value();

        // Routing-table export is independent of scenario execution.
        if (!cfg.routing_table_output_directory.empty()) {

            auto saved = RoutingTableExporter::save(res.value(), graph, cfg.routing_table_output_directory);

            if (!saved) {
                std::cerr << "Error: " << getError(saved).error().message << '\n';
                return 1;
            }

            std::cout << "Routing tables saved to: " << cfg.routing_table_output_directory << '\n';
        }

        // Execute capacity-upgrade scenarios independently.
        if (!cfg.capacity_upgrades.empty()) {

            ScenarioEngine scenario_engine;
            CapacityUpgrade capacity_upgrade(graph);

            auto scenario_results = scenario_engine.runScenario(graph, capacity_upgrade, cfg, ExecutionMode::STANDARD);

            if (!scenario_results) {
                std::cerr << "Scenario error: " << getError(scenario_results).error().message << '\n';
                return 1;
            }

            auto printed = printScenarioAnalysis(scenario_results.value());

            if (!printed) {
                std::cerr << "Error: " << getError(printed).error().message << '\n';
                return 1;
            }

            // After successful scenario execution:
            auto saved = ScenarioResultJsonWriter::write(scenario_results.value(), "results/scenarios/backbone_1221_capacity_upgrade.json");

            if (!saved) {
                std::cerr << "Scenario JSON export failed: " << getError(saved).error().message << '\n';
                return 1;
            }
        }
    }

    return 0;
}
