#include "../include/io/cli.h"
#include "../include/routing/routing_engine.h"
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

        auto& graph = *loadedGraph.value();

        auto imported = RoutingTableImporter::load(graph, cfg.routing_table_input_file);

        if (!imported) {
            std::cerr << "Error loading routing table: " << imported.error().message << '\n';
            return 1;
        }

        auto verified = RoutingSchemeVerifier::verify(graph, *imported.value());

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


        std::unordered_set<int> upgradedEdges;

        for (const auto& upgrade : cfg.capacity_upgrades) {
            const int edge = upgrade.edge_id;

            if (edge < 0 || edge >= graph.getNumDirectedEdges()) {
                std::cerr << "Invalid capacity-upgrade edge ID: " << edge << '\n';
                return 1;
            }

            const int reverse = graph.reverse(edge).id;
            const int canonical = std::min(edge, reverse);

            if (edge != canonical || !upgradedEdges.insert(canonical).second) {
                std::cerr << "Capacity upgrades must use distinct canonical physical edge IDs\n";
                return 1;
            }

            const auto oldData = graph.edgeData(edge);
            const double newCapacity = oldData.capacity * upgrade.factor;

            if (!std::isfinite(newCapacity) || newCapacity <= 0.0) {
                std::cerr << "Invalid upgraded capacity for edge " << edge << '\n';
                return 1;
            }

            auto newData = oldData;
            newData.capacity = newCapacity;
            graph.setEdgeData(edge, newData);
        }


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
        return 0;
    }

    RoutingEngine engine;
    auto res = engine.entry(cfg);
    if (!res) {
        std::cerr << "Error: " << getError(res).error().message << std::endl;
        return 1;
    }

    if (!cfg.routing_table_output_directory.empty()) {
        auto graphResult = engine.getGraph();

        if (!graphResult) {
            std::cerr << "Error: " << getError(graphResult).error().message << '\n';
            return 1;
        }

        auto saved = RoutingTableExporter::save(res.value(), *graphResult.value(), cfg.routing_table_output_directory);

        if (!saved) {
            std::cerr << "Error: " << getError(saved).error().message << '\n';
            return 1;
        }

        std::cout << "Routing tables saved to: " << cfg.routing_table_output_directory << '\n';
    }

    return 0; // something new
}
