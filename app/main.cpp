#include "../include/io/cli.h"
#include "../include/routing/routing_engine.h"

#include "analysis/capacity_upgrade/capacity_upgrade.h"
#include "analysis/failures/link_failures.h"
#include "analysis/resilience/resilience_analyzer.h"
#include "analysis/resilience/resilience_result_json.h"
#include "analysis/scenario_engine.h"
#include "analysis/scenario_result_json.h"
#include "analysis/workspace/workspace_configuration.h"
#include "analysis/workspace/workspace_preparer.h"
#include "io/result_io.h"
#include "routing/demand_preparation.h"
#include "routing/routing_runner.h"
#include "routing/storage/routing_scheme_verifier.h"
#include "routing/storage/routing_table_adapter.h"
#include "routing/storage/routing_table_exporter.h"
#include "routing/storage/routing_table_importer.h"
#include <cstdlib>
#include <string>

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
    case Action::AnalyzeN1:
    case Action::AnalyzeTraffic:
        break;
    }

    // Running engine
    auto cfg = cmd.value().config;
    // ------------------------------------------------------------
    // N-1 network resilience analysis
    // ------------------------------------------------------------
    if (cmd->action == Action::AnalyzeN1) {

        // Load the original network using the existing graph loader.
        auto loadedGraph = load_graph_optimized(cfg);

        if (!loadedGraph) {
            std::cerr << "Error loading N-1 graph: " << loadedGraph.error().message << '\n';
            return 1;
        }

        optimized::Graph<EdgeData> graph = *loadedGraph.value();

        // Construct the complete workspace from the ORIGINAL graph.
        WorkspaceConfiguration modifications;

        modifications.capacity_upgrades = cfg.capacity_upgrades;
        modifications.demand_scale = cfg.demand_scale;

        for (const auto& item : cmd->commodity_scales) {
            modifications.commodity_overrides.push_back({item.source, item.target, item.factor});
        }

        // Generate original demands first, then apply traffic
        // transformations and capacity upgrades.
        auto workspace = WorkspacePreparer::prepare(graph, cfg, modifications);

        if (!workspace) {
            std::cerr << "Workspace preparation failed: " << getError(workspace).error().message << '\n';
            return 1;
        }

        // Compute the N-1 baseline against the effective graph
        // and the externally prepared demand matrix.
        auto prepared = ScenarioEngine::prepare(graph, cfg, cmd->execution_mode, workspace->effective_demands);

        if (!prepared) {
            std::cerr << "N-1 preparation failed: " << getError(prepared).error().message << '\n';
            return 1;
        }

        // Reuse the tested N-1 overload.
        // Every failure is restored to the effective workspace.
        auto analysis = analyzeN1(graph, prepared.value());

        if (!analysis) {
            std::cerr << "N-1 execution failed: " << getError(analysis).error().message << '\n';
            return 1;
        }

        // Verify that analysis did not modify the effective graph.
        if (optimized::GraphFingerprint::compute(graph) != workspace->effective_fingerprint) {
            std::cerr << "N-1 failed to restore workspace graph\n";
            return 1;
        }

        auto written = ResilienceResultJsonWriter::write(analysis.value(), cfg.output_filename);

        if (!written) {
            std::cerr << "N-1 analysis error: " << getError(written).error().message << '\n';
            return 1;
        }

        std::cout << "N-1 report written to: " << cfg.output_filename << '\n';

        return 0;
    }

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

        if (cmd->action == Action::AnalyzeTraffic) {
            // Build the baseline independently of the requested growth.
            Config baselineConfig = cfg;
            baselineConfig.demand_scale = 1.0;

            auto generated = ScenarioEngine::computeTotalDemandVolume(graph, baselineConfig);

            if (!generated) {
                std::cerr << "Traffic demand generation error: " << getError(generated).error().message << '\n';
                return 1;
            }

            TrafficGrowth scenario;

            auto configured = scenario.setGlobalScale(cfg.demand_scale);

            if (!configured) {
                std::cerr << "Traffic configuration error: " << getError(configured).error().message << '\n';
                return 1;
            }

            // Apply directed source-target demand multipliers.
            for (const auto& override : cmd->commodity_scales) {
                auto added = scenario.growPerCommodity(override.source, override.target, override.factor);

                if (!added) {
                    std::cerr << "Invalid commodity override (" << override.source << ", " << override.target << "): " << getError(added).error().message << '\n';

                    return 1;
                }
            }

            auto analysis = ScenarioEngine::evaluateTrafficGrowth(graph, result, generated.value().first, scenario);

            if (!analysis) {
                std::cerr << "Traffic analysis error: " << getError(analysis).error().message << '\n';
                return 1;
            }

            // Record the transformation the user actually requested.
            // These fields describe inputs, not derived congestion metrics.
            analysis->requested_global_scale = cfg.demand_scale;

            analysis->commodity_overrides.clear();
            analysis->commodity_overrides.reserve(cmd->commodity_scales.size());

            for (const auto& item : cmd->commodity_scales) {
                analysis->commodity_overrides.push_back({item.source, item.target, item.factor});
            }

            // Canonicalize output order for reproducible reports.
            std::sort(analysis->commodity_overrides.begin(), analysis->commodity_overrides.end(),
                      [](const auto& a, const auto& b) { return std::tie(a.source, a.target) < std::tie(b.source, b.target); });

            auto written = TrafficGrowthResultJsonWriter::write(analysis.value(), cfg.output_filename);

            if (!written) {
                std::cerr << "Traffic report error: " << getError(written).error().message << '\n';
                return 1;
            }

            std::cout << "Traffic analysis report written to: " << cfg.output_filename << '\n';

            std::cout << "Baseline max congestion: " << analysis->baseline.max_congestion << '\n';

            std::cout << "Modified max congestion: " << analysis->modified.max_congestion << '\n';

            return 0;
        }

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

        std::cout << "Routing scheme loaded and verified. Optimization solver was not executed.\n";
        return 0;
    }

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
    // Persist newly computed routing schemes when requested.
    // ------------------------------------------------------------
    if (!cfg.routing_table_output_directory.empty()) {
        auto exported = RoutingTableExporter::save(res.value(), graph, cfg.routing_table_output_directory);

        if (!exported) {
            std::cerr << "Routing table export error: " << getError(exported).error().message << '\n';
            return 1;
        }

        std::cout << "Routing tables exported to: " << cfg.routing_table_output_directory << '\n';
    }

    // Obtain the graph when exporting routing tables
    // or executing a scenario.
    const bool hasCapacityUpgrades = !cfg.capacity_upgrades.empty();

    const bool hasLinkFailures = !cfg.failed_links.empty();

    const bool hasScenario = hasCapacityUpgrades || hasLinkFailures;

    // Execute requested topology/capacity scenarios.
    if (hasScenario) {

        ScenarioEngine scenario_engine;

        ExecutionMode scenario_mode = cmd->execution_mode;

        if (const char* value = std::getenv("ER_SCENARIO_EXECUTION_MODE")) {

            const std::string mode(value);

            if (mode == "standard") {
                scenario_mode = ExecutionMode::STANDARD;
            } else if (mode == "research") {
                scenario_mode = ExecutionMode::RESEARCH;
            }
        }

        // Configure, validate, execute and optionally
        // persist one individual scenario.
        auto executeScenario = [&](Scenario& scenario) -> Result<void> {
            auto configured = scenario.configure(cfg);

            if (!configured)
                return getError(configured);

            auto verified = scenario.verify(graph);

            if (!verified)
                return getError(verified);

            auto result = scenario_engine.runScenario(graph, scenario, cfg, scenario_mode);

            if (!result)
                return getError(result);

            std::filesystem::path reportPath;

            const bool isCapacity = scenario.kind() == ScenarioKind::CAPACITY_UPGRADE;

            // Default paths for the individual scenario types.
            if (isCapacity) {
                reportPath = "results/scenarios/capacity.json";
            } else {
                reportPath = "results/scenarios/failure.json";
            }

            // Allow callers and integration tests to override the path.
            if (const char* configuredPath = std::getenv("ER_SCENARIO_REPORT_PATH")) {

                if (*configuredPath != '\0') {
                    reportPath = configuredPath;

                    // Both scenarios may execute in one CLI invocation.
                    // Give each scenario its own output file.
                    if (hasCapacityUpgrades && hasLinkFailures) {
                        const std::string suffix = isCapacity ? "_capacity" : "_failure";

                        const auto parent = reportPath.parent_path();
                        const auto stem = reportPath.stem().string();
                        const auto extension = reportPath.extension().string();

                        reportPath = parent / (stem + suffix + (extension.empty() ? ".json" : extension));
                    }
                }
            }

            auto written = ScenarioResultJsonWriter::write(result.value(), reportPath);

            if (!written) {
                return getError(written);
            }

            std::cout << "Scenario JSON written to: " << reportPath << '\n';

            return {};
        };

        // Capacity upgrades.
        if (hasCapacityUpgrades) {
            CapacityUpgrade scenario;

            auto result = executeScenario(scenario);

            if (!result) {
                std::cerr << "Capacity scenario error: " << getError(result).error().message << '\n';

                return 1;
            }
        }

        // Link failures.
        if (hasLinkFailures) {
            LinkFailure scenario;

            auto result = executeScenario(scenario);

            if (!result) {
                std::cerr << "Link failure scenario error: " << getError(result).error().message << '\n';

                return 1;
            }
        }
    }

    return 0;
}
