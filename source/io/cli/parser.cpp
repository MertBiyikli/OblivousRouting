//
// Created by Mert on 30.10.26.
//

#include "io/cli/parser.h"

#include "io/cli/options.h"
#include "io/cli/validator.h"
#include "io/parse_argument_io.h"
inline Result<std::vector<SolverType>> parseSolvers(const std::string& value) {
    auto parsed = parse_solver_list(value);

    if (!parsed) {
        return makeErrorMessage(ErrorCode::InvalidSolver, "Unknown solver list: " + value);
    }
    return parsed.value();
}

inline Result<std::vector<DemandModelType>> parseDemands(const std::string& value) {
    auto parsed = parse_demand_model_list(value);
    if (!parsed) {
        return makeErrorMessage(ErrorCode::InvalidDemand, "Unknown demand list: " + value);
    }
    return parsed.value();
}

inline std::optional<OutputFormat> inferOutputFormat(const std::string& filename) {
    const auto dot = filename.find_last_of('.');

    if (dot == std::string::npos || dot + 1 >= filename.size()) {
        return std::nullopt;
    }

    return parse_output_format(filename.substr(dot + 1));
}

inline Result<GraphFormat> parseGraphFormat(const std::string& value) {
    auto parsed = parse_graph_format_token(value);
    if (!parsed) {
        return makeErrorMessage(ErrorCode::InvalidGraphFormat, "Unknown graph format: " + value);
    }
    return parsed.value();
}

inline Result<OutputFormat> parseOutputFormat(const std::string& value) {
    auto parsed = parse_output_format(value);
    if (!parsed) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Unknown output format: " + value);
    }
    return parsed.value();
}

Result<Command> parse(int argc, char** argv) {
    if (argc <= 1) {
        return Command{.action = Action::Help};
    }

    const std::string first = argv[1];
    if (first == "--help" || first == "-h" || first == "help") {
        return Command{.action = Action::Help};
    }

    if (first == "--list-solvers" || first == "list-solvers") {
        return Command{.action = Action::ListSolvers};
    }

    const bool isN1 = first == "analyze-n1";
    const bool isTraffic = first == "analyze-traffic";
    const bool isScenario = first == "analyze-scenario";

    if (first != "solve" && !isN1 && !isTraffic && !isScenario) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Unknown command: " + first);
    }

    try {
        const auto options = solveOptions();
        boost::program_options::variables_map variables;

        /*
         * Boost.Program_options treats argv[0] as the program name.
         *
         * We pass the selected subcommand as that program name,
         * which causes Boost to start parsing at the first actual
         * option (--solver, --graph, ...).
         */
        boost::program_options::store(boost::program_options::command_line_parser(argc - 1, argv + 1).options(options).run(), variables);

        boost::program_options::notify(variables);

        if (variables.count("help")) {
            return Command{.action = Action::Help};
        }

        std::optional<AnalysisKind> analysisKind;

        if (variables.count("analysis")) {
            if (!isScenario) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "--analysis requires analyze-scenario");
            }

            const auto value = variables["analysis"].as<std::string>();

            if (value == "link-failure") {
                analysisKind = AnalysisKind::LinkFailure;
            } else if (value == "capacity") {
                analysisKind = AnalysisKind::CapacityComparison;
            } else {
                return makeErrorMessage(ErrorCode::InvalidArgument, "Unsupported analysis kind: " + value);
            }
        }

        /*
         * Required parameters.
         */
        if (!variables.count("solver")) {
            return makeErrorMessage(ErrorCode::InputNotFound, "Missing required option: --solver");
        }

        if (!variables.count("graph")) {
            return makeErrorMessage(ErrorCode::InputNotFound, "Missing required option: --graph");
        }

        /*
         * Threads.
         */
        const int threads = variables["threads"].as<int>();

        if (threads <= 0) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "--threads must be a positive integer.");
        }

        /*
         * Seed.
         */
        const int seed = variables["seed"].as<int>();

        if (seed < 0) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "--seed must be a non-negative integer.");
        }

        ExecutionMode executionMode = ExecutionMode::FAST;

        if (variables.count("execution-mode")) {
            if (!isN1 && !isScenario) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "--execution-mode is currently supported only by analyze-n1 and analyze-scenario");
            }

            const auto mode = variables["execution-mode"].as<std::string>();

            if (mode == "fast") {
                executionMode = ExecutionMode::FAST;
            } else if (mode == "standard") {
                executionMode = ExecutionMode::STANDARD;
            } else if (mode == "research") {
                return makeErrorMessage(ErrorCode::InvalidArgument, "N-1 RESEARCH mode requires multi-solver support "
                                                                    "and is not implemented yet");
            } else {
                return makeErrorMessage(ErrorCode::InvalidArgument, "Unknown execution mode: " + mode);
            }
        }

        /*
         * Solver list.
         */
        auto solvers = parseSolvers(variables["solver"].as<std::string>());

        if (!solvers) {
            return getError(solvers);
        }

        /*
         * Optional demand models.
         */
        std::vector<DemandModelType> demands;

        if (variables.count("demand")) {
            auto parsed = parseDemands(variables["demand"].as<std::string>());

            if (!parsed) {
                return getError(parsed);
            }

            demands = std::move(parsed.value());
        }

        /*
         * Graph representation.
         */
        auto graphFormat = parseGraphFormat(variables["graph-format"].as<std::string>());

        if (!graphFormat) {
            return getError(graphFormat);
        }

        /*
         * Output filename.
         */
        std::string outputFilename;

        if (variables.count("output")) {

            outputFilename = variables["output"].as<std::string>();

            if (outputFilename.empty()) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "--output must not be empty.");
            }
        }

        /*
         * Output format.
         *
         * No output file -> json.
         */
        OutputFormat outputFormat = OutputFormat::JSON;

        if (variables.count("output-format")) {

            auto parsed = parseOutputFormat(variables["output-format"].as<std::string>());

            if (!parsed) {
                return getError(parsed);
            }

            outputFormat = parsed.value();

        } else if (!outputFilename.empty()) {

            if (auto inferred = inferOutputFormat(outputFilename)) {
                outputFormat = *inferred;
            } else {

                return makeErrorMessage(ErrorCode::FormatNotFound, "Could not infer output format from '" + outputFilename +
                                                                       "'. Use .json/.txt or pass "
                                                                       "--output-format.");
            }
        }

        const double demandScale = variables["demand-scale"].as<double>();

        if (!std::isfinite(demandScale) || demandScale <= 0.0) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "--demand-scale must be finite and positive");
        }

        std::vector<_CapacityUpgrade> upgrades;

        if (variables.count("capacity-upgrade")) {
            for (const auto& token : variables["capacity-upgrade"].as<std::vector<std::string>>()) {
                const auto colon = token.find(':');

                if (colon == std::string::npos) {
                    return makeErrorMessage(ErrorCode::InvalidArgument, "Expected capacity upgrade as edge_id:factor");
                }

                try {
                    std::size_t usedId = 0;
                    std::size_t usedFactor = 0;
                    const std::string idText = token.substr(0, colon);
                    const std::string factorText = token.substr(colon + 1);

                    const int edgeId = std::stoi(idText, &usedId);
                    const double factor = std::stod(factorText, &usedFactor);

                    if (usedId != idText.size() || usedFactor != factorText.size() || edgeId < 0 || !std::isfinite(factor) || factor < 1.0) {
                        return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid capacity upgrade: " + token);
                    }

                    upgrades.push_back({edgeId, factor});
                } catch (const std::exception&) {
                    return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid capacity upgrade: " + token);
                }
            }
        }

        std::vector<int> failedLinks;

        if (variables.count("link-failure")) {
            const auto& tokens = variables["link-failure"].as<std::vector<std::string>>();
            std::unordered_set<int> seen;

            for (const auto& token : tokens) {
                try {
                    std::size_t used = 0;

                    const int edgeId = std::stoi(token, &used);

                    if (used != token.size() || edgeId < 0) {
                        return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid link failure edge ID: " + token);
                    }

                    if (!seen.insert(edgeId).second) {
                        return makeErrorMessage(ErrorCode::InvalidArgument, "Duplicate link failure edge ID: " + token);
                    }

                    failedLinks.push_back(edgeId);
                } catch (const std::exception&) {
                    return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid link failure edge ID: " + token);
                }
            }
        }

        /*
                if ((demandScale != 1.0 || !upgrades.empty()) && !variables.count("load-routing-table")) {
                    return makeErrorMessage(ErrorCode::InvalidArgument, "Traffic/capacity scenarios currently require --load-routing-table");
                }*/

        /*
         * Visualization output.
         */
        std::string visualizationDirectory;

        if (variables.count("visualization")) {

            visualizationDirectory = variables["visualization"].as<std::string>();

            if (visualizationDirectory.empty()) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "--visualization must not be empty.");
            }
        }

        /*
         * Failure analysis and recovery.
         */
        bool failure_analysis = variables["failure-analysis"].as<bool>();

        bool failure_recovery = variables["failure-recovery"].as<bool>();

        if (failure_recovery && !failure_analysis) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "--failure-recovery requires --failure-analysis");
        }

        std::string routingTableInput;
        std::string routingTableDirectory;

        if (variables.count("load-routing-table")) {
            routingTableInput = variables["load-routing-table"].as<std::string>();

            if (routingTableInput.empty()) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "--load-routing-table requires a non-empty path");
            }

            if (solvers.value().size() != 1) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "Loading a routing table requires exactly one --solver");
            }

            if (failure_recovery) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "Layer-2 recomputation is not supported when loading a saved routing table");
            }

            if (!routingTableDirectory.empty()) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "Cannot use --load-routing-table and --save-routing-table together");
            }
        }

        if (variables.count("save-routing-table")) {
            routingTableDirectory = variables["save-routing-table"].as<std::string>();

            if (routingTableDirectory.empty()) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "--save-routing-table requires a non-empty directory");
            }
        }

        const CliValidationContext validation{
            .action = isN1         ? Action::AnalyzeN1
                      : isTraffic  ? Action::AnalyzeTraffic
                      : isScenario ? Action::AnalyzeScenario
                                   : Action::Run,

            .solver_count = solvers.value().size(),
            .demand_count = demands.size(),

            .valid_json_output = !outputFilename.empty() && outputFormat == OutputFormat::JSON && inferOutputFormat(outputFilename) == OutputFormat::JSON,

            .has_routing_table_input = !routingTableInput.empty(),

            .has_capacity_upgrades = !upgrades.empty(),
            .has_link_failures = !failedLinks.empty(),

            .has_load_routing_table = variables.count("load-routing-table") != 0,

            .has_save_routing_table = variables.count("save-routing-table") != 0,

            .has_visualization = variables.count("visualization") != 0,

            .has_execution_mode = variables.count("execution-mode") != 0,

            .failure_analysis = failure_analysis,
            .failure_recovery = failure_recovery,
            .has_analysis_option = variables.count("analysis") != 0,

            .has_exactly_one_failure_target = failedLinks.size() == 1,

            .has_invalid_failure_target_count = failedLinks.size() != 1,
            .analysis_kind = analysisKind,
            .capacity_upgrade_count = upgrades.size(),
            .link_failure_count = failedLinks.size(),
        };

        auto compatible = validateCommandCompatibility(validation);

        if (!compatible) {
            return getError(compatible);
        }

        std::vector<CommodityScale> commodityScales;

        if (variables.count("commodity-scale")) {
            // The solve command also supports workspace demand modifications.
            if (!isTraffic && !isN1 & !isScenario) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "--commodity-scale is supported only by "
                                                                    "analyze-traffic,analyze-n1 and analyze-scenario");
            }

            std::set<std::pair<int, int>> seen;

            for (const auto& token : variables["commodity-scale"].as<std::vector<std::string>>()) {

                const auto first = token.find(':');
                const auto second = token.find(':', first == std::string::npos ? first : first + 1);

                if (first == std::string::npos || second == std::string::npos || token.find(':', second + 1) != std::string::npos) {
                    return makeErrorMessage(ErrorCode::InvalidArgument, "Expected --commodity-scale source:target:factor: " + token);
                }

                try {
                    const auto sourceText = token.substr(0, first);
                    const auto targetText = token.substr(first + 1, second - first - 1);
                    const auto factorText = token.substr(second + 1);

                    std::size_t sourceUsed = 0;
                    std::size_t targetUsed = 0;
                    std::size_t factorUsed = 0;

                    const int source = std::stoi(sourceText, &sourceUsed);
                    const int target = std::stoi(targetText, &targetUsed);
                    const double factor = std::stod(factorText, &factorUsed);

                    if (sourceUsed != sourceText.size() || targetUsed != targetText.size() || factorUsed != factorText.size() || source < 0 || target < 0 || source == target ||
                        !std::isfinite(factor) || factor <= 0.0) {
                        return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid commodity scale: " + token);
                        }

                    if (!seen.insert({source, target}).second) {
                        return makeErrorMessage(ErrorCode::InvalidArgument, "Duplicate commodity override: " + token);
                    }

                    commodityScales.push_back({source, target, factor});

                } catch (const std::exception&) {
                    return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid commodity scale: " + token);
                }
            }
        }

        std::optional<ScenarioAnalysisRequest> scenarioRequest;

        if (isScenario) {
            ScenarioAnalysisRequest request;

            request.kind = *analysisKind;
            request.execution_mode = executionMode;

            request.workspace.demand_scale = demandScale;

            for (const auto& item : commodityScales) {
                request.workspace.commodity_overrides.push_back({item.source, item.target, item.factor});
            }

            switch (request.kind) {
            case AnalysisKind::LinkFailure: {
                // Upgrades establish the baseline workspace.
                request.workspace.capacity_upgrades = upgrades;
                request.failed_edge_id = failedLinks.front();
                break;
            }

            case AnalysisKind::CapacityComparison: {
                // Capacity changes are scenario events, not baseline changes.
                request.workspace.capacity_upgrades.clear();
                request.failed_edge_id.reset();

                request.capacity_upgrade_events.reserve(upgrades.size());

                for (const auto& upgrade : upgrades) {
                    request.capacity_upgrade_events.push_back(
                        CapacityUpgradeEvent{
                            upgrade.edge_id,
                            upgrade.factor
                        }
                    );
                }

                break;
            }
            }

            scenarioRequest = std::move(request);
        } // if (isScenario)
            /*
             * The CLI's job ends here.
             *
             * Everything below the CLI receives a normal
             * Config object and does not need to know about
             * argv / flags / subcommands.
             */
            Config config{.solvers = std::move(solvers.value()),
                          .filename = variables["graph"].as<std::string>(),
                          .evaluate_demand_models = !demands.empty(),
                          .demand_models = std::move(demands),
                          .offline_opt_per_model = {},
                          .demand_maps = {},
                          .graph_format = graphFormat.value(),
                          .num_threads = threads,
                          .output_filename = std::move(outputFilename),
                          .output_format = outputFormat,
                          .seed = seed,
                          .failure_analysis = failure_analysis,
                          .failure_recovery = failure_recovery,
                          .routing_table_input_file = std::move(routingTableInput),
                          .routing_table_output_directory = std::move(routingTableDirectory),
                          .visualization_output_directory = std::move(visualizationDirectory),
                          .demand_scale = demandScale,
                          .capacity_upgrades = std::move(upgrades),
                          .failed_links = std::move(failedLinks)};
            return Command{
                .action = isN1         ? Action::AnalyzeN1
                          : isTraffic  ? Action::AnalyzeTraffic
                          : isScenario ? Action::AnalyzeScenario
                                       : Action::Run,
                .config = std::move(config),
                .execution_mode = executionMode,
                .commodity_scales = std::move(commodityScales),
                .scenario_request = std::move(scenarioRequest),

            };
        }
        catch (const boost::program_options::error& error) {
            return makeErrorMessage(ErrorCode::InvalidArgument, error.what());
        }
        catch (const std::exception& error) {
            return fromStdException(error, ErrorCode::InvalidArgument);
        }
    }


