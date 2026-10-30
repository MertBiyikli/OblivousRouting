//
// Created by Mert on 05.10.26.
//

#ifndef E_ROUTING_CLI_H
#define E_ROUTING_CLI_H

#include "analysis/scenario_results.h"
#include "parse_argument_io.h"

#include "core/config.h"
#include <boost/program_options.hpp>
#include <unordered_set>

enum class Action {
    Run,
    AnalyzeN1,
    AnalyzeTraffic,
    Help,
    ListSolvers
};

struct CommodityScale {
    int source;
    int target;
    double factor;
};

struct Command {
    Action action = Action::Run;
    Config config{};
    ExecutionMode execution_mode = ExecutionMode::FAST;
    std::vector<CommodityScale> commodity_scales;
};


inline Result<std::vector<SolverType>> parseSolvers(const std::string& value) {
    auto parsed = parse_solver_list(value);

    if (!parsed) {
        return makeErrorMessage(ErrorCode::InvalidSolver,"Unknown solver list: " + value);
    }
    return parsed.value();

}


inline Result<std::vector<DemandModelType>> parseDemands(const std::string& value) {
    auto parsed = parse_demand_model_list(value);
    if (!parsed) {
        return makeErrorMessage(ErrorCode::InvalidDemand,"Unknown demand list: " + value);
    }
    return parsed.value();
}

inline std::optional<OutputFormat> inferOutputFormat(const std::string& filename) {
    const auto dot = filename.find_last_of('.');

    if (dot == std::string::npos ||dot + 1 >= filename.size()) {
        return std::nullopt;
    }

    return parse_output_format(filename.substr(dot + 1));
}

inline Result<GraphFormat> parseGraphFormat(const std::string& value) {
    auto parsed = parse_graph_format_token(value);
    if (!parsed) {
        return makeErrorMessage(ErrorCode::InvalidGraphFormat,"Unknown graph format: " + value);
    }
    return parsed.value();
}

inline Result<OutputFormat> parseOutputFormat(const std::string& value) {
    auto parsed = parse_output_format(value);
    if (!parsed) {
        return makeErrorMessage(ErrorCode::InvalidArgument,"Unknown output format: " + value);
    }
    return parsed.value();
}

boost::program_options::options_description solveOptions() {
    boost::program_options::options_description options("Solve options");

    options.add_options()

        ("help,h",
         "Show solve help")

        ("solver,s",
         boost::program_options::value<std::string>(),
         "Solver name or comma-separated solver list")

        ("graph,g",
         boost::program_options::value<std::string>(),
         "Input graph file")

        ("load-routing-table", boost::program_options::value<std::string>(),
            "Load a serialized routing scheme instead of executing the solver")

        ("demand,d",
         boost::program_options::value<std::string>(),
         "Demand models: uniform,gravity,bimodal,gaussian")

        ("graph-format",
         boost::program_options::value<std::string>()->default_value("csr"),
         "Graph representation: csr or adjlist")

        ("threads,j",
         boost::program_options::value<int>()->default_value(1),
         "Positive worker-thread count")

        ("output,o",
         boost::program_options::value<std::string>(),
         "Output file")

        ("output-format",
         boost::program_options::value<std::string>(),
         "Output format: json, txt, or cout")

        ("seed",
         boost::program_options::value<int>()->default_value(42),
         "Non-negative random seed")

        ("execution-mode",
        boost::program_options::value<std::string>(),
        "N-1 execution mode: fast, standard, research")


        ("failure-analysis",
         boost::program_options::bool_switch()->default_value(false),
         "Run Layer-1 static N-1 link-failure exposure analysis")

        ("failure-recovery",
         boost::program_options::bool_switch()->default_value(false),
         "Run Layer-2 N-1 routing recomputation and recovery analysis")


        ("save-routing-table",
         boost::program_options::value<std::string>(),
         "Save computed routing tables as binary files in the specified directory")



        ("demand-scale", boost::program_options::value<double>()->default_value(1.0), "Demand multiplier for persisted routing analysis")
        ("commodity-scale",boost::program_options::value<std::vector<std::string>>()->multitoken()->composing(),"Directed commodity multiplier: source:target:factor")
        ("capacity-upgrade", boost::program_options::value<std::vector<std::string>>()->multitoken()->composing(), "Physical link capacity multiplier as edge_id:factor")
        ("link-failure",boost::program_options::value<std::vector<std::string>>()->multitoken()->composing(),"Simulate physical link failures using edge IDs, ""e.g. --link-failure 14 20 24")


        ("visualization",
         boost::program_options::value<std::string>(),
         "Visualization JSON output directory");

    return options;
}

inline std::string solverList() {
    std::ostringstream out;

    out
        << "Available solvers:\n"
        << "\n"
        << "  electrical       "
        << "Electrical Flow (sketching)\n"

        << "  electrical_naive "
        << "Electrical Flow (naive/exact loads)\n"

        << "  raecke_frt       "
        << "Raecke MWU + FRT\n"

        << "  raecke_ckr       "
        << "Raecke MWU + Fast-CKR\n"

        << "  raecke_mst       "
        << "Raecke MWU + randomized MST\n"

        << "  cohen            "
        << "Applegate-Cohen LP\n"

        << "  expander         "
        << "Electrified Expander Hierarchy\n"

        << "  expander_mwu     "
        << "Expander hierarchy + MWU\n"

        << "  semi_electrical  "
        << "Semi-oblivious routing, electrical base\n"

        << "  semi_tree        "
        << "Semi-oblivious routing, tree base\n"

        << "  semi_expander    "
        << "Semi-oblivious routing, expander base\n";

    return out.str();
}

std::string help(std::string_view executable) {
    std::ostringstream out;

    const auto options = solveOptions();

    out
        << "E-Routing: Network Resilience and Capacity Planning Engine.\n"
        << "\n"

        << "Usage:\n"

        << "  "
        << executable
        << " solve "
        << "--solver <name[,name...]> "
        << "--graph <file> "
        << "[options]\n"

        << "  "
        << executable
        << " analyze-n1 --solver <name> --graph <file> --demand <model> --output <report.json> [options]\n"

        << "  "
        << executable
        << " --list-solvers\n"

        << "  "
        << executable
        << " --help\n"

        << "\n"

        << options

        << "\n"

        << "Examples:\n"

        << "  "
        << executable
        << " solve "
        << "--solver electrical "
        << "--graph "
        << "experiments/datasets/small/Backbone/1221.lgf"
        << "\n"

        << "\n"

        << "  "
        << executable
        << " solve "
        << "-s electrical,raecke_ckr "
        << "-g experiments/datasets/small/Backbone/1221.lgf "
        << "-d uniform "
        << "-j 4 "
        << "-o result/backbone.json"
        << "\n";

    return out.str();
}

inline Result<Command> parse(int argc, char** argv)
{
    if (argc <= 1) {
        return Command{ .action = Action::Help};
    }

    const std::string first = argv[1];
    if (first == "--help" ||first == "-h" ||first == "help") {
        return Command{.action = Action::Help};
    }

    if (first == "--list-solvers" || first == "list-solvers") {
        return Command{.action = Action::ListSolvers};
    }

    const bool isN1 = first == "analyze-n1";
    const bool isTraffic = first == "analyze-traffic";

    if (first != "solve" && !isN1 && !isTraffic) {
        return makeErrorMessage(
            ErrorCode::InvalidArgument,
            "Unknown command: " + first
        );
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
        boost::program_options::store(
            boost::program_options::command_line_parser(
                argc - 1,
                argv + 1
            )
            .options(options)
            .run(),
            variables
        );

        boost::program_options::notify(variables);

        if (variables.count("help")) {
            return Command{.action = Action::Help};
        }


        /*
         * Required parameters.
         */
        if (!variables.count("solver")) {
            return makeErrorMessage(ErrorCode::InputNotFound,"Missing required option: --solver");
        }


        if (!variables.count("graph")) {
            return makeErrorMessage(ErrorCode::InputNotFound,"Missing required option: --graph");
        }


        /*
         * Threads.
         */
        const int threads = variables["threads"].as<int>();

        if (threads <= 0) {
            return makeErrorMessage(ErrorCode::InvalidArgument,"--threads must be a positive integer.");
        }


        /*
         * Seed.
         */
        const int seed = variables["seed"].as<int>();

        if (seed < 0) {
            return makeErrorMessage(ErrorCode::InvalidArgument,"--seed must be a non-negative integer.");
        }

        ExecutionMode executionMode = ExecutionMode::FAST;

        if (variables.count("execution-mode")) {
            if (!isN1) {
                return makeErrorMessage(
                    ErrorCode::InvalidArgument,
                    "--execution-mode is currently supported only by analyze-n1"
                );
            }

            const auto mode = variables["execution-mode"].as<std::string>();

            if (mode == "fast") {
                executionMode = ExecutionMode::FAST;
            } else if (mode == "standard") {
                executionMode = ExecutionMode::STANDARD;
            } else if (mode == "research") {
                return makeErrorMessage(
                    ErrorCode::InvalidArgument,
                    "N-1 RESEARCH mode requires multi-solver support "
                    "and is not implemented yet"
                );
            } else {
                return makeErrorMessage(
                    ErrorCode::InvalidArgument,
                    "Unknown execution mode: " + mode
                );
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

        if (variables.count("demand"))
        {
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
                return makeErrorMessage(ErrorCode::InvalidArgument,"--output must not be empty.");
            }
        }


        /*
         * Output format.
         *
         * No output file -> json.
         */
        OutputFormat outputFormat = OutputFormat::JSON;


        if (variables.count("output-format")) {

            auto parsed =parseOutputFormat(variables["output-format"].as<std::string>());

            if (!parsed) {
                return getError(parsed);
            }

            outputFormat =parsed.value();

        } else if (!outputFilename.empty()) {

            if (auto inferred = inferOutputFormat(outputFilename)) {
                outputFormat = *inferred;
            } else {

                return makeErrorMessage(ErrorCode::FormatNotFound,
                    "Could not infer output format from '" +
                    outputFilename +
                    "'. Use .json/.txt or pass "
                    "--output-format."
                );
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

                    if (edgeId < 0) {
                        return makeErrorMessage(ErrorCode::InvalidArgument,"Invalid link failure edge ID: " + token);
                    }

                    if (!seen.insert(edgeId).second) {
                        return makeErrorMessage(ErrorCode::InvalidArgument,"Duplicate link failure edge ID: " + token);
                    }

                    failedLinks.push_back(edgeId);
                } catch (const std::exception&) {
                    return makeErrorMessage(ErrorCode::InvalidArgument,"Invalid link failure edge ID: " + token);
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
                return makeErrorMessage(ErrorCode::InvalidArgument,"--visualization must not be empty.");
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


        // N-1 is a separate analysis operation. Reject irrelevant solve-only flags
        // rather than silently ignoring user input.
        if (isN1) {
            if (solvers.value().size() != 1) {
                return makeErrorMessage(ErrorCode::InvalidArgument,"analyze-n1 requires exactly one --solver.");
            }
            if (demands.size() != 1) {
                return makeErrorMessage(ErrorCode::InvalidArgument,"analyze-n1 requires exactly one --demand model.");
            }
            if (outputFilename.empty() || outputFormat != OutputFormat::JSON ||
                inferOutputFormat(outputFilename) != OutputFormat::JSON) {
                return makeErrorMessage(ErrorCode::InvalidArgument,"analyze-n1 requires --output <filename.json>.");
            }
            if (variables.count("load-routing-table") ||
                variables.count("save-routing-table") ||
                variables.count("link-failure") ||
                variables.count("visualization") ||
                failure_analysis || failure_recovery) {
                            return makeErrorMessage(
                                ErrorCode::InvalidArgument,
                                "analyze-n1 does not accept routing-table, "
                                "link-failure, visualization, "
                                "or legacy failure-analysis flags."
                            );
                }
        }

        if (isTraffic) {
            if (solvers.value().size() != 1 || demands.size() != 1) {
                return makeErrorMessage(
                    ErrorCode::InvalidArgument,
                    "analyze-traffic requires exactly one solver and demand model."
                );
            }

            if (routingTableInput.empty()) {
                return makeErrorMessage(
                    ErrorCode::InvalidArgument,
                    "analyze-traffic requires --load-routing-table."
                );
            }

            if (outputFilename.empty() ||
                inferOutputFormat(outputFilename) != OutputFormat::JSON ||
                outputFormat != OutputFormat::JSON) {
                return makeErrorMessage(
                    ErrorCode::InvalidArgument,
                    "analyze-traffic requires --output <filename.json>."
                );
                }

            if (!upgrades.empty() || !failedLinks.empty() ||
                failure_analysis || failure_recovery ||
                variables.count("save-routing-table") ||
                variables.count("visualization") ||
                variables.count("execution-mode")) {
                return makeErrorMessage(
                    ErrorCode::InvalidArgument,
                    "analyze-traffic does not accept topology, failure, "
                    "routing-export, visualization, or execution-mode options."
                );
                }
        }


        std::vector<CommodityScale> commodityScales;

        if (variables.count("commodity-scale")) {
            // The solve command also supports workspace demand modifications.
            if (!isTraffic && !isN1 && first != "solve") {
                return makeErrorMessage(
                    ErrorCode::InvalidArgument,
                    "--commodity-scale requires solve, analyze-traffic or analyze-n1"
                );
            }

            std::set<std::pair<int, int>> seen;

            for (const auto& token :
                 variables["commodity-scale"].as<std::vector<std::string>>()) {

                const auto first = token.find(':');
                const auto second = token.find(
                    ':',
                    first == std::string::npos ? first : first + 1
                );

                if (first == std::string::npos ||
                    second == std::string::npos ||
                    token.find(':', second + 1) != std::string::npos) {
                    return makeErrorMessage(
                        ErrorCode::InvalidArgument,
                        "Expected --commodity-scale source:target:factor: " + token
                    );
                }

                try {
                    const auto sourceText = token.substr(0, first);
                    const auto targetText =
                        token.substr(first + 1, second - first - 1);
                    const auto factorText = token.substr(second + 1);

                    std::size_t sourceUsed = 0;
                    std::size_t targetUsed = 0;
                    std::size_t factorUsed = 0;

                    const int source = std::stoi(sourceText, &sourceUsed);
                    const int target = std::stoi(targetText, &targetUsed);
                    const double factor = std::stod(factorText, &factorUsed);

                    if (sourceUsed != sourceText.size() ||
                        targetUsed != targetText.size() ||
                        factorUsed != factorText.size() ||
                        source < 0 ||
                        target < 0 ||
                        source == target ||
                        !std::isfinite(factor) ||
                        factor <= 0.0) {
                        return makeErrorMessage(
                            ErrorCode::InvalidArgument,
                            "Invalid commodity scale: " + token
                        );
                    }

                    if (!seen.insert({source, target}).second) {
                        return makeErrorMessage(
                            ErrorCode::InvalidArgument,
                            "Duplicate commodity override: " + token
                        );
                    }

                    commodityScales.push_back({
                        source, target, factor
                    });

                } catch (const std::exception&) {
                    return makeErrorMessage(
                        ErrorCode::InvalidArgument,
                        "Invalid commodity scale: " + token
                    );
                }
            }
        }


        /*
         * The CLI's job ends here.
         *
         * Everything below the CLI receives a normal
         * Config object and does not need to know about
         * argv / flags / subcommands.
         */
        Config config{ .solvers = std::move(solvers.value()),
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
            .visualization_output_directory =std::move(visualizationDirectory),
            .demand_scale = demandScale,
            .capacity_upgrades = std::move(upgrades),
            .failed_links = std::move(failedLinks)
        };
        return Command{
        .action = isN1 ? Action::AnalyzeN1
                   : (isTraffic ? Action::AnalyzeTraffic : Action::Run),
            .config = std::move(config),
            .execution_mode = executionMode,
            .commodity_scales = std::move(commodityScales)
        };
    } catch (const boost::program_options::error& error) {
        return makeErrorMessage(ErrorCode::InvalidArgument, error.what());
    } catch (const std::exception& error) {
        return fromStdException(error,ErrorCode::InvalidArgument);
    }

}



#endif //E_ROUTING_CLI_H


