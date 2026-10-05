//
// Created by Mert on 05.10.26.
//

#ifndef E_ROUTING_CLI_H
#define E_ROUTING_CLI_H


#include "parse_argument_io.h"
#include <boost/program_options.hpp>
#include "core/config.h"


enum class Action {
    Run,
    Help,
    ListSolvers
};

struct Command {
    Action action = Action::Run;
    Config config{};
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
        << "E-Routing\n"
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

    if (first != "solve") {
        return makeErrorMessage(ErrorCode::InvalidArgument,"Unknown command: " + first +". Use 'solve', '--help', ""or '--list-solvers'.");
    }


    try {
        const auto options = solveOptions();
        boost::program_options::variables_map variables;


        /*
         * Boost.Program_options treats argv[0] as the program name.
         *
         * We pass the "solve" subcommand as that program name,
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
            .visualization_output_directory =std::move(visualizationDirectory)
        };
        return Command{ .action = Action::Run,.config = std::move(config)};
    } catch (const boost::program_options::error& error) {
        return makeErrorMessage(ErrorCode::InvalidArgument, error.what());
    } catch (const std::exception& error) {
        return fromStdException(error,ErrorCode::InvalidArgument);
    }
}



#endif //E_ROUTING_CLI_H
