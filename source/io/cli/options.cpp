//
// Created by Mert on 30.10.26.
//

#include "io/cli/options.h"

boost::program_options::options_description solveOptions() {
    boost::program_options::options_description options("Solve options");

    options.add_options()

        ("help,h", "Show solve help")

            ("solver,s", boost::program_options::value<std::string>(), "Solver name or comma-separated solver list")

                ("graph,g", boost::program_options::value<std::string>(), "Input graph file")

                    ("load-routing-table", boost::program_options::value<std::string>(), "Load a serialized routing scheme instead of executing the solver")

                        ("demand,d", boost::program_options::value<std::string>(), "Demand models: uniform,gravity,bimodal,gaussian")

                            ("graph-format", boost::program_options::value<std::string>()->default_value("csr"), "Graph representation: csr or adjlist")

                                ("threads,j", boost::program_options::value<int>()->default_value(1), "Positive worker-thread count")

                                    ("output,o", boost::program_options::value<std::string>(), "Output file")

                                        ("output-format", boost::program_options::value<std::string>(), "Output format: json, txt, or cout")

                                            ("seed", boost::program_options::value<int>()->default_value(42), "Non-negative random seed")

                                                ("execution-mode", boost::program_options::value<std::string>(), "N-1 execution mode: fast, standard, research")

                                                    ("failure-analysis", boost::program_options::bool_switch()->default_value(false), "Run Layer-1 static N-1 link-failure exposure analysis")

                                                        ("failure-recovery", boost::program_options::bool_switch()->default_value(false), "Run Layer-2 N-1 routing recomputation and recovery analysis")

                                                            ("save-routing-table", boost::program_options::value<std::string>(),
                                                             "Save computed routing tables as binary files in the specified directory")

                                                                ("demand-scale", boost::program_options::value<double>()->default_value(1.0),
                                                                 "Demand multiplier for persisted routing analysis")("commodity-scale",
                                                                                                                     boost::program_options::value<std::vector<std::string>>()
                                                                                                                         ->multitoken()
                                                                                                                         ->composing(),
                                                                                                                     "Directed commodity multiplier: source:target:factor")(
                                                                    "capacity-upgrade", boost::program_options::value<std::vector<std::string>>()->multitoken()->composing(),
                                                                    "Physical link capacity multiplier as edge_id:factor")("link-failure",
                                                                                                                           boost::program_options::value<std::vector<std::string>>()
                                                                                                                               ->multitoken()
                                                                                                                               ->composing(),
                                                                                                                           "Simulate physical link failures using edge IDs, "
                                                                                                                           "e.g. --link-failure 14 20 24")

                                                                    ("visualization", boost::program_options::value<std::string>(),
                                                                     "Visualization JSON output directory")("analysis", boost::program_options::value<std::string>(),
                                                                                                            "Scenario analysis: link-failure (capacity comparison planned)");

    return options;
}
