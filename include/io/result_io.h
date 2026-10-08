//
// Created by Mert Biyikli on 06.07.26.
//

#ifndef OBLIVIOUSROUTING_RESULT_IO_H
#define OBLIVIOUSROUTING_RESULT_IO_H

#include "../routing/routing_result.h"
#include "core/errors.h"
#include <cctype>
#include <filesystem>
#include <format>
#include <fstream>

class RoutingResultWriter {
  public:
    static Result<void> write(const RoutingExperimentResult& exp, const Config& cfg) {
        /*
         * stdout requires no file handling.
         */
        if (cfg.output_format == OutputFormat::COUT) {
            for (const auto& result : exp.solver_results) {
                writeCout(result);
                std::cout << "\n";
            }
            return {};
        }

        /*
         * Determine extension from the selected output format.
         */
        std::string extension;

        switch (cfg.output_format) {
        case OutputFormat::JSON:
            extension = ".json";
            break;

        case OutputFormat::TEXT:
            extension = ".txt";
            break;

        default:
            return makeErrorMessage(ErrorCode::FormatNotFound, "Unsupported output format.");
        }

        /*
         * One experiment now produces one output file,
         * regardless of the number of solvers.
         */
        std::filesystem::path output_path;

        if (cfg.output_filename.empty()) {
            output_path = std::filesystem::path("result") / ("run" + extension);
        } else {
            output_path = std::filesystem::path(cfg.output_filename);
            output_path.replace_extension(extension);
        }

        if (output_path.has_parent_path()) {
            std::error_code ec;

            std::filesystem::create_directories(output_path.parent_path(), ec);

            if (ec) {
                return makeErrorMessage(ErrorCode::FileNotFound, "Failed to create output directory: " + output_path.parent_path().string() + " (" + ec.message() + ")");
            }
        }

        std::ofstream file(output_path, std::ios::out | std::ios::trunc);

        if (!file.is_open()) {
            return makeErrorMessage(ErrorCode::FileNotFound, "Failed to open output file: " + output_path.string());
        }

        switch (cfg.output_format) {
        case OutputFormat::TEXT:
            for (const auto& result : exp.solver_results) {
                writeText(result, file);
                file << '\n';
            }
            break;

        case OutputFormat::JSON:
            writeJson(exp, cfg, file);
            break;

        default:
            return makeErrorMessage(ErrorCode::FormatNotFound, "Unsupported output format.");
        }

        file.flush();

        if (!file.good()) {
            return makeErrorMessage(ErrorCode::RuntimeError, "Failed while writing output file: " + output_path.string());
        }

        return {};
    }

  private:
    static void writeText(const IRoutingResult& result, std::ostream& out) {
        out << "Date: " << std::chrono::system_clock::now() << '\n';
        out << "Solver: " << result.solver_name << '\n';
        out << "Graph: " << result.graph_name << '\n';
        out << "Nodes: " << result.nodes << '\n';
        out << "Edges: " << result.edges << '\n';

        if (!result.routing_base.empty()) {
            out << "Routing base: " << result.routing_base << '\n';
        }

        out << "Status: " << static_cast<int>(result.status) << '\n';
        out << "Total runtime (microseconds): " << result.total_runtime_microseconds << '\n';

        if (result.preprocessing_runtime_microseconds >= 0.0) {
            out << "Preprocessing runtime (microseconds): " << result.preprocessing_runtime_microseconds << '\n';
        }

        if (result.solve_runtime_microseconds >= 0.0) {
            out << "Solve runtime (microseconds): " << result.solve_runtime_microseconds << '\n';
        }

        out << "Oblivious ratio: " << result.oblivious_ratio << '\n';

        if (result.candidate_paths > 0) {
            out << "Candidate paths: " << result.candidate_paths << '\n';
            out << "Average paths per pair: " << result.average_paths_per_pair << '\n';
        }

        for (const auto& eval : result.demand_evaluations) {
            out << "Demand [" << demandModelName(eval.demand_type) << "]\n";

            out << "  Congestion: " << eval.congestion << '\n';

            if (eval.runtime_microseconds >= 0.0) {
                out << "  Runtime (microseconds): " << eval.runtime_microseconds << '\n';
            }
            out << "  Link failure analysis:\n";
            out << "    Tested links: " << eval.failure_tested_links << '\n';
            out << "    Most critical edge ID: " << eval.failure_most_critical_edge_id << '\n';
            out << "    Most critical link: " << eval.failure_most_critical_source << " -> " << eval.failure_most_critical_target << '\n';
            out << "    Maximum lost traffic fraction: " << eval.failure_maximum_lost_traffic_fraction << '\n';
            out << "    Average lost traffic fraction: " << eval.failure_average_lost_traffic_fraction << '\n';
            out << "    Median lost traffic fraction: " << eval.failure_median_lost_traffic_fraction << '\n';
            out << "    Maximum affected demand fraction: " << eval.failure_maximum_affected_demand_fraction << '\n';
            out << "    Traffic-carrying links: " << eval.failure_traffic_carrying_links << '\n';
            out << "    Critical links >= 10%: " << eval.failure_critical_links_10_percent << '\n';
            out << "    Critical links >= 25%: " << eval.failure_critical_links_25_percent << '\n';
            out << "    Critical links >= 50%: " << eval.failure_critical_links_50_percent << '\n';
        }

        if (!result.mwu_metrics.empty()) {
            const auto& mwu = result.mwu_metrics;

            out << "MWU metrics:\n";
            out << "  Iterations: " << mwu.iteration_count << '\n';
            out << "  Solve time (microseconds): " << mwu.solve_time << '\n';
            out << "  Transformation time (microseconds): " << mwu.transformation_time << '\n';
            out << "  Load computation time (microseconds): " << mwu.load_computation_time << '\n';
            out << "  Weight update time (microseconds): " << mwu.mwu_weight_update_time << '\n';

            if (!mwu.oracle_running_times.empty()) {
                out << "  Average oracle time (microseconds): " << mwu.averageOracleTime() << '\n';
            }
        }
        if (!result.expander_metrics.empty()) {
            const auto& metrics = result.expander_metrics;

            out << "Expander hierarchy metrics:\n";

            out << "  Hierarchy runtime (microseconds): " << metrics.hierarchy_runtime_microseconds << '\n';

            out << "  Tree construction runtime (microseconds): " << metrics.tree_runtime_microseconds << '\n';

            out << "  Basis-flow runtime (microseconds): " << metrics.basis_flow_runtime_microseconds << '\n';

            out << "  Hierarchy levels: " << metrics.hierarchy_levels << '\n';

            out << "  Hierarchy clusters: " << metrics.hierarchy_clusters << '\n';

            out << "  Clusters per level: ";

            for (std::size_t index = 0; index < metrics.clusters_per_level.size(); ++index) {
                if (index > 0) {
                    out << ", ";
                }

                out << metrics.clusters_per_level[index];
            }

            out << '\n';

            out << "  Maximum cluster vertices: " << metrics.max_cluster_vertices << '\n';

            out << "  Average cluster vertices: " << metrics.average_cluster_vertices << '\n';

            out << "  Tree nodes: " << metrics.tree_nodes << '\n';

            out << "  Tree edges: " << metrics.tree_edges << '\n';

            out << "  Tree depth: " << metrics.tree_depth << '\n';

            out << "  Basis flows: " << metrics.basis_flows << '\n';

            out << "  Total electrical solves: " << metrics.total_electrical_solves << '\n';

            out << "  Average electrical solves per basis flow: " << metrics.averageElectricalSolves() << '\n';

            out << "  Maximum basis embedding congestion: " << metrics.max_basis_embedding_congestion << '\n';

            out << "  Maximum conservation error: " << metrics.max_conservation_error << '\n';
        }
    }

    static void writeCout(const IRoutingResult& result) {
        std::cout << "Date: " << std::chrono::system_clock::now() << '\n';
        std::cout << "Solver: " << result.solver_name << '\n';
        std::cout << "Graph: " << result.graph_name << '\n';
        std::cout << "Nodes: " << result.nodes << '\n';
        std::cout << "Edges: " << result.edges << '\n';

        if (!result.routing_base.empty()) {
            std::cout << "Routing base: " << result.routing_base << '\n';
        }

        std::cout << "Status: " << static_cast<int>(result.status) << '\n';
        std::cout << "Total runtime (microseconds): " << result.total_runtime_microseconds << '\n';

        if (result.preprocessing_runtime_microseconds >= 0.0) {
            std::cout << "Preprocessing runtime (microseconds): " << result.preprocessing_runtime_microseconds << '\n';
        }

        if (result.solve_runtime_microseconds >= 0.0) {
            std::cout << "Solve runtime (microseconds): " << result.solve_runtime_microseconds << '\n';
        }

        std::cout << "Oblivious ratio: " << result.oblivious_ratio << '\n';

        if (result.candidate_paths > 0) {
            std::cout << "Candidate paths: " << result.candidate_paths << '\n';
            std::cout << "Average paths per pair: " << result.average_paths_per_pair << '\n';
        }
        for (const auto& eval : result.demand_evaluations) {
            std::cout << "Demand [" << demandModelName(eval.demand_type) << "]\n";
            std::cout << "  Congestion: " << eval.congestion << '\n';

            if (eval.runtime_microseconds >= 0.0) {
                std::cout << "  Runtime (microseconds): " << eval.runtime_microseconds << '\n';
            }
            std::cout << "  Link failure analysis:\n";
            std::cout << "    Tested links: " << eval.failure_tested_links << '\n';
            std::cout << "    Most critical edge ID: " << eval.failure_most_critical_edge_id << '\n';
            std::cout << "    Most critical link: " << eval.failure_most_critical_source << " -> " << eval.failure_most_critical_target << '\n';
            std::cout << "    Maximum lost traffic fraction: " << eval.failure_maximum_lost_traffic_fraction << '\n';
            std::cout << "    Average lost traffic fraction: " << eval.failure_average_lost_traffic_fraction << '\n';
            std::cout << "    Median lost traffic fraction: " << eval.failure_median_lost_traffic_fraction << '\n';
            std::cout << "    Maximum affected demand fraction: " << eval.failure_maximum_affected_demand_fraction << '\n';
            std::cout << "    Traffic-carrying links: " << eval.failure_traffic_carrying_links << '\n';
            std::cout << "    Critical links >= 10%: " << eval.failure_critical_links_10_percent << '\n';
            std::cout << "    Critical links >= 25%: " << eval.failure_critical_links_25_percent << '\n';
            std::cout << "    Critical links >= 50%: " << eval.failure_critical_links_50_percent << '\n';
        }

        if (!result.mwu_metrics.empty()) {
            const auto& mwu = result.mwu_metrics;

            std::cout << "MWU metrics:\n";
            std::cout << "  Iterations: " << mwu.iteration_count << '\n';
            std::cout << "  Solve time (microseconds): " << mwu.solve_time << '\n';
            std::cout << "  Transformation time (microseconds): " << mwu.transformation_time << '\n';
            std::cout << "  Load computation time (microseconds): " << mwu.load_computation_time << '\n';
            std::cout << "  Weight update time (microseconds): " << mwu.mwu_weight_update_time << '\n';

            if (!mwu.oracle_running_times.empty()) {
                std::cout << "  Average oracle time (microseconds): " << mwu.averageOracleTime() << '\n';

                std::cout << "  Oracle calls: " << mwu.oracle_running_times.size() << '\n';
            }
        }
        if (!result.expander_metrics.empty()) {
            const auto& metrics = result.expander_metrics;

            std::cout << "Expander hierarchy metrics:\n";

            std::cout << "  Hierarchy runtime (microseconds): " << metrics.hierarchy_runtime_microseconds << '\n';

            std::cout << "  Tree construction runtime (microseconds): " << metrics.tree_runtime_microseconds << '\n';

            std::cout << "  Basis-flow runtime (microseconds): " << metrics.basis_flow_runtime_microseconds << '\n';

            std::cout << "  Hierarchy levels: " << metrics.hierarchy_levels << '\n';

            std::cout << "  Hierarchy clusters: " << metrics.hierarchy_clusters << '\n';

            std::cout << "  Clusters per level: ";

            for (std::size_t index = 0; index < metrics.clusters_per_level.size(); ++index) {
                if (index > 0) {
                    std::cout << ", ";
                }

                std::cout << metrics.clusters_per_level[index];
            }

            std::cout << '\n';

            std::cout << "  Maximum cluster vertices: " << metrics.max_cluster_vertices << '\n';

            std::cout << "  Average cluster vertices: " << metrics.average_cluster_vertices << '\n';

            std::cout << "  Tree nodes: " << metrics.tree_nodes << '\n';

            std::cout << "  Tree edges: " << metrics.tree_edges << '\n';

            std::cout << "  Tree depth: " << metrics.tree_depth << '\n';

            std::cout << "  Basis flows: " << metrics.basis_flows << '\n';

            std::cout << "  Total electrical solves: " << metrics.total_electrical_solves << '\n';

            std::cout << "  Average electrical solves per basis flow: " << metrics.averageElectricalSolves() << '\n';

            std::cout << "  Maximum basis embedding congestion: " << metrics.max_basis_embedding_congestion << '\n';

            std::cout << "  Maximum conservation error: " << metrics.max_conservation_error << '\n';
        }
    }

    static std::string jsonEscape(const std::string& value) {
        std::string escaped;

        escaped.reserve(value.size());

        for (const char c : value) {
            switch (c) {
            case '"':
                escaped += "\\\"";
                break;

            case '\\':
                escaped += "\\\\";
                break;

            case '\b':
                escaped += "\\b";
                break;

            case '\f':
                escaped += "\\f";
                break;

            case '\n':
                escaped += "\\n";
                break;

            case '\r':
                escaped += "\\r";
                break;

            case '\t':
                escaped += "\\t";
                break;

            default:
                escaped += c;
                break;
            }
        }

        return escaped;
    }

    static void writeJson(const RoutingExperimentResult& experiment, const Config& cfg, std::ostream& out) {
        out << "{\n";

        /*
         * Schema metadata
         */
        out << "  \"schema_version\": \"1.0\",\n";

        out << "  \"timestamp\": \"" << jsonEscape(std::format("{}", std::chrono::system_clock::now())) << "\",\n";

        /*
         * Graph metadata
         */
        out << "  \"graph\": {\n";

        out << "    \"name\": \"" << jsonEscape(experiment.graph_name) << "\",\n";

        out << "    \"path\": \"" << jsonEscape(experiment.graph_path) << "\",\n";

        out << "    \"nodes\": " << experiment.nodes << ",\n";

        out << "    \"edges\": " << experiment.edges << "\n";

        out << "  },\n";

        /*
         * Run configuration
         */
        out << "  \"configuration\": {\n";

        out << "    \"seed\": " << cfg.seed << ",\n";

        out << "    \"threads\": " << cfg.num_threads << ",\n";

        out << "    \"demand_models\": [";

        for (std::size_t i = 0; i < cfg.demand_models.size(); ++i) {
            if (i > 0) {
                out << ", ";
            }

            out << "\"" << jsonEscape(demandModelName(cfg.demand_models[i])) << "\"";
        }

        out << "],\n";

        out << "    \"failure_recovery\": " << (cfg.failure_recovery ? "true" : "false") << "\n";

        out << "  },\n";

        /*
         * Solver results
         */
        out << "  \"solver_results\": [\n";

        for (std::size_t solver_index = 0; solver_index < experiment.solver_results.size(); ++solver_index) {
            const auto& result = experiment.solver_results[solver_index];

            out << "    {\n";

            /*
             * Solver identity
             */
            out << "      \"solver\": \"" << jsonEscape(result.solver_name) << "\",\n";

            out << "      \"solver_type\": \"" << jsonEscape(getSolverTypeName(result.type)) << "\",\n";

            out << "      \"status\": \"" << resultStatusName(result.status) << "\"";

            if (!result.routing_base.empty()) {
                out << ",\n";

                out << "      \"routing_base\": \"" << jsonEscape(result.routing_base) << "\"";
            }

            out << ",\n";

            /*
             * Runtime
             */
            out << "      \"runtime\": {\n";

            out << "        \"total_microseconds\": " << result.total_runtime_microseconds;

            if (result.preprocessing_runtime_microseconds >= 0.0) {
                out << ",\n";

                out << "        \"preprocessing_microseconds\": " << result.preprocessing_runtime_microseconds;
            }

            if (result.solve_runtime_microseconds >= 0.0) {
                out << ",\n";

                out << "        \"solve_microseconds\": " << result.solve_runtime_microseconds;
            }

            out << "\n";
            out << "      },\n";

            /*
             * Quality metrics
             */
            out << "      \"quality\": {\n";

            out << "        \"oblivious_ratio\": " << result.oblivious_ratio << "\n";

            out << "      },\n";

            /*
             * Routing scheme information
             */
            out << "      \"routing_scheme\": {\n";

            out << "        \"candidate_paths\": " << result.candidate_paths << ",\n";

            out << "        \"average_paths_per_pair\": " << result.average_paths_per_pair << "\n";

            out << "      },\n";

            /*
             * Demand evaluations
             */
            out << "      \"demand_evaluations\": [\n";

            for (std::size_t i = 0; i < result.demand_evaluations.size(); ++i) {
                const auto& eval = result.demand_evaluations[i];

                out << "        {\n";

                out << "          \"demand_model\": \"" << jsonEscape(demandModelName(eval.demand_type)) << "\",\n";

                out << "          \"congestion\": " << eval.congestion << ",\n";

                out << "          \"runtime_microseconds\": " << eval.runtime_microseconds << ",\n";

                /*
                 * ------------------------------------------------------------
                 * Layer 1: static link-failure exposure
                 * ------------------------------------------------------------
                 */
                out << "          \"failure_analysis\": {\n";

                out << "            \"runtime_microseconds\": " << eval.failure_analysis_runtime_microseconds << ",\n";

                out << "            \"tested_links\": " << eval.failure_tested_links << ",\n";

                out << "            \"most_critical_edge_id\": " << eval.failure_most_critical_edge_id << ",\n";

                out << "            \"most_critical_source\": " << eval.failure_most_critical_source << ",\n";

                out << "            \"most_critical_target\": " << eval.failure_most_critical_target << ",\n";

                out << "            \"maximum_lost_traffic_fraction\": " << eval.failure_maximum_lost_traffic_fraction << ",\n";

                out << "            \"average_lost_traffic_fraction\": " << eval.failure_average_lost_traffic_fraction << ",\n";

                out << "            \"median_lost_traffic_fraction\": " << eval.failure_median_lost_traffic_fraction << ",\n";

                out << "            \"maximum_affected_demand_fraction\": " << eval.failure_maximum_affected_demand_fraction << ",\n";

                out << "            \"traffic_carrying_links\": " << eval.failure_traffic_carrying_links << ",\n";

                out << "            \"critical_links_10_percent\": " << eval.failure_critical_links_10_percent << ",\n";

                out << "            \"critical_links_25_percent\": " << eval.failure_critical_links_25_percent << ",\n";

                out << "            \"critical_links_50_percent\": " << eval.failure_critical_links_50_percent << "\n";

                out << "          }";

                /*
                 * ------------------------------------------------------------
                 * Layer 2: actual N-1 failure + solver recomputation
                 * ------------------------------------------------------------
                 */
                if (eval.failure_recovery_available) {
                    out << ",\n";

                    out << "          \"failure_recovery\": {\n";

                    out << "            \"tested_links\": " << eval.recovery_tested_links << ",\n";

                    out << "            \"disconnected_failures\": " << eval.recovery_disconnected_failures << ",\n";

                    out << "            \"successful_recomputations\": " << eval.recovery_successful_recomputations << ",\n";

                    out << "            \"failed_recomputations\": " << eval.recovery_failed_recomputations << ",\n";

                    out << "            \"maximum_unroutable_demand_fraction\": " << eval.recovery_maximum_unroutable_demand_fraction << ",\n";

                    out << "            \"average_unroutable_demand_fraction\": " << eval.recovery_average_unroutable_demand_fraction << ",\n";

                    out << "            \"maximum_post_failure_congestion\": " << eval.recovery_maximum_post_failure_congestion << ",\n";

                    out << "            \"average_post_failure_congestion\": " << eval.recovery_average_post_failure_congestion << ",\n";

                    out << "            \"maximum_congestion_increase_factor\": " << eval.recovery_maximum_congestion_increase_factor << ",\n";

                    out << "            \"average_congestion_increase_factor\": " << eval.recovery_average_congestion_increase_factor << ",\n";

                    out << "            \"average_recomputation_runtime_microseconds\": " << eval.recovery_average_recomputation_runtime_microseconds << ",\n";

                    out << "            \"maximum_recomputation_runtime_microseconds\": " << eval.recovery_maximum_recomputation_runtime_microseconds << ",\n";

                    out << "            \"worst_failed_edge_id\": " << eval.recovery_worst_failed_edge_id << ",\n";

                    out << "            \"worst_failed_source\": " << eval.recovery_worst_failed_source << ",\n";

                    out << "            \"worst_failed_target\": " << eval.recovery_worst_failed_target << ",\n";

                    /*
                     * Worst disconnecting failure.
                     */
                    out << "            \"worst_disconnect_edge_id\": " << eval.recovery_worst_disconnect_edge_id << ",\n";

                    out << "            \"worst_disconnect_source\": " << eval.recovery_worst_disconnect_source << ",\n";

                    out << "            \"worst_disconnect_target\": " << eval.recovery_worst_disconnect_target << ",\n";

                    out << "            \"worst_disconnect_unroutable_demand_fraction\": " << eval.recovery_worst_disconnect_unroutable_demand_fraction << ",\n";

                    /*
                     * Worst survivable congestion failure.
                     */
                    out << "            \"worst_congestion_edge_id\": " << eval.recovery_worst_congestion_edge_id << ",\n";

                    out << "            \"worst_congestion_source\": " << eval.recovery_worst_congestion_source << ",\n";

                    out << "            \"worst_congestion_target\": " << eval.recovery_worst_congestion_target << ",\n";

                    out << "            \"worst_congestion_baseline\": " << eval.recovery_worst_congestion_baseline << ",\n";

                    out << "            \"worst_congestion_post_failure\": " << eval.recovery_worst_congestion_post_failure << ",\n";

                    out << "            \"worst_congestion_increase_factor\": " << eval.recovery_worst_congestion_increase_factor << ",\n";

                    /*
                     * Slowest successful recovery.
                     */
                    out << "            \"slowest_recovery_edge_id\": " << eval.recovery_slowest_recovery_edge_id << ",\n";

                    out << "            \"slowest_recovery_source\": " << eval.recovery_slowest_recovery_source << ",\n";

                    out << "            \"slowest_recovery_target\": " << eval.recovery_slowest_recovery_target << ",\n";

                    out << "            \"slowest_recovery_runtime_microseconds\": " << eval.recovery_slowest_recovery_runtime_microseconds << "\n";

                    out << "          }";
                }

                /*
                 * Close the complete demand-evaluation object only AFTER
                 * both Layer 1 and optional Layer 2 have been written.
                 */
                out << "\n";
                out << "        }";

                if (i + 1 < result.demand_evaluations.size()) {
                    out << ",";
                }

                out << "\n";
            }

            out << "      ],\n";

            /*
             * Solver-specific metrics.
             */
            out << "      \"algorithm_metrics\": {";

            bool wrote_algorithm_metrics = false;

            if (!result.mwu_metrics.empty()) {
                const auto& mwu = result.mwu_metrics;

                out << "\n";
                out << "        \"mwu\": {\n";

                out << "          \"iteration_count\": " << mwu.iteration_count << ",\n";

                out << "          \"solve_time_microseconds\": " << mwu.solve_time << ",\n";

                out << "          \"transformation_time_microseconds\": " << mwu.transformation_time << ",\n";

                out << "          \"load_computation_time_microseconds\": " << mwu.load_computation_time << ",\n";

                out << "          \"weight_update_time_microseconds\": " << mwu.mwu_weight_update_time;

                if (!mwu.oracle_running_times.empty()) {
                    out << ",\n";

                    out << "          \"average_oracle_time_microseconds\": " << mwu.averageOracleTime() << ",\n";

                    out << "          \"oracle_calls\": " << mwu.oracle_running_times.size() << "\n";
                } else {
                    out << "\n";
                }

                out << "        }";

                wrote_algorithm_metrics = true;
            }

            if (!result.expander_metrics.empty()) {
                const auto& metrics = result.expander_metrics;

                if (wrote_algorithm_metrics) {
                    out << ",";
                }

                out << "\n";
                out << "        \"expander\": {\n";

                out << "          \"hierarchy_runtime_microseconds\": " << metrics.hierarchy_runtime_microseconds << ",\n";

                out << "          \"tree_runtime_microseconds\": " << metrics.tree_runtime_microseconds << ",\n";

                out << "          \"basis_flow_runtime_microseconds\": " << metrics.basis_flow_runtime_microseconds << ",\n";

                out << "          \"hierarchy_levels\": " << metrics.hierarchy_levels << ",\n";

                out << "          \"hierarchy_clusters\": " << metrics.hierarchy_clusters << ",\n";

                out << "          \"clusters_per_level\": [";

                for (std::size_t i = 0; i < metrics.clusters_per_level.size(); ++i) {
                    if (i > 0) {
                        out << ", ";
                    }

                    out << metrics.clusters_per_level[i];
                }

                out << "],\n";

                out << "          \"max_cluster_vertices\": " << metrics.max_cluster_vertices << ",\n";

                out << "          \"average_cluster_vertices\": " << metrics.average_cluster_vertices << ",\n";

                out << "          \"tree_nodes\": " << metrics.tree_nodes << ",\n";

                out << "          \"tree_edges\": " << metrics.tree_edges << ",\n";

                out << "          \"tree_depth\": " << metrics.tree_depth << ",\n";

                out << "          \"basis_flows\": " << metrics.basis_flows << ",\n";

                out << "          \"total_electrical_solves\": " << metrics.total_electrical_solves << ",\n";

                out << "          \"average_electrical_solves\": " << metrics.averageElectricalSolves() << ",\n";

                out << "          \"max_basis_embedding_congestion\": " << metrics.max_basis_embedding_congestion << ",\n";

                out << "          \"max_conservation_error\": " << metrics.max_conservation_error << "\n";

                out << "        }";

                wrote_algorithm_metrics = true;
            }

            if (wrote_algorithm_metrics) {
                out << "\n";
                out << "      }\n";
            } else {
                out << "}\n";
            }

            out << "    }";

            if (solver_index + 1 < experiment.solver_results.size()) {
                out << ",";
            }

            out << "\n";
        }

        out << "  ]\n";
        out << "}\n";
    }

    static void writeExpanderMetricsJson(const ExpanderMetrics& metrics, std::ostream& out) {
        out << "  \"expander_metrics\": {\n";

        out << "    \"hierarchy_runtime_microseconds\": " << metrics.hierarchy_runtime_microseconds << ",\n";

        out << "    \"tree_runtime_microseconds\": " << metrics.tree_runtime_microseconds << ",\n";

        out << "    \"basis_flow_runtime_microseconds\": " << metrics.basis_flow_runtime_microseconds << ",\n";

        out << "    \"hierarchy_levels\": " << metrics.hierarchy_levels << ",\n";

        out << "    \"hierarchy_clusters\": " << metrics.hierarchy_clusters << ",\n";

        out << "    \"clusters_per_level\": [";

        for (std::size_t index = 0; index < metrics.clusters_per_level.size(); ++index) {
            if (index > 0) {
                out << ", ";
            }

            out << metrics.clusters_per_level[index];
        }

        out << "],\n";

        out << "    \"max_cluster_vertices\": " << metrics.max_cluster_vertices << ",\n";

        out << "    \"average_cluster_vertices\": " << metrics.average_cluster_vertices << ",\n";

        out << "    \"tree_nodes\": " << metrics.tree_nodes << ",\n";

        out << "    \"tree_edges\": " << metrics.tree_edges << ",\n";

        out << "    \"tree_depth\": " << metrics.tree_depth << ",\n";

        out << "    \"basis_flows\": " << metrics.basis_flows << ",\n";

        out << "    \"total_electrical_solves\": " << metrics.total_electrical_solves << ",\n";

        out << "    \"average_electrical_solves\": " << metrics.averageElectricalSolves() << ",\n";

        out << "    \"max_basis_embedding_congestion\": " << metrics.max_basis_embedding_congestion << ",\n";

        out << "    \"max_conservation_error\": " << metrics.max_conservation_error << "\n";

        out << "  }\n";
    }

    // Helpers
    static std::string safeFileName(std::string name) {
        for (char& c : name) {
            const auto value = static_cast<unsigned char>(c);

            if (!std::isalnum(value) && c != '-' && c != '_') {
                c = '_';
            }
        }

        return name;
    }

    static std::string resultStatusName(ResultStatus status) {
        switch (status) {
        case ResultStatus::OK:
            return "ok";
        case ResultStatus::ERROR_INVALID_SOLVER:
            return "error: invalid solver";
        case ResultStatus::ERROR_INVALID_ROUTING_SCHEME:
            return "error: invalid routing scheme";
        case ResultStatus::ERROR_MISSING_DEMAND_MODELS:
            return "error: missing demand models";
        default:
            return "error";
        }
    }
};

#endif // OBLIVIOUSROUTING_RESULT_IO_H
