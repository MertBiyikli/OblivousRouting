//
// Created by Mert Biyikli on 06.07.26.
//

#ifndef OBLIVIOUSROUTING_RESULT_IO_H
#define OBLIVIOUSROUTING_RESULT_IO_H

#include "core/errors.h"
#include <fstream>
#include "../routing/routing_result.h"




class RoutingResultWriter {
public:
    static Result<void> write(
        const IRoutingResult& result,
        const std::string& path,
        OutputFormat format
    ) {

        if (format == OutputFormat::COUT) {
            writeCout(result);
            return {};
        }

        // else we write into a file
        std::ofstream file(path);

        if (!file.is_open()) {
            return makeErrorMessage(ErrorCode::FileNotFound, "Failed to open output file: "+ path);
        }

        switch (format) {
            case OutputFormat::TEXT:
                writeText(result, file);
                return {};
            case OutputFormat::JSON:
                writeJson(result, file);
                return {};
            default:
                return makeErrorMessage(ErrorCode::FormatNotFound, "Unsupported output format.");
        }

    }

private:
    static void writeText(
        const IRoutingResult& result,
        std::ostream& out
    ) {
        out << "Date: " << std::chrono::system_clock::now() << '\n';
        out << "Solver: " << result.solver_name << '\n';
        out << "Graph: " << result.graph_name << '\n';
        out << "Nodes: " << result.nodes << '\n';
        out << "Edges: " << result.edges << '\n';

        if (!result.routing_base.empty()) {
            out << "Routing base: " << result.routing_base << '\n';
        }

        out << "Status: " << static_cast<int>(result.status) << '\n';
        out << "Total runtime (microseconds): "
            << result.total_runtime_microseconds << '\n';

        if (result.preprocessing_runtime_microseconds >= 0.0) {
            out << "Preprocessing runtime (microseconds): "
                << result.preprocessing_runtime_microseconds << '\n';
        }

        if (result.solve_runtime_microseconds >= 0.0) {
            out << "Solve runtime (microseconds): "
                << result.solve_runtime_microseconds << '\n';
        }

        out << "Oblivious ratio: " << result.oblivious_ratio << '\n';

        if (result.candidate_paths > 0) {
            out << "Candidate paths: " << result.candidate_paths << '\n';
            out << "Average paths per pair: "
                << result.average_paths_per_pair << '\n';
        }

        for (const auto& eval : result.demand_evaluations) {
            out << "Demand [" << demandModelName(eval.demand_type) << "]\n";
            out << "  Congestion: " << eval.congestion << '\n';


            if (eval.runtime_microseconds >= 0.0) {
                out << "  Runtime (microseconds): "
                    << eval.runtime_microseconds << '\n';
            }
        }

        if (!result.mwu_metrics.empty()) {
            const auto& mwu = result.mwu_metrics;

            out << "MWU metrics:\n";
            out << "  Iterations: " << mwu.iteration_count << '\n';
            out << "  Solve time (microseconds): " << mwu.solve_time << '\n';
            out << "  Transformation time (microseconds): "
                << mwu.transformation_time << '\n';
            out << "  Load computation time (microseconds): "
                << mwu.load_computation_time << '\n';
            out << "  Weight update time (microseconds): "
                << mwu.mwu_weight_update_time << '\n';

            if (!mwu.oracle_running_times.empty()) {
                out << "  Average oracle time (microseconds): "
                    << mwu.averageOracleTime() << '\n';
            }
        }
        if (!result.expander_metrics.empty()) {
            const auto& metrics = result.expander_metrics;

            out << "Expander hierarchy metrics:\n";

            out << "  Hierarchy runtime (microseconds): "
                << metrics.hierarchy_runtime_microseconds
                << '\n';

            out << "  Tree construction runtime (microseconds): "
                << metrics.tree_runtime_microseconds
                << '\n';

            out << "  Basis-flow runtime (microseconds): "
                << metrics.basis_flow_runtime_microseconds
                << '\n';

            out << "  Hierarchy levels: "
                << metrics.hierarchy_levels
                << '\n';

            out << "  Hierarchy clusters: "
                << metrics.hierarchy_clusters
                << '\n';

            out << "  Clusters per level: ";

            for (std::size_t index = 0; index < metrics.clusters_per_level.size(); ++index) {
                if (index > 0) {
                    out << ", ";
                }

                out << metrics.clusters_per_level[index];
                 }

            out << '\n';

            out << "  Maximum cluster vertices: "
                << metrics.max_cluster_vertices
                << '\n';

            out << "  Average cluster vertices: "
                << metrics.average_cluster_vertices
                << '\n';

            out << "  Tree nodes: "
                << metrics.tree_nodes
                << '\n';

            out << "  Tree edges: "
                << metrics.tree_edges
                << '\n';

            out << "  Tree depth: "
                << metrics.tree_depth
                << '\n';

            out << "  Basis flows: "
                << metrics.basis_flows
                << '\n';

            out << "  Total electrical solves: "
                << metrics.total_electrical_solves
                << '\n';

            out << "  Average electrical solves per basis flow: "
                << metrics.averageElectricalSolves()
                << '\n';

            out << "  Maximum basis embedding congestion: "
                << metrics.max_basis_embedding_congestion
                << '\n';

            out << "  Maximum conservation error: "
                << metrics.max_conservation_error
                << '\n';
        }
    }

    static void writeCout(
        const IRoutingResult& result) {
        std::cout << "Date: " << std::chrono::system_clock::now() << '\n';
        std::cout << "Solver: " << result.solver_name << '\n';
        std::cout << "Graph: " << result.graph_name << '\n';
        std::cout << "Nodes: " << result.nodes << '\n';
        std::cout << "Edges: " << result.edges << '\n';

        if (!result.routing_base.empty()) {
            std::cout << "Routing base: " << result.routing_base << '\n';
        }

        std::cout << "Status: " << static_cast<int>(result.status) << '\n';
        std::cout << "Total runtime (microseconds): "
            << result.total_runtime_microseconds << '\n';

        if (result.preprocessing_runtime_microseconds >= 0.0) {
            std::cout << "Preprocessing runtime (microseconds): "
                << result.preprocessing_runtime_microseconds << '\n';
        }

        if (result.solve_runtime_microseconds >= 0.0) {
            std::cout << "Solve runtime (microseconds): "
                << result.solve_runtime_microseconds << '\n';
        }

        std::cout << "Oblivious ratio: " << result.oblivious_ratio << '\n';

        if (result.candidate_paths > 0) {
            std::cout << "Candidate paths: " << result.candidate_paths << '\n';
            std::cout << "Average paths per pair: "
                << result.average_paths_per_pair << '\n';
        }
        for (const auto& eval : result.demand_evaluations) {
            std::cout << "Demand [" << demandModelName(eval.demand_type) << "]\n";
            std::cout << "  Congestion: " << eval.congestion << '\n';


            if (eval.runtime_microseconds >= 0.0) {
                std::cout << "  Runtime (microseconds): "
                    << eval.runtime_microseconds << '\n';
            }
        }

        if (!result.mwu_metrics.empty()) {
            const auto& mwu = result.mwu_metrics;

            std::cout << "MWU metrics:\n";
            std::cout << "  Iterations: " << mwu.iteration_count << '\n';
            std::cout << "  Solve time (microseconds): " << mwu.solve_time << '\n';
            std::cout << "  Transformation time (microseconds): "
                << mwu.transformation_time << '\n';
            std::cout << "  Load computation time (microseconds): "
                << mwu.load_computation_time << '\n';
            std::cout << "  Weight update time (microseconds): "
                << mwu.mwu_weight_update_time << '\n';

            if (!mwu.oracle_running_times.empty()) {
                std::cout << "  Average oracle time (microseconds): "
                    << mwu.averageOracleTime() << '\n';

                std::cout << "  Oracle calls: "
                    << mwu.oracle_running_times.size() << '\n';
            }
        }
        if (!result.expander_metrics.empty()) {
            const auto& metrics = result.expander_metrics;

            std::cout << "Expander hierarchy metrics:\n";

            std::cout << "  Hierarchy runtime (microseconds): "
                << metrics.hierarchy_runtime_microseconds
                << '\n';

            std::cout << "  Tree construction runtime (microseconds): "
                << metrics.tree_runtime_microseconds
                << '\n';

            std::cout << "  Basis-flow runtime (microseconds): "
                << metrics.basis_flow_runtime_microseconds
                << '\n';

            std::cout << "  Hierarchy levels: "
                << metrics.hierarchy_levels
                << '\n';

            std::cout << "  Hierarchy clusters: "
                << metrics.hierarchy_clusters
                << '\n';

            std::cout << "  Clusters per level: ";

            for (std::size_t index = 0; index < metrics.clusters_per_level.size(); ++index) {
                if (index > 0) {
                    std::cout << ", ";
                }

                std::cout << metrics.clusters_per_level[index];
                 }

            std::cout << '\n';

            std::cout << "  Maximum cluster vertices: "
                << metrics.max_cluster_vertices
                << '\n';

            std::cout << "  Average cluster vertices: "
                << metrics.average_cluster_vertices
                << '\n';

            std::cout << "  Tree nodes: "
                << metrics.tree_nodes
                << '\n';

            std::cout << "  Tree edges: "
                << metrics.tree_edges
                << '\n';

            std::cout << "  Tree depth: "
                << metrics.tree_depth
                << '\n';

            std::cout << "  Basis flows: "
                << metrics.basis_flows
                << '\n';

            std::cout << "  Total electrical solves: "
                << metrics.total_electrical_solves
                << '\n';

            std::cout << "  Average electrical solves per basis flow: "
                << metrics.averageElectricalSolves()
                << '\n';

            std::cout << "  Maximum basis embedding congestion: "
                << metrics.max_basis_embedding_congestion
                << '\n';

            std::cout << "  Maximum conservation error: "
                << metrics.max_conservation_error
                << '\n';
        }
    }

    static std::string jsonEscape(const std::string& s) {
        std::string escaped;
        escaped.reserve(s.size());

        for (char c : s) {
            switch (c) {
                case '"': escaped += "\\\""; break;
                case '\\': escaped += "\\\\"; break;
                case '\n': escaped += "\\n"; break;
                case '\t': escaped += "\\t"; break;
                default: escaped += c; break;
            }
        }

        return escaped;
    }

    static void writeJson(const IRoutingResult& result,std::ostream& out) {
        out << "{\n";
        out << R"(  "date": ")" << std::chrono::system_clock::now() << "\",\n";
        out << R"(  "solver": ")" << (result.solver_name) << "\",\n";
        out << "  \"graph\": \"" << (result.graph_name) << "\",\n";
        out << "  \"nodes\": \"" << (result.nodes) << "\",\n";
        out << "  \"edges\": \"" << (result.edges) << "\",\n";
        out << "  \"status\": " << static_cast<int>(result.status) << ",\n";

        if (!result.routing_base.empty()) {
            out << "  \"routing_base\": \""
                << jsonEscape(result.routing_base) << "\",\n";
        }

        out << "  \"total_runtime_microseconds\": "
            << result.total_runtime_microseconds << ",\n";
        out << "  \"preprocessing_runtime_microseconds\": "
            << result.preprocessing_runtime_microseconds << ",\n";
        out << "  \"solve_runtime_microseconds\": "
            << result.solve_runtime_microseconds << ",\n";
        out << "  \"oblivious_ratio\": "
            << result.oblivious_ratio << ",\n";
        out << "  \"candidate_paths\": "
            << result.candidate_paths << ",\n";
        out << "  \"average_paths_per_pair\": "
            << result.average_paths_per_pair << ",\n";

        if (!result.mwu_metrics.empty()) {
            const auto& mwu = result.mwu_metrics;

            out << "  \"mwu_metrics\": {\n";
            out << "    \"iteration_count\": " << mwu.iteration_count << ",\n";
            out << "    \"solve_time_microseconds\": " << mwu.solve_time << ",\n";
            out << "    \"transformation_time_microseconds\": "
                << mwu.transformation_time << ",\n";
            out << "    \"load_computation_time_microseconds\": "
                << mwu.load_computation_time << ",\n";
            out << "    \"weight_update_time_microseconds\": "
                << mwu.mwu_weight_update_time << ",\n";
            out << "    \"average_oracle_time_microseconds\": "
                << mwu.averageOracleTime() << "\n";

            out << "  },\n";
        }

        if (!result.expander_metrics.empty()) {
            writeExpanderMetricsJson(
                result.expander_metrics,
                out
            );
        }


        out << "  \"demand_evaluations\": [\n";

        for (std::size_t i = 0; i < result.demand_evaluations.size(); ++i) {
            const auto& eval = result.demand_evaluations[i];

            out << "    {\n";
            out << "      \"demand_model\": \""
                << jsonEscape(demandModelName(eval.demand_type)) << "\",\n";
            out << "      \"congestion\": "
                << eval.congestion << ",\n";
            out << "      \"runtime_microseconds\": "
                << eval.runtime_microseconds << "\n";
            out << "    }";

            if (i + 1 < result.demand_evaluations.size()) {
                out << ",";
            }

            out << "\n";
        }

        out << "  ]\n";


        out << "}\n";


    }

    static void writeExpanderMetricsJson(
    const ExpanderMetrics& metrics,
    std::ostream& out
) {
    out << "  \"expander_metrics\": {\n";

    out << "    \"hierarchy_runtime_microseconds\": "
        << metrics.hierarchy_runtime_microseconds
        << ",\n";

    out << "    \"tree_runtime_microseconds\": "
        << metrics.tree_runtime_microseconds
        << ",\n";

    out << "    \"basis_flow_runtime_microseconds\": "
        << metrics.basis_flow_runtime_microseconds
        << ",\n";

    out << "    \"hierarchy_levels\": "
        << metrics.hierarchy_levels
        << ",\n";

    out << "    \"hierarchy_clusters\": "
        << metrics.hierarchy_clusters
        << ",\n";

    out << "    \"clusters_per_level\": [";

    for (std::size_t index = 0;index < metrics.clusters_per_level.size();++index) {
        if (index > 0) {
            out << ", ";
        }

        out << metrics.clusters_per_level[index];
    }

    out << "],\n";

    out << "    \"max_cluster_vertices\": "
        << metrics.max_cluster_vertices
        << ",\n";

    out << "    \"average_cluster_vertices\": "
        << metrics.average_cluster_vertices
        << ",\n";

    out << "    \"tree_nodes\": "
        << metrics.tree_nodes
        << ",\n";

    out << "    \"tree_edges\": "
        << metrics.tree_edges
        << ",\n";

    out << "    \"tree_depth\": "
        << metrics.tree_depth
        << ",\n";

    out << "    \"basis_flows\": "
        << metrics.basis_flows
        << ",\n";

    out << "    \"total_electrical_solves\": "
        << metrics.total_electrical_solves
        << ",\n";

    out << "    \"average_electrical_solves\": "
        << metrics.averageElectricalSolves()
        << ",\n";

    out << "    \"max_basis_embedding_congestion\": "
        << metrics.max_basis_embedding_congestion
        << ",\n";

    out << "    \"max_conservation_error\": "
        << metrics.max_conservation_error
        << "\n";

    out << "  },\n";
}
};

#endif //OBLIVIOUSROUTING_RESULT_IO_H