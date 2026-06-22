//
// Created by Mert Biyikli on 17.06.26.
//

#ifndef OBLIVIOUSROUTING_UTILS_H
#define OBLIVIOUSROUTING_UTILS_H

#include <fstream>
#include "../io/solver_io.h"
#include "../io/demand_io.h"
enum class OutputFormat {
    TEXT,
    JSON
};

enum class ResultStatus {
    OK,
    ERROR_INVALID_SOLVER,
    ERROR_INVALID_ROUTING_SCHEME,
    ERROR_MISSING_DEMAND_MODELS,
    ERROR_EXCEPTION
};

struct DemandEvaluationResult {
    DemandModelType demand_type{};
    double congestion = -1.0;

    double runtime_microseconds = -1.0;
};

struct RoutingRunResult {
    ResultStatus status = ResultStatus::OK;
    SolverType type{};

    std::string solver_name;
    std::string routing_base; // "Tree", "Electrical", or empty for normal solvers

    double total_runtime_microseconds = 0.0;
    double preprocessing_runtime_microseconds = 0.0;
    double solve_runtime_microseconds = 0.0;

    double oblivious_ratio = -1.0;
    int mwu_iterations = -1;

    std::size_t candidate_paths = 0;
    double average_paths_per_pair = 0.0;

    std::vector<DemandEvaluationResult> demand_evaluations;

    // Optional: only filled for normal solvers or last scheme if needed
    std::unique_ptr<RoutingScheme> scheme;

    // Semi-oblivious: one demand-specific scheme per demand model
    std::unordered_map<std::string, std::unique_ptr<RoutingScheme>> demand_schemes;
};



class RoutingResultWriter {
public:
    static bool write(
        const RoutingRunResult& result,
        const std::string& path,
        OutputFormat format
    ) {
        std::ofstream file(path);

        if (!file.is_open()) {
            std::cerr << "[ERROR] Failed to open output file: " << path << '\n';
            return false;
        }

        switch (format) {
            case OutputFormat::TEXT:
                writeText(result, file);
                return true;

            case OutputFormat::JSON:
                writeJson(result, file);
                return true;
            default:
                std::cerr << "[ERROR] Unsupported output format\n";
                return false;
        }
    }

private:
    static void writeText(
        const RoutingRunResult& result,
        std::ostream& out
    ) {
        out << "Solver: " << result.solver_name << '\n';

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
        out << "MWU iterations: " << result.mwu_iterations << '\n';

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

    static void writeJson(
        const RoutingRunResult& result,
        std::ostream& out
    ) {
        out << "{\n";
        out << "  \"solver\": \"" << (result.solver_name) << "\",\n";
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
        out << "  \"mwu_iterations\": "
            << result.mwu_iterations << ",\n";
        out << "  \"candidate_paths\": "
            << result.candidate_paths << ",\n";
        out << "  \"average_paths_per_pair\": "
            << result.average_paths_per_pair << ",\n";

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
};
#endif //OBLIVIOUSROUTING_UTILS_H