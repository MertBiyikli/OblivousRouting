//
// Created by Mert Biyikli on 17.06.26.
//

#ifndef OBLIVIOUSROUTING_UTILS_H
#define OBLIVIOUSROUTING_UTILS_H

#pragma once

#include <fstream>
#include "../routing/routing_result.h"


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
        out << "Date: " << std::chrono::system_clock::now() << '\n';
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
        out << "  \"date\": \"" << (std::chrono::system_clock::now()) << "\",\n";
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


class MWUMetrics{
public:
    MWUMetrics() {
        iteration_count = 0;
        solve_time = 0;
        transformation_time = 0;
        mwu_weight_update_time = 0;
    }

    virtual ~MWUMetrics() = default;
    std::vector<double> oracle_running_times;
    int iteration_count;
    double solve_time;
    double transformation_time;
    double mwu_weight_update_time;
    double load_computation_time{};


    void printTimeStats() {
        std::cout << "Solve time: " << this->solve_time << " micro seconds\n";
        std::cout << "Transformation time: " << transformation_time << " micro seconds\n";
        std::cout << "MWU iterations: " << this->iteration_count << "\n";
        std::cout << "MWU load computation: " << this->load_computation_time << " micro seconds\n";
        double average_oracle_time = 0.0;
        for (double t : this->oracle_running_times) {
            average_oracle_time += t;
        }
        std::cout << "Average oracle time: " << (average_oracle_time/static_cast<double>(this->iteration_count)) << " micro seconds\n";
        std::cout << "Total MWU weight update time: " << mwu_weight_update_time << " micro seconds\n";
    }

    const int getIterationCount() const {
        return iteration_count;
    }
};

#endif //OBLIVIOUSROUTING_UTILS_H