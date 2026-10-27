//
// Created by Mert on 27.10.26.
//

#ifndef E_ROUTING_ANALYSIS_ENGINE_H
#define E_ROUTING_ANALYSIS_ENGINE_H


#include <analysis/analysis_result.h>
#include <analysis/scenario_analysis.h>

#include <filesystem>
#include <fstream>
#include <iomanip>

class ScenarioEngine {
public:
    static Result<ScenarioResult> runScenario(optimized::Graph<EdgeData>& graph, Scenario& scenario, const Config& cfg,ExecutionMode mode, double demand_scale = 1.0);
    static Result<std::pair<demands, double>> computeTotalDemandVolume(const optimized::Graph<EdgeData>& graph, const Config& cfg);
};



inline auto printBottlenecks = [](const NetworkData& data,
                           const std::string& name) {
    auto links = data.links;

    std::sort(
        links.begin(),
        links.end(),
        [](const LinkData& a, const LinkData& b) {
            return a.utilization > b.utilization;
        }
    );

    std::cout << "\n" << name << " — Top bottlenecks\n";

    for (std::size_t i = 0;
         i < std::min<std::size_t>(5, links.size());
         ++i) {

        const auto& link = links[i];

        std::cout
            << "Edge " << link.edge_id
            << " | Load: " << link.load
            << " | Capacity: " << link.capacity
            << " | Utilization: " << link.utilization
            << '\n';
         }
};



static Result<void> printScenarioAnalysis(const ScenarioResult& result) {
    for (const auto& solver : result.solver_results) {

        if (!solver.baseline_Data || !solver.fixed_routing_Data || !solver.modified_Data || !solver.fixed_routing_comparison || !solver.comparison) {
            return makeErrorMessage(ErrorCode::SolverFailed, "Incomplete solver analysis.");
        }

        const auto& baseline = *solver.baseline_Data;

        const auto& fixed = *solver.fixed_routing_Data;

        const auto& reoptimized = *solver.modified_Data;


        std::cout << "Solver algorithms: " << std::endl;
        for (const auto& algo_names : solver.solver_algorithms) {
            std::cout << "Algorithm: " << algo_names << std::endl;
        }

        std::cout << "Baseline congestion: " << baseline.max_congestion << '\n';

        std::cout << "Fixed routing congestion: " << fixed.max_congestion << '\n';

        std::cout << "Reoptimized congestion: " << reoptimized.max_congestion << '\n';

        const auto fixed_reduction = solver.fixed_routing_comparison->relative_reduction;

        const auto reoptimized_reduction = solver.comparison->relative_reduction;

        if (fixed_reduction) {
            std::cout << "Fixed routing improvement: " << *fixed_reduction * 100.0 << "%\n";
        }

        if (reoptimized_reduction) {
            std::cout << "Reoptimized improvement: " << *reoptimized_reduction * 100.0 << "%\n";
        }

        std::cout << "Baseline solver time: " << solver.baseline_solver_seconds << " s\n";

        std::cout << "Modified solver time: " << solver.modified_solver_seconds << " s\n";

        printBottlenecks(*solver.baseline_Data, "Baseline");
        printBottlenecks(*solver.fixed_routing_Data, "Fixed routing");
        printBottlenecks(*solver.modified_Data, "Reoptimized");


        double max_load_difference = 0.0;
        double total_load_difference = 0.0;
        std::size_t changed_edges = 0;

        std::unordered_map<int, double> baseline_loads;

        for (const auto& link : baseline.links) {
            baseline_loads.emplace(link.edge_id, link.load);
        }

        for (const auto& link : reoptimized.links) {
            const auto it = baseline_loads.find(link.edge_id);

            if (it == baseline_loads.end()) {
                continue;
            }

            const double difference =
                std::abs(it->second - link.load);

            max_load_difference =
                std::max(max_load_difference, difference);

            total_load_difference += difference;

            if (difference > 1e-8) {
                ++changed_edges;
            }
        }

        std::cout << "\nRouting redistribution analysis:\n";
        std::cout << "Changed physical-link loads: "
                  << changed_edges << '\n';
        std::cout << "Maximum load difference: "
                  << max_load_difference << '\n';
        std::cout << "Total absolute load difference: "
                  << total_load_difference << '\n';


        for (const auto& link : reoptimized.links) {
            auto it = baseline_loads.find(link.edge_id);
            if (it == baseline_loads.end()) continue;

            const double delta = link.load - it->second;

            if (std::abs(delta) > 1e-8) {
                std::cout
                    << "Edge " << link.edge_id
                    << " | Before: " << it->second
                    << " | After: " << link.load
                    << " | Delta: " << delta
                    << '\n';
            }
        }

    }
    return {};
}

// Export the exact scenario demands and one entry per physical link.
// Assumes physical links have symmetric capacity and reverse-arc pairing.
class ScenarioLPExporter {
public:
    static Result<void> save(
        const optimized::Graph<EdgeData>& graph,
        const demands& fixed_demands,
        const std::filesystem::path& output_path
    ) {
        namespace fs = std::filesystem;
        std::error_code ec;
        if (!output_path.parent_path().empty()) {
            fs::create_directories(output_path.parent_path(), ec);
            if (ec) {
                return makeErrorMessage(ErrorCode::RuntimeError,
                    "Cannot create LP export directory: " + ec.message());
            }
        }

        std::ofstream out(output_path, std::ios::out | std::ios::trunc);
        if (!out) {
            return makeErrorMessage(ErrorCode::RuntimeError,
                "Cannot open LP export: " + output_path.string());
        }

        out << std::setprecision(17);
        out << "{\n  \"nodes\": " << graph.getNumNodes() << ",\n  \"edges\": [\n";

        bool first = true;
        const int m = graph.getNumDirectedEdges();
        for (int e = 0; e < m; ++e) {
            const int reverse = graph.reverse(e).id;
            if (reverse < 0 || reverse >= m || graph.reverse(reverse).id != e) {
                return makeErrorMessage(ErrorCode::InvalidArgument,
                    "Invalid reverse-edge pairing in LP export.");
            }
            if (reverse == e) {
                return makeErrorMessage(ErrorCode::InvalidArgument,
                    "Self-reversing arc not supported by physical-link LP export.");
            }
            if (e > reverse) continue;

            const auto [u, v] = graph.getEdgeEndpoints(e);
            const double c = graph.edgeData(e).capacity;
            const double reverse_c = graph.edgeData(reverse).capacity;
            if (u == v || !std::isfinite(c) || c <= 0.0 ||
                !std::isfinite(reverse_c) || reverse_c <= 0.0 ||
                std::abs(c - reverse_c) > 1e-9 * std::max({1.0, c, reverse_c})) {
                return makeErrorMessage(ErrorCode::InvalidArgument,
                    "LP exporter needs positive, symmetric capacities for each physical link.");
            }
            if (!first) out << ",\n";
            first = false;
            out << "    {\"id\": " << e << ", \"u\": " << u
                << ", \"v\": " << v << ", \"capacity\": " << c << "}";
        }
        out << "\n  ],\n  \"demands\": [\n";

        for (std::size_t i = 0; i < fixed_demands.demand_values.size(); ++i) {
            const auto [s, t] = fixed_demands.getDemandPair(i);
            const double d = fixed_demands.getDemandValue(i);
            if (s < 0 || s >= graph.getNumNodes() ||
                t < 0 || t >= graph.getNumNodes() || s == t ||
                !std::isfinite(d) || d < 0.0) {
                return makeErrorMessage(ErrorCode::InvalidDemand,
                    "Invalid commodity in LP export.");
            }
            if (i) out << ",\n";
            out << "    {\"source\": " << s << ", \"target\": " << t
                << ", \"value\": " << d << "}";
        }
        out << "\n  ]\n}\n";
        if (!out) {
            return makeErrorMessage(ErrorCode::RuntimeError,
                "Failed writing LP export: " + output_path.string());
        }
        return {};
    }
};

#endif // E_ROUTING_ANALYSIS_ENGINE_H
