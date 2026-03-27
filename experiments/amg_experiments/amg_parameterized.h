//
// Created by Mert Biyikli on 09.02.26.
//

#ifndef OBLIVIOUSROUTING_AMG_PARAMETERIZED_H
#define OBLIVIOUSROUTING_AMG_PARAMETERIZED_H

//
// Created by Mert Biyikli on 09.02.26.
//
#include <iostream>
#include <fstream>
#include "../../src/datastructures/IGraph.h"
#include "../../src/parse_parameter.h"

// ---- small helpers ----
static inline std::string csvEscape(std::string s) {
    // wrap in quotes if needed; escape internal quotes
    bool need = false;
    for (char c : s) if (c == ',' || c == '"' || c == '\n' || c == '\r') { need = true; break; }
    if (!need) return s;

    std::string out;
    out.reserve(s.size() + 2);
    out.push_back('"');
    for (char c : s) {
        if (c == '"') out.push_back('"');
        out.push_back(c);
    }
    out.push_back('"');
    return out;
}

template <class T>
static inline T ptGet(const boost::property_tree::ptree& p, const std::string& key, const T& def) {
    auto opt = p.get_optional<T>(key);
    return opt ? *opt : def;
}

static inline double msSince(const std::chrono::steady_clock::time_point& t0) {
    using namespace std::chrono;
    return duration_cast<duration<double, std::milli>>(steady_clock::now() - t0).count();
}

class ElectricalFlowWithManyParams {
public:
    void init(const std::string& config_file, IGraph& g, const std::string& csv_file = "amg_parameterized_results.csv") {

        // Prepare CSV (overwrite)
        std::ofstream out(csv_file);
        if (!out) {
            // Create a new csv file
            out.open(csv_file, std::ios::out);
            if (!out) {
                std::cerr << "Error: Unable to create CSV file: " << csv_file << "\n";
                return;
            }
        }

        // Header
        out
            << "graph_nodes,graph_edges,config_file,run_name,"
            << "solver_type,coarsening_type,relax_type,solver_tol,solver_maxiter,"
            << "init_ms,run_ms,scale_ms,total_ms\n";

        // Load parameters from config file
        boost::property_tree::ptree prm;
        boost::property_tree::read_json(config_file, prm);
        LinearRoutingTable table;
        table.init(g);

        for (auto& run : prm.get_child("runs")) {
            std::string run_name = run.second.get<std::string>("name");
            auto solver_params = run.second.get_child("params");

            // Extract a few common params for CSV (safe defaults if missing)
            const std::string solver_type   = ptGet<std::string>(solver_params, "solver.type", "<unset>");
            const std::string coarsen_type  = ptGet<std::string>(solver_params, "precond.coarsening.type", "<unset>");
            const std::string relax_type    = ptGet<std::string>(solver_params, "precond.relax.type", "<unset>");
            const double solver_tol         = ptGet<double>(solver_params, "solver.tol", -1.0);
            const int solver_maxiter        = ptGet<int>(solver_params, "solver.maxiter", -1);

            std::cout << "Initializing solver: " << run_name << "\n";

            // Measure timings
            double init_ms = 0.0, run_ms = 0.0, scale_ms = 0.0;

            auto t_total = std::chrono::steady_clock::now();

            {
                auto t = std::chrono::steady_clock::now();
                std::unique_ptr<ElectricalFlowOptimized> solver =
                    std::make_unique<ElectricalFlowOptimized>(g, 0);

                // your init signature: init(bool, string, ptree)
                solver->init(false, solver_params);
                init_ms = msSince(t);

                t = std::chrono::steady_clock::now();
                solver->run(table);
                run_ms = msSince(t);

                t = std::chrono::steady_clock::now();
                solver->scaleFlowDown(table);
                scale_ms = msSince(t);

                // keep your console stats (optional)
                solver->printTimeStats();
            }

            const double total_ms = msSince(t_total);

            // Write CSV row
            out
                << g.getNumNodes() << "," << g.getNumDirectedEdges() << ","
                << csvEscape(config_file) << ","
                << csvEscape(run_name) << ","
                << csvEscape(solver_type) << ","
                << csvEscape(coarsen_type) << ","
                << csvEscape(relax_type) << ","
                << std::setprecision(17) << solver_tol << ","
                << solver_maxiter << ","
                << init_ms << ","
                << run_ms << ","
                << scale_ms << ","
                << total_ms
                << "\n";

            out.flush();

            // reset table for next run
            table.init(g);
        }

        std::cout << "AMG experiment CSV written to: " << csv_file << "\n";
    }
};

int AMG_Experiment(IGraph& G, const std::string csv) {
    ElectricalFlowWithManyParams ef_experiments;
    ef_experiments.init("../configs/amg_configs_multiple_params.json", G, csv);

    return 0;
}

#endif //OBLIVIOUSROUTING_AMG_PARAMETERIZED_H