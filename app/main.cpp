#include "../include/data_structures/graph/graph_csr.h"
#include "../include/algorithms/lp/lp_ac.h"
#include "../include/algorithms/mwu/electrical_mwu.h"
#include "../include/algorithms/mwu/tree_mwu.h"
#include "../include/algorithms/mwu/oracle/tree/tree_oracle.h"
#include "../include/algorithms/mwu/oracle/tree/mst/mst_oracle.h"
#include "../include/algorithms/mwu/oracle/tree/frt/frt.h"
#include "../include/algorithms/mwu/oracle/tree/fast_ckr/fast_ckr.h"

#include "../include/io/parse_argurment_io.h"
#include "../include/io/solver_io.h"
#include "../include/io/demand_io.h"
#include "../include/io/graph_io.h"
#include "../include/core/routing_scheme.h"


// Helper to get solver name from type
inline std::string getSolverName(SolverType type) {
    static const std::map<SolverType, std::string> names{
        {SolverType::ELECTRICAL, "Electrical Flow"},
        {SolverType::RAECKE_FRT_FLAT, "Raecke FRT (Flat HST)"},
        {SolverType::RAECKE_CKR_FLAT, "Raecke CKR (Flat HST)"},
        {SolverType::RAECKE_RANDOM_MST_FLAT, "Random MST (Flat HST)"},
        {SolverType::RAECKE_FRT_MENDELSCALING_FLAT, "Raecke FRT + MendelScaling (Flat HST)"},
        {SolverType::RAECKE_CKR_MENDELSCALING_FLAT, "Raecke CKR + MendelScaling (Flat HST)"},
        {SolverType::LP_APPLEGATE_COHEN, "LP Applegate-Cohen"},
        {SolverType::ELECTRICAL_PARALLEL_BATCHES, "Electrical Flow (Parallel)"},
        {SolverType::RAECKE_FRT_POINTER, "Raecke FRT (Pointer HST)"},
        {SolverType::RAECKE_CKR_POINTER, "Raecke CKR (Pointer HST)"},
        {SolverType::RAECKE_RANDOM_MST_POINTER, "Random MST (Pointer HST)"},
        {SolverType::RAECKE_FRT_MENDELSCALING_POINTER, "Raecke FRT + MendelScaling (Pointer HST)"},
        {SolverType::RAECKE_CKR_MENDELSCALING_POINTER, "Raecke CKR + MendelScaling (Pointer HST)"}
    };
    auto it = names.find(type);
    return (it != names.end()) ? it->second : "Unknown Solver";
}

int main(int argc, char **argv) {
    std::string err;
    auto cfg = parse_parameter(argc, argv, &err);
    if (!cfg) {
        std::cerr << err << std::endl;
        return -1;
    }

    // Load or create graph
    auto graph = makegraph(cfg->graph_format);
    if (!cfg->filename.empty()) {
        readLGFFile(*graph, cfg->filename);
    }
    graph->finalize();
    std::cout << "Graph loaded: " << graph->getNumNodes() << " nodes, " << graph->getNumUndirectedEdges() <<" edges.\n";
    // Precompute offline optimal congestion if needed
    std::map<std::string, double> offline_opt_per_model;
    std::map<std::string, demands> demand_maps;

    if (!cfg->demand_models.empty()) {
        HandleDemandModels(cfg, *graph,
            [&](const std::string& model_name, const demands& dmap) {
                demand_maps[model_name] = dmap;
                offline_opt_per_model[model_name] = computeOfflineOptimalCongestion(*graph, dmap);
            });
    }

    // Run each solver
    for (SolverType type : cfg->solvers) {
        std::cout << "\n=== Running solver: " << getSolverName(type) << " ===\n";

        auto solver_opt = makeSolver(type, *graph, cfg->cycle_strategy);
        if (!solver_opt) {
            std::cerr << "Failed to create solver of type " << static_cast<int>(type) << "\n";
            continue;
        }
        auto& solver = *solver_opt;

        auto t0 = timeNow();
        auto scheme = solver->solve();
        auto t1 = timeNow();
        double total_time = duration(t1-t0);

        std::cout << "Total time: " << total_time << " micro seconds\n";

        //scheme->printRoutingTable();

        // Print time statistics if available
        if (auto mwu = dynamic_cast<MWUFramework*>(solver.get())) {
            mwu->printTimeStats();
        }

        // Compute oblivious ratio for linear schemes
        if (auto linear_scheme = dynamic_cast<LinearRoutingScheme*>(scheme.get())) {
            double oblivious_ratio = linear_scheme->computeObliviousRatio();
            std::cout << "Oblivious ratio: " << oblivious_ratio << "\n";
        }

        // Evaluate demand models if provided
        if (cfg->demand_models.empty()) continue;

        for (const auto& [model_name, offline_cong] : offline_opt_per_model) {
            const demands& dmap = demand_maps.at(model_name);
            double scheme_cong = computeRoutingSchemeCongestion(*graph, scheme, dmap);
            printStatsForDemandModel(model_name, {offline_cong, scheme_cong});
        }
    }

    return 0;
}