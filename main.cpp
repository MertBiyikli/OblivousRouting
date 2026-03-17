#include <iostream>
#include "src/datastructures/IGraph.h"
#include "src/datastructures/GraphCSR.h"
#include "src/datastructures/GraphADJList.h"
#include "src/datastructures/GraphADJMatrix.h"
#include "src/parse_parameter.h"


int main(int argc, char **argv) {

    std::string err;
    auto cfg = parse_parameter(argc, argv, &err);
    if (!cfg) { std::cerr << err; return -1; }


    auto G = makegraph(cfg->graph_format);
    readLGFFile(*G, cfg->filename,  true);
    G->finalize();

    std::cout << "Graph loaded: " << G->getNumNodes() << " nodes, " << G->getNumEdges() << " edges.\n";
    G->printGraphType();

    // Pre-compute demand maps and offline optimal congestions once — outside the solver loop.
    std::unordered_map<std::string, DemandMap> demand_maps;
    std::unordered_map<std::string, double>    offline_opt_per_model;
    HandleDemandModels(cfg, *G, [&](const std::string& model_name, const DemandMap& dmap) {
        demand_maps[model_name]         = dmap;
        offline_opt_per_model[model_name] = computeOfflineOptimalCongestion(*G, dmap);
    });

    for (SolverType type : cfg->solvers) {
        std::cout << "\n=== Running solver: " << ::solverNames[static_cast<int>(type)] << " ===\n";

        auto solver = makeSolver(type, *G);

        auto t0 = timeNow();
        if (cfg->num_threads > 1) {
            solver->setNumThreads(cfg->num_threads);
        }
        auto scheme = solver->solve();
        double solve_time = duration(timeNow() - t0);

        //scheme->printRoutingTable();
        std::cout << "Total running time: " << solve_time << " micro seconds\n";

        if (auto mwu = dynamic_cast<MWUFramework*>(solver.get())) {
            mwu->printTimeStats();
        }

        double total_scales = 0;
        if (auto tree_solver_flat = dynamic_cast<TreeMWU<FlatHST>*>(solver.get())) {
            for (int& scales : tree_solver_flat->getScales()) {
                total_scales += scales;
            }
            std::cout << "Average scales: " << total_scales/tree_solver_flat->getScales().size() << "\n";
        }else if (auto tree_solver_pointer = dynamic_cast<TreeMWU<std::shared_ptr<HSTNode>>*>(solver.get())) {
            for (int& scales : tree_solver_pointer->getScales()) {
                total_scales += scales;
            }
            std::cout << "Average scales: " << total_scales/tree_solver_pointer->getScales().size() << "\n";
        }

        if (scheme) {
            if ( auto linear_scheme = dynamic_cast<LinearRoutingScheme*>(scheme.get()) ) {
                double oblivious_ratio = linear_scheme->computeObliviousRatio();
                std::cout << "Oblivious ratio of the linear routing scheme: " << oblivious_ratio << "\n";
            }
        }

        if (cfg->demand_models.empty()) continue;

        for (const auto& [model_name, offline_cong] : offline_opt_per_model) {
            const DemandMap& dmap = demand_maps.at(model_name);
            double scheme_cong = computeRoutingSchemeCongestion(*G, scheme, dmap);
            printStatsForDemandModel(model_name, {offline_cong, scheme_cong});
        }
    }
    return 0;
}
