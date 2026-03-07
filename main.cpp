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
        auto scheme = solver->solve();
        double solve_time = duration(timeNow() - t0);

        // scheme->printRoutingTable();
        std::cout << "Total running time: " << solve_time << " ms\n";

        if (auto mwu = dynamic_cast<MWUFramework*>(solver.get())) {
            mwu->printTimeStats();
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
