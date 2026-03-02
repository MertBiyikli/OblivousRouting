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

    // --- optional: demand model evaluation ---
    DemandMap demand_map;
    double offline_opt_cong = 0;
    HandleDemandModel(argc, argv, cfg, *G, demand_map);
    if (cfg->demand_model != DemandModelType::NONE) {
        offline_opt_cong = computeOfflineOptimalCongestion(*G, demand_map);
    }

    for (SolverType type : cfg->solvers) {
        std::cout << "\n=== Running solver: " << ::solverNames[static_cast<int>(type)] << " ===\n";

        auto solver = makeSolver(type, *G);

        // --- solve ---
        auto t0 = timeNow();
        auto scheme = solver->solve();
        double solve_time = duration(timeNow() - t0);

        //scheme->printRoutingTable();

        std::cout << "Total running time: " << solve_time << " ms\n";

        // if the algorithm is an MWU framework, print iteration stats
        if (auto mwu = dynamic_cast<MWUFramework*>(solver.get())) {
            mwu->printTimeStats();
        }

        // compute and print the congestion of the computed routing scheme
        double scheme_congestion = computeRoutingSchemeCongestion(*G, scheme, demand_map);
        if (cfg->demand_model != DemandModelType::NONE) {
            printStatsForDemandModel(argv, {offline_opt_cong, scheme_congestion});
        }
    }
    return 0;
}
