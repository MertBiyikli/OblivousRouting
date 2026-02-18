#include <iostream>
#include "src/datastructures/IGraph.h"
#include "src/datastructures/GraphCSR.h"
#include "src/datastructures/GraphADJ.h"
#include "src/parse_parameter.h"


int main(int argc, char **argv) {

    std::string err;
    auto cfg = parse_parameter(argc, argv, &err);
    if (!cfg) { std::cerr << err; return -1; }

#ifdef _OPENMP
    std::cout << "Number of threads active: " << omp_get_max_threads() << "\n";
#endif

    GraphCSR G;
    readLGFFile(G, cfg->filename,  true);
    G.finalize();
    std::cout << "Graph loaded: " << G.getNumNodes() << " nodes, "
              << G.getNumEdges() << " edges.\n";

    if ( !G.isConnected()) {
        std::cerr << "Input graph is not connected!\n";
        // return -1;
    }


    for (SolverType type : cfg->solvers) {
        std::cout << "\n=== Running solver: " << ::solverNames[static_cast<int>(type)] << " ===\n";

        auto solver = makeSolver(type, G);

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

        // --- optional: demand model evaluation ---
        auto result = HandleDemandModel(argc, argv, cfg, G, scheme);
        printStatsForDemandModel(argv, result);
    }
    return 0;
}
