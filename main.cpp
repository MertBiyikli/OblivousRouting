#include <iostream>
#include "src/datastructures/GraphADJ.h"
#include "src/electrical/electrical_flow_optimized.h"
#include "src/datastructures/GraphCSR.h"
#include "src/parse_parameter.h"
#include "src/expanders/xcut_expander_hierarchy.h"
#include "src/tree_based/tree_mwu.h"


#include "experiments/amg_experiments/amg_parameterized.h"
#include "src/tree_based/fast_ckr.h"

#include "src/tree_based/frt.h"
#include "src/tree_based/random_mst.h"


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
    //G.print();

    //ExpanderRouting expander_routing;
    // expander_routing.init(G);

    // Run the AMG parameterized experiment

    // first define the output.csv
/*
    if (argv[2] == nullptr) {
        std::cerr << "Graph input is missing.\n";
        return -1;
    }
    std::string graphFile = argv[2];
    auto it = graphFile.find(".lgf");
    // remove the extension
    if (it != std::string::npos) {
        graphFile = graphFile.substr(0, it);
    }
    std::string output = "../results/amg_parameterized_results/"+graphFile.substr(graphFile.find_last_of("/\\") + 1)+".csv";
    AMG_Experiment(G, output);
*/

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
        if (result.first > 0.0 && result.second > 0.0) {
            std::cout << "Ratio off the optimal offline solution "
                      << (argv[3] ? argv[3] : "<empty>")
                      << " model demand: "
                      << ( result.second  / result.first ) * 100.0 << "% "
                      << " ("
                      << result.first << " / " << result.second
                      << ")\n";
        }
    }
    return 0;
}
