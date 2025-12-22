#include <iostream>
#include "src/datastructures/GraphADJ.h"
#include "src/electrical/electrical_flow_optimized.h"
#include "src/datastructures/GraphCSR.h"
#include "src/parse_parameter.h"



#define duration(a) std::chrono::duration_cast<std::chrono::milliseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

int main(int argc, char **argv) {

    std::string err;
    auto cfg = parse_parameter(argc, argv, &err);
    if (!cfg) { std::cerr << err; return -1; } // EX_USAGE


    GraphCSR g_csr;// = std::make_unique<GraphCSR>();
    readLGFFile(g_csr, cfg->filename,  true);
    g_csr.finalize();
    std::cout << "Graph loaded: " << g_csr.getNumNodes() << " nodes, "
              << g_csr.getNumEdges() << " edges.\n";

    for (SolverType type : cfg->solvers) {
        std::cout << "\n=== Running solver: " << ::solverNames[static_cast<int>(type)] << " ===\n";

        auto solver = makeSolver(type, g_csr);

        // --- solve ---
        auto t0 = timeNow();
        auto scheme = solver->solve();
        double solve_time = duration(timeNow() - t0);

        std::cout << "Solve time: " << solve_time << " ms\n";

        // --- optional: demand model evaluation ---
        auto result = HandleDemandModel(argc, argv, cfg, g_csr, scheme);
        if (result.first > 0.0 && result.second > 0.0) {
            std::cout << "Ratio off the optimal offline solution "
                      << (argv[3] ? argv[3] : "<empty>")
                      << " model demand: "
                      << ( result.second - result.first  / result.second ) * 100.0 << "% "
                      << " ("
                      << result.first << " / " << result.second
                      << ")\n";
        }

    }

    return 0;
}
