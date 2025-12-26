#include <iostream>
#include "src/datastructures/GraphADJ.h"
#include "src/electrical/electrical_flow_optimized.h"
#include "src/datastructures/GraphCSR.h"
#include "src/parse_parameter.h"





int main(int argc, char **argv) {

    std::string err;
    auto cfg = parse_parameter(argc, argv, &err);
    if (!cfg) { std::cerr << err; return -1; }


    GraphCSR G;
    readLGFFile(G, cfg->filename,  true);
    G.finalize();
    std::cout << "Graph loaded: " << G.getNumNodes() << " nodes, "
              << G.getNumEdges() << " edges.\n";

    for (SolverType type : cfg->solvers) {
        std::cout << "\n=== Running solver: " << ::solverNames[static_cast<int>(type)] << " ===\n";

        auto solver = makeSolver(type, G);

        // --- solve ---
        auto t0 = timeNow();
        auto scheme = solver->solve();
        double solve_time = duration(timeNow() - t0);

        // scheme->printRoutingTable();

        std::cout << "Total running time: " << solve_time << " ms\n";

        // if the algorithm is an MWU framework, print iteration stats
        if (auto mwu = dynamic_cast<MWUFramework*>(solver.get())) {
            std::cout << "Solve time: " << mwu->solve_time << " ms\n";
            std::cout << "Transformation time: " << (solve_time - mwu->solve_time) << " ms\n";
            std::cout << "MWU iterations: " << mwu->iteration_count << "\n";
            double average_oracle_time = 0.0;
            for (double t : mwu->oracle_running_times) {
                average_oracle_time += t;
            }
            std::cout << "Average oracle time: " << (average_oracle_time/static_cast<double>(mwu->oracle_running_times.size())) << " ms\n";
        }

        // --- optional: demand model evaluation ---
        auto result = HandleDemandModel(argc, argv, cfg, G, scheme);
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
