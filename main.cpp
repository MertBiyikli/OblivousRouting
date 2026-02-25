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
    G->printGraphType();
    readLGFFile(*G, cfg->filename,  true);
    G->finalize();
/*
    for (int e = 0; e<G->getNumEdges(); e++) {
        auto [u, v] = G->edgeEndpoints(e);
        if (G->getEdgeDistance(e) <= 0.0 || std::isnan(G->getEdgeDistance(e))) {
            std::cerr << "Warning: edge (" << u << ", " << v << ") has non-positive or NaN distance. Setting to 1.0.\n";
            G->updateEdgeDistance(e, 1.0);
        }else {
            std::cout << "Edge (" << u << ", " << v << ") has distance: " << G->getEdgeDistance(e) << "\n";
        }
         if (G->getEdgeCapacity(e) <= 0.0 || std::isnan(G->getEdgeCapacity(e))) {
             std::cerr << "Warning: edge (" << u << ", " << v << ") has non-positive or NaN capacity. Setting to 1.0.\n";
         }else {
             std::cout << "Edge (" << u << ", " << v << ") has capacity: " << G->getEdgeCapacity(e) << "\n";
         }
    }
    std::cout << "Graph loaded: " << G->getNumNodes() << " nodes, " << G->getNumEdges() << " edges.\n";

    // compute shortest path distances and diameter for later use in some algorithms
    std::cout << "Computing shortest path distances and diameter...\n";
    for (auto c : G->getShortestPath(0, 3)) {
        std::cout << c << " ";
    }
    std::cout << "\nDiameter: " << G->GetDiameter() << "\n";
*/
    if ( !G->isConnected()) {
        std::cerr << "Input graph is not connected!\n";
    }

    for (SolverType type : cfg->solvers) {
        std::cout << "\n=== Running solver: " << ::solverNames[static_cast<int>(type)] << " ===\n";

        auto solver = makeSolver(type, *G);

        // --- solve ---
        auto t0 = timeNow();
        auto scheme = solver->solve();
        double solve_time = duration(timeNow() - t0);

        // scheme->printRoutingTable();


        std::cout << "Total running time: " << solve_time << " ms\n";

        // if the algorithm is an MWU framework, print iteration stats
        if (auto mwu = dynamic_cast<MWUFramework*>(solver.get())) {
            mwu->printTimeStats();
        }

        // --- optional: demand model evaluation ---
        auto result = HandleDemandModel(argc, argv, cfg, *G, scheme);
        printStatsForDemandModel(argv, result);
    }
    return 0;
}
