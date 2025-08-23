//
// Created by Mert Biyikli on 18.08.25.
//

#include "raecke_solver.h"

void RaeckeSolver::solve(const Graph& graph) {
    // Initialize the RaeckeFRT object
    raeckeFRT.setGraph(graph);

    // Run the FRT algorithm
    raeckeFRT.run();
}

void RaeckeSolver::storeFlow() {
     double sumOfLambdas = 0.0;
    for(size_t i = 0; i < raeckeFRT.getTrees().size(); ++i) {
        double lambda_i = raeckeFRT.getLambdas()[i];
        sumOfLambdas += lambda_i;
        double normalized_lambda = lambda_i / sumOfLambdas; // Normalize by the last lambda (which should be 1.0)
        auto t = raeckeFRT.getTrees()[i];
        auto copyGraph = raeckeFRT.getGraphs()[i];
        raeckeTransform.addTree(raeckeFRT.getTrees()[i], normalized_lambda, raeckeFRT.getGraphs()[i]);

        if (debug) {
            // print out everything
            std::cout << "Tree " << i << ": Lambda = " << normalized_lambda << ", Graph = " << copyGraph.getNumNodes() << " nodes, "
                      << copyGraph.getNumEdges() << " edges.\n";
            std::cout << "Tree " << i << ": Max Rload = " << raeckeFRT.getMaxRload(i, t) << std::endl;
        }

    }


    // compute the congestion for the Raecke solution
    // given the demand map
    auto const& routingRaecke = raeckeTransform.getRoutingTable();

    for (const auto& [edge, demandMap] : routingRaecke) {
        for (const auto& [d, fraction] : demandMap) {

            if (debug)
                std::cout << "Edge (" << edge.first << ", " << edge.second << ") with demand (" << d.first << ", " << d.second << ") has flow fraction: " << fraction << "\n";

            f_e_st[edge][d]=fraction;

        }
    }
}
