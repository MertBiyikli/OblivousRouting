//
// Created by Mert Biyikli on 18.09.25.
//

#include "raecke_mst_solver.h"

/*
void RaeckeMSTSolver::solve(const Graph& graph) {


    raeckeMST.setGraph(graph);
    raeckeMST.run();

    this->iteration_count = raeckeMST.getTrees().size();

    this->oracle_running_times = raeckeMST.oracle_running_times;
}

void RaeckeMSTSolver::storeFlow() {
     double sumOfLambdas = 0.0;
    for(size_t i = 0; i < raeckeMST.getTrees().size(); ++i) {
        double lambda_i = raeckeMST.getLambdas()[i];
        sumOfLambdas += lambda_i;
        double normalized_lambda = lambda_i / sumOfLambdas; // Normalize by the last lambda (which should be 1.0)
        auto t = raeckeMST.getTrees()[i];
        auto copyGraph = raeckeMST.getGraphs()[i];
        raeckeTransform.addTree(raeckeMST.getTrees()[i], normalized_lambda, raeckeMST.getGraphs()[i]);

        if (debug) {
            // print out everything
            std::cout << "Tree " << i << ": Lambda = " << normalized_lambda << ", Graph = " << copyGraph.getNumNodes() << " nodes, "
                      << copyGraph.getNumEdges() << " edges.\n";
            std::cout << "Tree " << i << ": Max Rload = " << raeckeMST.getMaxRload(i) << std::endl;
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
*/