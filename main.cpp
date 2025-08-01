#include <iostream>
#include "src/graph.h"
#include "test/randomgraphs/random_graph_generator.h"
#include "src/lp_solver/LPSolver.h"
#include "src/lp_solver/MCCF_lp_solver.h"
#include "src/tree_based/mcct/mcct_derandomized_weighted_solver.h"
#include "src/tree_based/raecke_tree_decomp.h"
#include "src/tree_based/raecke_transform.h"
#include "src/electrical/electrical_flow_naive.h"



int main() {

    // Julia context
    jl_init();

    // Possibly preload modules
    jl_eval_string("using LapSolve");
    jl_eval_string("include(\"/Users/halilibrahim/Desktop/Thesis/ObliviousRouting/src/electrical/LaplacianSolver.jl\")");
    jl_eval_string("using .LaplacianSolver");
    // Graph G = RandomGraphGenerator::generate(10, 20, 1.0, 10.0);

    auto g = std::make_shared<RaeckeGraph>();
    g->readLFGFile("../experiments/random/Deltacom.lgf", true);
    int n = g->getNumNodes();
/*
    double c = 1;
    for(int i = 0; i < n; ++i) {
        g->addEdge(i, (i+1)%n, c);

    }*/


    //auto raeckGraph = g->getRaeckeGraph();

    g->print();

    std::unordered_map<std::pair<int, int>, double> edgeWeights;
    for(int i = 0; i < n; ++i) {
        for (const auto& u : g->neighbors(i)) {
            if (i < u) { // Avoid double setting weights for undirected edges
                edgeWeights[{i, u}] = static_cast<double>(rand() % 10 + 1);
            }
        }
    }

    /*
    sol.init(raeckGraph, edgeWeights, true);
    Eigen::VectorXd x;

    x = Eigen::VectorXd::Random(raeckGraph.getNumNodes());
    auto result = sol.solve(x, true);

    std::cout << "Result of Laplacian solver:\n" << result << std::endl;

    sol.updateEdgeWeights({0, 1}, 5.0);
    std::cout << "Result after updating edge :\n " << sol.solve(x, true) << std::endl;

*/

    auto start_time = std::chrono::high_resolution_clock::now();
    ElectricalFlowNaive efSolver;
    efSolver.init(*g, false);
    efSolver.run();
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Electrical Flow Naive solver took " << elapsed.count() << " seconds.\n";
    /*
    // get the load for commodity (0, 1), (0, 2)
    std::vector<std::pair<int, int>> commodity;
    for(int i = 0; i < n; ++i) {
        for(int j = i+1; j < n; ++j) {
            commodity.emplace_back(i, j);
        }
    }

    auto routing = efSolver.getRoutingForCommodity(commodity);

    for (const auto& [edge, flow] : routing) {
        std::cout << "commodity (" << edge.first << ", " << edge.second << "): Flow = \n";
        for(int edge_id = 0; edge_id < flow.size(); edge_id++) {
            std::cout << " edge (" << efSolver.edges[edge_id].first << " / " << efSolver.edges[edge_id].second << " ) : " << flow.transpose()(edge_id) << std::endl;
        }
        std::cout << std::endl;
    }
    */


    std::cout << "Max congestion electrical flow: " << efSolver.getMaximumCongestion() << std::endl;


    jl_atexit_hook(0);

    /*
    // set random edge distances for debugging purposes
    for (int& v : raeckGraph.getVertices()) {
        for (const auto& u : raeckGraph.neighbors(v)) {
            if (v < u) { // Avoid double setting distances for undirected edges
                raeckGraph.updateEdgeDistance(v, u,  static_cast<double>(rand() % 10 + 1));
            }
        }
    }
    raeckGraph.print();
    */
/*
    RaeckeFRT solver;
    solver.setGraph(raeckGraph);
    solver.run();

    for(double& beta : solver.getLambdas()) {
        std::cout << "Lambda: " << beta << std::endl;
    }

    RaeckeSolutionTransform transform;
    double sumOfLambdas = 0.0;
    for(size_t i = 0; i < solver.getTrees().size(); ++i) {
        double lambda_i = solver.getLambdas()[i];
        sumOfLambdas += lambda_i;
        double normalized_lambda = lambda_i / sumOfLambdas; // Normalize by the last lambda (which should be 1.0)
        auto t = solver.getTrees()[i];
        auto copyGraph = solver.getGraphs()[i];
        transform.addTree(solver.getTrees()[i], normalized_lambda, solver.getGraphs()[i]);
    }

    auto const& routing = transform.getRoutingTable();
    std::cout << "Räcke's oblivious routing table:" << std::endl;
    for (const auto& [edge, demandMap] : routing) {
        std::cout << "Edge (" << edge.first << ", " << edge.second << "): ";
        for (const auto& [demand, fraction] : demandMap) {
            std::cout << "[" << demand.first << " → " << demand.second << "] = "
                      << fraction << ", ";
        }
        std::cout << std::endl;
    }

    verifyFlowConservation(routing, n);
*/


    // this will enumerate and print every commodity’s flow‐paths:

    /*
    g->addEdge(0,1, 100);
    g->addEdge(0,2, 100);
    g->addEdge(1,2, 100);
    g->addEdge(0,3, 100);
     */
/*
    for(int i = 0; i < n; ++i) {
        for (auto& edge : g->neighbors(i)) {
            edge->distance = 1;
        }
    }*/
    /*
    // Step 3: Set edge (0, 2) to have higher distance (20)
    int u = 0, v = 2;
    auto edge = std::find_if(
            g->neighbors(u).begin(), g->neighbors(u).end(),
            [v](const std::shared_ptr<Edge>& e) { return e->target == v; });
    if (edge != g->neighbors(u).end()) {
        (*edge)->distance = 20;
    }

    auto rev_edge = std::find_if(
            g->neighbors(v).begin(), g->neighbors(v).end(),
            [u](const std::shared_ptr<Edge>& e) { return e->target == u; });
    if (rev_edge != g->neighbors(v).end()) {
        (*rev_edge)->distance = 20;
    }*/
    /*
    MCCTDerandomizedWeightedSolver solver;

    solver.setGraph(g);
    solver.addDemand(0, 1, 1.0);
    auto t = solver.getBestTree();
    t.ExportToDot("DecompositionTree.dot");
*/


    /*
    RackeObliviousRoutingSolver rackeSolver;
    rackeSolver.setGraph(g);
    rackeSolver.createTreesAndLambda();

    rackeSolver.printObliviousRoutingTable();
    std::cout << "Räcke's oblivious routing table:" << std::endl;
    RaeckeTransform transform;
    DemandMap routing;

    double lambda_sum = 0.0;
    for (size_t i = 0; i < rackeSolver.getTrees().size(); ++i) {
        double lambda_i = rackeSolver.getLambdas()[i];
        lambda_sum += lambda_i;
        double normalized_lambda = lambda_i / lambda_sum;

        routing = transform.addTree(rackeSolver.getTrees()[i], normalized_lambda, rackeSolver.getGraphs()[i]);
    }

    // Step 4: Done. `routing` contains the oblivious routing table.
    printRouting(routing);
     */


    // uncomment for full view
/*
    std::cout << "Running CMF Solver..." << std::endl;
    int c = 1.0;
    CMMF_Solver lp_solver;
    lp_solver.init(*g);
    lp_solver.setDebug(false);
    for(int i = 0; i<g->getNumNodes(); i++) {
        for(int j = i+1; j<g->getNumNodes(); j++) {
            lp_solver.AddDemands(Demand{i, j},c);
        }
    }

    start_time = std::chrono::high_resolution_clock::now();
    lp_solver.solve(*g);
    end_time = std::chrono::high_resolution_clock::now();
    std::cout << "CMF LP solver took " << std::chrono::duration<double>(end_time-start_time).count() << " seconds. " << std::endl;


    std::cout << "Running Cohen LP Solver..." << std::endl;
    LPSolver LP;
    LP.setDebug(false);
    start_time = std::chrono::high_resolution_clock::now();
    LP.solve(*g);
    end_time = std::chrono::high_resolution_clock::now();
    std::cout << "Cohen LP solver took " << std::chrono::duration<double>(end_time-start_time).count() << " seconds. " << std::endl;
    std::cout << "Maximum congestion using Applegate & Cohen: " << LP.getMaximumCongestion(*g) << std::endl;


    std::cout << "Running Raecke solver..." << std::endl;
    start_time = std::chrono::high_resolution_clock::now();
    RaeckeFRT raecke;
    raecke.setGraph(*g);
    raecke.run();

    end_time = std::chrono::high_resolution_clock::now();
    std::cout << "Raecke (FRT) solver took " << std::chrono::duration<double>(end_time-start_time).count() << " seconds. " << std::endl;

    RaeckeSolutionTransform transform;
    double sumOfLambdas = 0.0;
    for(size_t i = 0; i < raecke.getTrees().size(); ++i) {
        double lambda_i = raecke.getLambdas()[i];
        sumOfLambdas += lambda_i;
        double normalized_lambda = lambda_i / sumOfLambdas; // Normalize by the last lambda (which should be 1.0)
        auto t = raecke.getTrees()[i];
        auto copyGraph = raecke.getGraphs()[i];
        transform.addTree(raecke.getTrees()[i], normalized_lambda, raecke.getGraphs()[i]);
    }

    auto const& routing = transform.getRoutingTable();
    std::unordered_map<std::pair<int, int>, double> totalFlow;
    for (const auto& [edge, demandMap] : routing) {
        for (const auto& [demand, fraction] : demandMap) {
            int u = edge.first, v = edge.second;

            if(u > v) {

                totalFlow[{v, u}] += fraction;
            }else{

                totalFlow[edge] += fraction;
            }

        }
    }
    double max_congestion_raecke = 0;

    for(auto& [edge, value] : totalFlow) {
        double cong = value/g->getEdgeCapacity(edge.first, edge.second);
        if(max_congestion_raecke < cong) {
            max_congestion_raecke = cong;
        }
    }


    std::cout << "Raecke generated: " << max_congestion_raecke << std::endl;
*/
    return 0;
}
