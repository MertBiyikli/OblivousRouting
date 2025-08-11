#include <iostream>
#include "src/graph.h"
#include "test/randomgraphs/random_graph_generator.h"
#include "src/lp_solver/LPSolver.h"
#include "src/lp_solver/MCCF_lp_solver.h"
#include "src/tree_based/mcct/mcct_derandomized_weighted_solver.h"
#include "src/tree_based/raecke_tree_decomp.h"
#include "src/tree_based/raecke_transform.h"
#include "src/electrical/electrical_flow_naive.h"
#include <source_location>
#include <filesystem>

int main() {

    std::cout << std::filesystem::current_path() << "\n";

    auto g = std::make_shared<Graph>();
    g->readLFGFile("experiments/random/triangle.lgf", true);
    int n = g->getNumNodes();

    g->print();



    auto start_time = std::chrono::high_resolution_clock::now();
    ElectricalFlowNaive efSolver;
    efSolver.init(*g, false);
    efSolver.run();
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Electrical Flow Naive solver took " << elapsed.count() << " seconds.\n";

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
        double flow_sum = flow.cwiseAbs().sum();
        for(int edge_id = 0; edge_id < flow.size(); edge_id++) {
            std::cout << " edge (" << efSolver.edges[edge_id].first << " / " << efSolver.edges[edge_id].second << " ) : " << flow.transpose()(edge_id) << std::endl;
        }
        std::cout << "Total flow for commodity (" << edge.first << ", " << edge.second << ") = " << flow_sum << std::endl;
        std::cout << std::endl;
    }

    // print the average flow values
    std::cout << "Average flow values (first 5 edges):\n";
    for (auto &[e_u, flow]: efSolver.f_e_u) {
        if (flow != 0.0) { // Only print significant flows
            std::cout << "Edge (" << efSolver.edges[e_u.first].first << ", " << efSolver.edges[e_u.first].second << ") commodity ( "
                      << e_u.second << " -> " << efSolver.x_fixed << " ): Flow = " << flow << "\n";
        }
    }

    // also sum up the flow for each commodity and list the flow sum
    std::cout << "Total flow for each commodity:\n";
    std::unordered_map<std::pair<int, int>, double> flow_sum;
    for (auto &[e_u, flow]: efSolver.f_e_u) {
        flow_sum[{e_u.second, efSolver.x_fixed}] += std::abs(flow); // Sum up the flow for each commodity
    }

    for (const auto& [commodity, total_flow] : flow_sum) {
        std::cout << "Commodity (" << commodity.first << ", " << commodity.second << "): Total Flow = " << total_flow << "\n";
    }



    std::cout << "Max congestion electrical flow: " << efSolver.getMaximumCongestion() << std::endl;

    // uncomment for full view
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

    auto const& routingRaecke = transform.getRoutingTable();
    std::unordered_map<std::pair<int, int>, double> totalFlow;
    for (const auto& [edge, demandMap] : routingRaecke) {
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

    return 0;
}
