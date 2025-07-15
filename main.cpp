#include <iostream>
#include "src/graph.h"
#include "test/randomgraphs/random_graph_generator.h"
#include "src/lp_solver/LPSolver.h"
#include "src/lp_solver/MCCF_lp_solver.h"
#include "src/tree_based/mcct/mcct_derandomized_weighted_solver.h"
#include "src/tree_based/raecke_tree_decomp.h"
#include "src/tree_based/raecke_transform.h"
#include "src/electrical/electrical_flow_naive.h"

// Verifies for each directed demand s→t that net‐out(s)=+1, net‐in(t)=+1, intermediates 0
void verifyFlowConservation(
        const EdgeDemandMap& routing,
        int N,
        double tol = 1e-9)
{
    // collect all directed demands
    std::set<std::pair<int,int>> demands;
    for (auto const& [e, netMap] : routing)
        for (auto const& [tp, f] : netMap)
            if (std::fabs(f) > tol)
                demands.insert({tp.first, tp.second});

    for (auto [s,t] : demands) {
        std::vector<double> net(N, 0.0);
        // for each arc (u→v), decompose signed f into u→v or v→u
        for (auto const& [e, netMap] : routing) {
            auto it = netMap.find({std::min(s,t),std::max(s,t)});
            if (it==netMap.end()) continue;
            double f = it->second;
            // if (u<v) f is for u→v; if u=>v, f is negative for (u,v) in the undirected graph,
            // but we recorded directed so:
            if (e.first < e.second) {
                // only positive on (u→v)
                net[e.first]  -= std::min(0.0, f);
                net[e.second] += std::min(0.0, f);
                net[e.first]  += std::max(0.0, f);
                net[e.second] -= std::max(0.0, f);
            } else {
                // invert sign
                net[e.first]  -= std::min(0.0, -f);
                net[e.second] += std::min(0.0, -f);
                net[e.first]  += std::max(0.0, -f);
                net[e.second] -= std::max(0.0, -f);
            }
        }

        bool ok = true;
        if (std::fabs(net[s] - 1.0) > tol) {
            std::cout << "[BAD] demand "<<s<<"→"<<t<<": net-out@"<<s<<"="<<net[s]<<" (≠1)\n";
            ok = false;
        }
        if (std::fabs(net[t] + 1.0) > tol) {
            std::cout << "[BAD] demand "<<s<<"→"<<t<<": net-out@"<<t<<"="<<net[t]<<" (≠–1)\n";
            ok = false;
        }
        for (int v = 0; v < N; ++v) {
            if (v==s||v==t) continue;
            if (std::fabs(net[v]) > tol) {
                std::cout << "[BAD] demand "<<s<<"→"<<t<<": node "<<v<<" net="<<net[v]<<" (≠0)\n";
                ok = false;
            }
        }
        if (ok) {
            std::cout << "[OK] demand "<<s<<"→"<<t<<" conserves flow.\n";
        }
    }
}



int main() {

    // Julia context
    jl_init();

    // Possibly preload modules
    jl_eval_string("using LapSolve");
    jl_eval_string("include(\"/Users/halilibrahim/Desktop/Thesis/ObliviousRouting/src/electrical/LaplacianSolver.jl\")");
    jl_eval_string("using .LaplacianSolver");
    // Graph G = RandomGraphGenerator::generate(10, 20, 1.0, 10.0);

    int n = 5;
    auto g = std::make_shared<Graph>(n);
    //G.readGraph("../test/randomgraphs/small/graph_10_6.dimacs");

    double c = 1;
    for(int i = 0; i < n; ++i) {
/*
        if(i == 2) {
            g->addEdge(i, (i+1)%n, c/2);
            continue;
        }*/

        g->addEdge(i, (i+1)%n, c);

    }


    auto raeckGraph = g->getRaeckeGraph();

    LaplacianSolver sol;
    std::unordered_map<std::pair<int, int>, double> edgeWeights;
    for(int i = 0; i < n; ++i) {
        for (const auto& u : raeckGraph.neighbors(i)) {
            if (i < u) { // Avoid double setting weights for undirected edges
                edgeWeights[{i, u}] = static_cast<double>(rand() % 10 + 1);
            }
        }
    }
    sol.init(raeckGraph, edgeWeights, true);
    Eigen::VectorXd x;
    // init x randomly
    x = Eigen::VectorXd::Random(raeckGraph.getNumNodes());
    auto result = sol.solve(x, true);

    std::cout << "Result of Laplacian solver:\n" << result << std::endl;

    sol.updateEdgeWeights({0, 1}, 5.0);
    std::cout << "Result after updating edge :\n " << sol.solve(x, true) << std::endl;


    jl_atexit_hook(0);
/*
    ElectricalFlowNaive efSolver;
    efSolver.init(raeckGraph);
    efSolver.run();
*/
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

/*
    std::cout << "Running CMMF Solver..." << std::endl;
    CMMF_Solver lp_solver;
    lp_solver.AddDemands(Demand{0, 4},c);
    lp_solver.AddDemands(Demand{2, 6},c);
    lp_solver.solve(*g);

    std::cout << "Running Cohen LP Solver..." << std::endl;
    LPSolver LP;
    LP.solve(*g);

*/

    return 0;
}
