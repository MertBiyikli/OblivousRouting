#include <iostream>
#include "src/graph.h"
#include "src/lp_solver/LPSolver.h"
#include "src/lp_solver/MCCF_lp_solver.h"
#include "src/tree_based/frt/raecke_frt_solver.h"
#include "src/tree_based/random_mst/raecke_mst_solver.h"
#include "src/electrical/electrical_flow_naive.h"
#include "src/electrical/electrical_flow_optimized.h"
#include "src/electrical/electrical_flow_parallel.h"
#include "experiments/performance/LP_Oblivious_Ratio.h"
#include "experiments/performance/linear_oblivious_routing_ratio.h"
#include "src/parse_parameter.h"
#include "src/electrical/electrical_flow_parallel_on_the_fly.h"
#include "src/tree_based/ckr/ckr_partition.h"
#include "src/tree_based/ckr/ckr_tree_decomposer.h"
#include "src/utils/lca_datstructure.h"
#include "src/tree_based/ckr/raecke_ckr_solver.h"
#include "src/tree_based/ckr/raecke_ckr_optimized.h"


int main(int argc, char **argv) {

    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = std::chrono::high_resolution_clock::now();

    std::string err;
    auto cfg = parse_parameter(argc, argv, &err);
    if (!cfg) { std::cerr << err; return -1; } // EX_USAGE

    Graph g;
    g.readLFGFile(cfg->filename, /*undirected?*/ true);
    std::cout << "Graph loaded: " << g.getNumNodes() << " nodes, " << g.getNumEdges() << " edges.\n";
    //g.print();

    std::unique_ptr<ObliviousRoutingSolver> solver;
    switch (cfg->solver) {
        case SolverType::ELECTRICAL_NAIVE:
            std::cout << "Running Electrical Flow (naive)...\n";
            solver = std::make_unique<ElectricalFlowNaive>();
            break;
        case SolverType::RAECKE_FRT:
            std::cout << "Running Tree-based (Raecke/FRT)...\n";
            solver = std::make_unique<RaeckeFRTSolver>();
            break;
        case SolverType::LP_APPLEGATE_COHEN:
            std::cout << "Running LP (Applegate–Cohen)...\n";
            solver = std::make_unique<LPSolver>();
            break;
        case SolverType::RAECKE_RANDOM_MST:
            std::cout << "Running Tree-based (Raecke/MST)...\n";
            solver = std::make_unique<RaeckeMSTSolver>();
            break;
        case SolverType::ELECTRICAL_OPTIMIZED:
            std::cout << "Running Electrical Flow (optimized)...\n";
            solver = std::make_unique<ElectricalFlowOptimized>();
            break;
        case SolverType::ELECTRICAL_PARALLEL_BATCHES:
            std::cout << "Running Electrical Flow (parallel - batches)...\n";
            solver = std::make_unique<ElectricalFlowParallel>();
            break;
        case SolverType::ELECTRICAL_PARALLEL_ONTHEFLY:
            std::cout << "Running Electrical Flow (parallel - on the fly)...\n";
            solver = std::make_unique<ElectricalFlowParallelOnTheFly>();
            break;
        case SolverType::RAECKE_CKR:
            std::cout << "Running Tree-based (Raecke/CKR)... \n";
            solver = std::make_unique<RaeckeCKROptimized>();
        default:
            break;
    }

    if ( !solver ) {
        std::cerr << "Solver not implemented.\n";
        return -1; // EX_UNAVAILABLE
    }

    // run the solver
    start_time = std::chrono::high_resolution_clock::now();
    //solver->debug = true;
    solver->solve(g);
    end_time = std::chrono::high_resolution_clock::now();
    std::cout << "Running time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " [milliseconds]" << std::endl;

    // TODO: this store flow is a major bottleneck for the tree based experiments.
    // solver->storeFlow();
    // solver->printFlow_();

    // verify flow conservation

    /*
    std::cout << "Testing the CKR Partitioning module...\n";
    // init distances
    double d = 1;
    for (int s = 0; s<g.getNumNodes(); s++) {
        for (auto& v : g.neighbors(s)) {
            if (s > v) continue;
            g.updateEdgeDistance(s, v, d);
            d+=1.0;
        }
    }

    g.createDistanceMatrix();
    g.print();
    // init random subset of vertices
    std::vector<int> X;
    for (int i = 0; i<g.getNumNodes()/2; i++) {
        X.push_back(i);
    }
    CKRPartion ckr;
    ckr.init(g, false);
    double delta = g.GetDiameter();
    // print out the input set
    std::cout << "Input set X:\n";
    for (const auto& x : X) {
        std::cout << x << " ";
    }
    std::cout << "\nDelta: " << delta << "\n";
    auto cluster = ckr.computePartition(X, delta);

    // print out the cluster
    std::cout << "Cluster assignments:\n";
    for (size_t i = 0; i<cluster.size(); i++) {
        std::cout << "Node " << X[i] << " -> Cluster " << cluster[i] << "\n";
    }


    int root(0);

    std::vector<lca::Edge> edges;
    RandomMST mst_tree;
    auto mst_edges = mst_tree.build_mst(g);
    for (auto [u, v]:mst_edges) {
        edges.push_back({ u, v});
    }
    lca::LCA_n_1 tree(edges, root);

    int q = 3;
    for (int i = 0; i < q; ++i) {
        std::cout << "Query " << i << "th, (u, v) = "<< '\n';
        int u = rand()%g.getNumNodes(), v=rand()%g.getNumNodes();
        std::cout << "LCA(" << u << ", " << v << ") = " << tree.LCA(u, v) << '\n';
    }

    //delta /= 2;
    TreeDecomposer decomposer;
    std::vector<int> node_ids(g.getNumNodes());
    std::iota(node_ids.begin(), node_ids.end(), 0);
    TreeNode* DecompTree = decomposer.decompose(g, delta, node_ids);


    std::cout << "CKR Decomposition Tree:" << std::endl;
    print_tree(DecompTree);
    // compute worst case demand for linear oblivious routing

    RaeckeCKRSolver sol;
    sol.solve(g);
    sol.storeFlow();
    sol.printFlow_();*/
    /*
    LinearObliviousRatio ratio;
    ratio.init(g, solver->f_e_st);
    std::cout << "Linear oblivious routing worst case demand congestion: " << ratio.solve() << std::endl;
*/
/*
    // compute worst case demand set
    ObliviousRatio OR;
    OR.init(g, solver->f_e_st);
    std::cout << "Worst case demand congestion: " << OR.solve() << std::endl;
*/
    // if a demand model is provided, compute the oblivious ratio for that demand model
    // HandleDemandModel(argc, argv, cfg, g, solver);


    if (cfg->solver != SolverType::LP_APPLEGATE_COHEN) {
        std::cout << "MWU number of iterations: " << solver->GetIterationCount() << std::endl;

        // compute the average oracle running time
        double mean = 0;
        for (const auto& t : solver->oracle_running_times) {
            mean += t;
        }
        mean /= solver->oracle_running_times.size();
        std::cout << "Average oracle time: " << mean << " [milliseconds]" << std::endl;

        // compute the average pure oracle running time
        double pure_mean = 0;
        for (const auto& t : solver->pure_oracle_running_times) {
            pure_mean += t;
        }
        pure_mean /= solver->pure_oracle_running_times.size();
        std::cout << "Average pure oracle time: " << pure_mean << " [milliseconds]" << std::endl;
    }

    return 0;
}
