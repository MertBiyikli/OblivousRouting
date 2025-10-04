#include <iostream>
#include "src/graph.h"
#include "src/lp_solver/LPSolver.h"
#include "src/lp_solver/MCCF_lp_solver.h"
#include "src/tree_based/frt/raecke_frt_solver.h"
#include "src/tree_based/random_mst/raecke_mst_solver.h"
#include "src/electrical/electrical_flow_naive.h"
#include "src/electrical/electrical_flow_optimized.h"
#include "experiments/performance/LP_Oblivious_Ratio.h"
#include "experiments/performance/linear_oblivious_routing_ratio.h"
#include "src/parse_parameter.h"




int main(int argc, char **argv) {

    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = std::chrono::high_resolution_clock::now();

    std::string err;
    auto cfg = parse_parameter(argc, argv, &err);
    if (!cfg) { std::cerr << err; return -1; } // EX_USAGE

    Graph g;
    g.readLFGFile(cfg->filename, /*undirected?*/ true);
    std::cout << "Graph loaded: " << g.getNumNodes() << " nodes, " << g.getNumEdges() << " edges.\n";
    // g.print();

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
    //solver->storeFlow();
    //solver->printFlow();
/*
    // compute worst case demand for linear oblivious routing
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
    //HandleDemandModel(argc, argv, cfg, g, solver);

    std::cout << "MWU number of iterations: " << solver->GetIterationCount() << std::endl;

    // compute the average oracle running time
    double mean = 0;
    for (const auto& t : solver->oracle_running_times) {
        mean += t;
    }
    mean /= solver->oracle_running_times.size();
    std::cout << "Average oracle time: " << mean << " [milliseconds]" << std::endl;


    return 0;
}
