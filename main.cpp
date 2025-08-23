#include <iostream>
#include "src/graph.h"
#include "test/randomgraphs/random_graph_generator.h"
#include "src/lp_solver/LPSolver.h"
#include "src/lp_solver/MCCF_lp_solver.h"
#include "src/tree_based/raecke_solver.h"
#include "src/electrical/electrical_flow_naive.h"
#include <filesystem>
#include "experiments/performance/demands/GravityModel.h"
#include "experiments/performance/LP_Oblivious_Ratio.h"
#include "src/parse_parameter.h"




int main(int argc, char **argv) {

    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = std::chrono::high_resolution_clock::now();

    std::string err;
    auto cfg = parse_parameter(argc, argv, &err);
    if (!cfg) { std::cerr << err; return -1; } // EX_USAGE

    Graph g;
    g.readLFGFile(cfg->filename, /*undirected?*/ true);


    std::unique_ptr<ObliviousRoutingSolver> solver;
    switch (cfg->solver) {
        case SolverType::ELECTRICAL_NAIVE:
            std::cout << "Running Electrical Flow (naive)…\n";
            solver = std::make_unique<ElectricalFlowNaive>();
            break;
        case SolverType::RAECKE_FRT:
            std::cout << "Running Tree-based (Raecke/FRT)…\n";
            solver = std::make_unique<RaeckeSolver>();
            break;
        case SolverType::LP_APPLEGATE_COHEN:
            std::cout << "Running LP (Applegate–Cohen)…\n";
            solver = std::make_unique<LPSolver>();
            break;
    }

    // run the solver
    start_time = std::chrono::high_resolution_clock::now();
    solver->solve(g);
    end_time = std::chrono::high_resolution_clock::now();
    std::cout << "Running time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " [milliseconds]" << std::endl;
    solver->storeFlow();

    ObliviousRatio OR;
    OR.init(g, solver->f_e_st);
    std::cout << "Worst case demand congestion: " << OR.solve() << std::endl;

    return 0;
}
