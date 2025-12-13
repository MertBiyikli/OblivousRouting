#include <iostream>
#include "src/datastructures/GraphADJ.h"
#include "src/electrical/electrical_flow_optimized.h"
#include "src/datastructures/GraphCSR.h"
#include "src/parse_parameter.h"

/*
#include "src/lp_solver/MCCF_lp_solver.h"
#include "src/lp_solver/LPSolver.h"
#include "src/tree_based/frt/raecke_frt_solver.h"
#include "src/tree_based/random_mst/raecke_mst_solver.h"
#include "src/electrical/electrical_flow_naive.h"
#include "src/electrical/electrical_flow_parallel.h"
#include "experiments/performance/LP_Oblivious_Ratio.h"
#include "experiments/performance/linear_oblivious_routing_ratio.h"
#include "src/electrical/electrical_flow_parallel_on_the_fly.h"
#include "src/tree_based/ckr/ckr_partition.h"
#include "src/tree_based/ckr/ckr_tree_decomposer.h"
#include "src/utils/lca_datstructure.h"
#include "src/tree_based/ckr/raecke_ckr_solver.h"
// #include "src/tree_based/ckr/raecke_ckr_optimized.h"

// #include "src/tree_based/optimized_versions/frt/raecke_frt_opt.h"
#include "src/tree_based/optimized_versions/frt/mcct_optimized.h"
// #include "src/tree_based/optimized_versions/raecke_optimized_base.h"
#include "src/tree_based/ckr/utils/ultrametric_tree.h"
#include "src/tree_based/ckr/optimized/efficient_raecke_ckr.h"
*/
#include "src/lp_solver/LPSolver.h"
#include "src/lp_solver/MCCF_lp_solver.h"
#include "src/tree_based/ckr/raecke_mwu_ckr.h"
#include "src/tree_based/ckr/raecke_oracle_ckr.h"
#include "src/tree_based/frt/raecke_mwu_frt.h"
#include "src/tree_based/frt/raecke_oracle_frt.h"
#include "src/tree_based/random_mst/raecke_mwu_random.h"

typedef std::chrono::high_resolution_clock::time_point TimeVar;

#define duration(a) std::chrono::duration_cast<std::chrono::milliseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

template<typename F, typename... Args>
double funcTime(F func, Args&&... args){
    TimeVar t1= timeNow();
    func(std::forward<Args>(args)...);
    return duration(timeNow()-t1);
}

int main(int argc, char **argv) {

    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = std::chrono::high_resolution_clock::now();

    std::string err;
    auto cfg = parse_parameter(argc, argv, &err);
    if (!cfg) { std::cerr << err; return -1; } // EX_USAGE

    std::unique_ptr<IGraph> g = std::make_unique<GraphADJ>();
    readLGFFile(*g, cfg->filename, /*undirected?*/ true);
    std::cout << "Graph loaded: " << g->getNumNodes() << " nodes, " << g->getNumEdges() << " edges.\n";
    //g.print();

    GraphCSR g_csr;// = std::make_unique<GraphCSR>();
    readLGFFile(g_csr, cfg->filename,  true);
    g_csr.finalize();

    std::unique_ptr<ObliviousRoutingSolver> solver;
    switch (cfg->solver) {
        /*
        case SolverType::ELECTRICAL_NAIVE:
            std::cout << "Running Electrical Flow (naive)...\n";
            solver = std::make_unique<ElectricalFlowNaive>();
            break;
        case SolverType::RAECKE_FRT:
            std::cout << "Running Tree-based (Raecke/FRT)...\n";
            solver = std::make_unique<RaeckeFRTSolver>(g_csr);
            break;
        case SolverType::LP_APPLEGATE_COHEN:
            std::cout << "Running LP (Applegate–Cohen)...\n";
            solver = std::make_unique<LPSolver>(g_csr);
            break;
        case SolverType::RAECKE_RANDOM_MST:
            std::cout << "Running Tree-based (Raecke/MST)...\n";
            solver = std::make_unique<RaeckeMSTSolver>();
            break;
            */

        case SolverType::ELECTRICAL_OPTIMIZED:
            std::cout << "Running Electrical Flow (optimized)...\n";
            solver = std::make_unique<ElectricalFlowOptimized>(g_csr, 0, false);
            break;
            /*
        case SolverType::ELECTRICAL_PARALLEL_BATCHES:
            std::cout << "Running Electrical Flow (parallel - batches)...\n";
            solver = std::make_unique<ElectricalFlowParallel>(g_csr);
            break;
        case SolverType::ELECTRICAL_PARALLEL_ONTHEFLY:
            std::cout << "Running Electrical Flow (parallel - on the fly)...\n";
            solver = std::make_unique<ElectricalFlowParallelOnTheFly>(g_csr);
            break;
        case SolverType::RAECKE_CKR:
            std::cout << "Running Tree-based (Raecke/CKR)... \n";
            solver = std::make_unique<RaeckeCKRSolver>(g_csr);
            break;

        case SolverType::RAECKE_CKR_OPTIMIZED:
            std::cout << "Running Tree-based (Raecke/CKR Optimized)... \n";
            solver = std::make_unique<EfficientRaeckeCKR>(g_csr);
            break;
            */

/*
        case SolverType::OPTIMIZED_RAECKE_FRT:
            std::cout << "Running Optimized Raecke FRT... \n";
            solver = std::make_unique<RaeckeFRTOptimized>();
            break;*/
        default:
            break;
    }

    if ( !solver ) {
        std::cerr << "Solver not implemented.\n";
        return -1; // EX_UNAVAILABLE
    }

    // run the solver
    start_time = std::chrono::high_resolution_clock::now();
    // solver->debug = true;
    //solver->(*g_csr);
    auto electrical_scheme = solver->solve();
    end_time = std::chrono::high_resolution_clock::now();
    std::cout << "Electrical Running time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " [milliseconds]" << std::endl;

    DemandMap demand_map;
    for (int s = 0; s < g_csr.getNumNodes(); ++s) {
        for (int t = 0; t < g_csr.getNumNodes(); ++t) {
            if (s >= t) continue;
            demand_map.addDemand(s, t, 1.0);
        }
    }

    std::vector<double> cong_electrical;
    electrical_scheme->routeDemands(cong_electrical, demand_map);


    double max_congestion_electrical = electrical_scheme->getMaxCongestion(cong_electrical);
    std::cout << "Max congestion (electrical): " << max_congestion_electrical << std::endl;


    // std::unique_ptr<RaeckeOracle> ckrOracle = std::make_unique<CKROracle>(g_csr);
    auto ckrsolver = std::make_unique<RaeckeMWU_CKR>(g_csr, 0);//, std::move(ckrOracle));
    start_time = std::chrono::high_resolution_clock::now();
    auto ckr_scheme = ckrsolver->solve();
    end_time = std::chrono::high_resolution_clock::now();
    std::cout << "CKR Running time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " [milliseconds]" << std::endl;


    std::vector<double> congestion_ckr;
    ckr_scheme->routeDemands(congestion_ckr, demand_map);

    double max_congestion_ckr = ckr_scheme->getMaxCongestion(congestion_ckr);
    std::cout << "Max congestion (ckr): " << max_congestion_ckr << std::endl;

    auto frtsolver = std::make_unique<RaeckeMWU_FRT>(g_csr, 0);//, std::move(ckrOracle));
    start_time = std::chrono::high_resolution_clock::now();
    auto frt_scheme = frtsolver->solve();
    end_time = std::chrono::high_resolution_clock::now();
    std::cout << "FRT Running time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " [milliseconds]" << std::endl;

    std::vector<double> frt_congestion;
    frt_scheme->routeDemands(frt_congestion, demand_map);


    double max_congestion_frt = frt_scheme->getMaxCongestion(frt_congestion);
    std::cout << "Max congestion (frt): " << max_congestion_frt << std::endl;



    auto randomTreeSolver = std::make_unique<RaeckeMWU_Random>(g_csr, 0);//, std::move(ckrOracle));
    start_time = std::chrono::high_resolution_clock::now();
    auto random_scheme = randomTreeSolver->solve();
    end_time = std::chrono::high_resolution_clock::now();
    std::cout << "Random MST Running time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " [milliseconds]" << std::endl;

    std::vector<double> random_congestion;
    random_scheme->routeDemands(random_congestion, demand_map);
    for (int e = 0; e < g_csr.getNumEdges(); ++e) {
        //std::cout << "Edge " << e << " outflow: " << frt_congestion[e] << std::endl;
    }

    double max_congestion_random = random_scheme->getMaxCongestion(random_congestion);
    std::cout << "Max congestion (random): " << max_congestion_random << std::endl;


    // LP Applegate and Cohen
    auto lp_solver = std::make_unique<LPSolver>(g_csr);
    start_time = std::chrono::high_resolution_clock::now();
    auto lp_scheme = lp_solver->solve();
    end_time = std::chrono::high_resolution_clock::now();
    std::cout << "LP  Running time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " [milliseconds]" << std::endl;

    std::vector<double> lp_congestion;
    lp_scheme->routeDemands(lp_congestion, demand_map);
    double max_congestion_lp = lp_scheme->getMaxCongestion(lp_congestion);
    std::cout << "Max congestion (lp-applegate & cohen): " << max_congestion_lp << std::endl;


    // LP offline optimal solution
    auto offline_lp_solver = std::make_unique<CMMF_Solver>(g_csr);
    for (int i = 0; i < demand_map.size(); ++i) {
        offline_lp_solver->AddDemands(demand_map.getDemandPair(i), demand_map.getDemandValue(i));
    }
    start_time = std::chrono::high_resolution_clock::now();
    auto offline_scheme = offline_lp_solver->solve();
    end_time = std::chrono::high_resolution_clock::now();

    std::cout << "Offline LP Running time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " [milliseconds]" << std::endl;

    std::vector<double> offline_congestion_vector;
    offline_scheme->routeDemands(offline_congestion_vector, demand_map);
    double offline_congestion = offline_scheme->getMaxCongestion(offline_congestion_vector);
    std::cout << "Max congestion (lp-offline): " << offline_congestion << std::endl;


    // TODO: this store flow is a major bottleneck for the tree based experiments.
/*
    solver->storeFlow();
    solver->printFlow();
*/

    // verify flow conservation


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
    // HandleDemandModel(argc, argv, cfg, g_csr, solver);


/*
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
*/

    /* Developing the efficient ckr raecke shit
     *
     */
/*
    std::cout << "Now running the CKR partition..." << std::endl;
    EfficientRaeckeCKR efficient_raecke_ckr_solver;
    efficient_raecke_ckr_solver.debug = false;
    start_time = std::chrono::high_resolution_clock::now();
    efficient_raecke_ckr_solver.setGraph(g_csr);
    efficient_raecke_ckr_solver.solve();
    end_time = std::chrono::high_resolution_clock::now();

    efficient_raecke_ckr_solver.storeFlow();
    efficient_raecke_ckr_solver.printFlow_();


    std::cout << "CKR Solver number of iterations: " << efficient_raecke_ckr_solver.GetIterationCount() << std::endl;
    // print out statistics
    double total_time_msec = std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();
    std::cout << "Total running time: " << total_time_msec << " [milliseconds]" << std::endl;

    double mean_oracle_time = 0;
    for (const auto& t : efficient_raecke_ckr_solver.oracle_running_times) {
        mean_oracle_time += t;
    }
    mean_oracle_time /= efficient_raecke_ckr_solver.oracle_running_times.size();
    std::cout << "Average oracle time: " << mean_oracle_time << " [milliseconds]" << std::endl;

    double mean_pure_oracle_time = 0;
    for (const auto& t : efficient_raecke_ckr_solver.pure_oracle_running_times) {
        mean_pure_oracle_time += t;
    }
    mean_pure_oracle_time /= efficient_raecke_ckr_solver.pure_oracle_running_times.size();
    std::cout << "Average pure oracle time: " << mean_pure_oracle_time << " [milliseconds]" << std::endl;
    return 0;*/
}
