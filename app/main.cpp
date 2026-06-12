#include "../include/algorithms/mwu/tree_mwu.h"
#include "../include/io/parse_argurment_io.h"
#include "../include/core/routing_engine.h"
#include "../include/algorithms/parallel/mwu/par_electrical_flow.h"

#include "../include/algorithms/parallel/execution/openmp_static.h"
#include "../include/algorithms/parallel/execution/openmp_dynamic.h"
#include "../include/algorithms/parallel/execution/openmp_chunked.h"
#include "../include/algorithms/parallel/execution/openmp_guided.h"

int g_available_threads = 1;
int main(int argc, char **argv) {

    // Parse command line arguments
    Config cfg;
    auto graph = load_graph(cfg, argc, argv);

    // Precompute offline optimal congestion if needed
    offlineOptimal(graph, cfg);

#if _OPENMP
    std::cout << "Running with OpenMP support. Max threads: " << omp_get_max_threads() << "\n";
    g_available_threads = omp_get_max_threads();
#else
    std::cout << "Running without OpenMP support. Parallel execution will be disabled.\n
#endif

#if DEBUG
    cfg.debug = true;
#endif
    cfg.debug = true;

    double seq_run = 0;
    RoutingEngine engine;
    for (SolverType type : cfg.solvers ) {
        if (auto result = engine.solve(*graph, cfg, type)) {
            seq_run = result->total_runtime;
            for (const auto& [str, _type] : SOLVER_MAP) {
                if ( _type == type ) {
                    result->storeAsFile(str + "_result.json", OutPutFormat::JASON);
                    break;
                }
            }
        }
    }
    //return 0;

#if _OPENMP
    std::cout << "Running with OpenMP support. Max threads: " << omp_get_max_threads() << "\n";
#else
    std::cout << "Running without OpenMP support. Parallel execution will be disabled.\n
#endif

    auto t0 = timeNow();
    std::cout <<"OMP Electrical MWU:\n";
    ParElectricalFlowMWU omp_static(*graph, 0, true, OpenMPStaticExecution{g_available_threads});
    auto omp_scheme = omp_static.solve();
    std::cout << "Running time openmp (static): " << (duration(timeNow() - t0))/ 1e6 << " seconds\n";
    std::cout << "Speed up: " << (seq_run/ (duration(timeNow() - t0))) << "x\n";

    t0 = timeNow();
    ParElectricalFlowMWU omp_dynamic(*graph, 0, true, OpenMPDynamicExecution{g_available_threads});
    omp_scheme = omp_dynamic.solve();
    std::cout << "Running time openmp (dynamic): " << (duration(timeNow() - t0))/ 1e6 << " seconds\n";
    std::cout << "Speed up: " << (seq_run/ (duration(timeNow() - t0))) << "x\n";

    t0 = timeNow();
    ParElectricalFlowMWU omp_guided(*graph, 0, true, OpenMPGuidedExecution{g_available_threads});
    omp_scheme = omp_guided.solve();
    std::cout << "Running time openmp (guided): " << (duration(timeNow() - t0))/ 1e6 << " seconds\n";
    std::cout << "Speed up: " << (seq_run/ (duration(timeNow() - t0))) << "x\n";

    t0 = timeNow();
    ParElectricalFlowMWU omp_chunked(*graph, 0, true, OpenMPChunkedExecution{g_available_threads});
    omp_scheme = omp_chunked.solve();
    std::cout << "Running time openmp (chunked): " << (duration(timeNow() - t0))/ 1e6 << " seconds\n";
    std::cout << "Speed up: " << (seq_run/ (duration(timeNow() - t0))) << "x\n";
    //omp_scheme->printRoutingTable();
}
