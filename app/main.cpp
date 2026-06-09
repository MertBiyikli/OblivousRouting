#include "../include/algorithms/mwu/tree_mwu.h"
#include "../include/io/parse_argurment_io.h"
#include "../include/core/routing_engine.h"
#include "../include/algorithms/parallel/mwu/par_electrical_flow.h"
int main(int argc, char **argv) {

    // Parse command line arguments
    Config cfg;
    auto graph = load_graph(cfg, argc, argv);

    // Precompute offline optimal congestion if needed
    offlineOptimal(graph, cfg);
/*
    RoutingEngine engine;
    for (SolverType type : cfg.solvers ) {
        auto result = engine.solve(*graph, cfg, type);
        if (result) {
            for (const auto& [str, _type] : SOLVER_MAP) {
                if ( _type == type ) {
                    result->storeAsFile(str + "_result.json", OutPutFormat::JASON);
                    break;
                }
            }
        }
    }
    return 0;
    */

    std::cout <<"Sequential Electrical MWU:\n";
    ParElectricalFlowMWU<SequentialExecution> sequential(*graph, 0, true, SequentialExecution{});
    auto seq_scheme = sequential.solve();

    seq_scheme->printRoutingTable();

    std::cout <<"OMP Electrical MWU:\n";
    ParElectricalFlowMWU<OpenMPExecution> omp(*graph, 0, true, OpenMPExecution{});
    auto omp_scheme = omp.solve();

    omp_scheme->printRoutingTable();
}
