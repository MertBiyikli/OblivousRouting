#include "../include/algorithms/mwu/tree_mwu.h"
#include "../include/io/parse_argurment_io.h"
#include "../include/core/routing_engine.h"
#include "../include/core/utils.h"
#include "../include/algorithms/semi_oblivious/or_tools_linear_optimizer.h"

int main(int argc, char **argv) {
    // Parse command line arguments
    Config cfg;
    auto graph = load_graph(cfg, argc, argv);

    offlineOptimal(graph, cfg);
    for (auto [type, cong] : cfg.offline_opt_per_model) {
        std::cout << "Offline opt ["+type+"]: " << cong << std::endl;
    }

    RoutingResultWriter parser;
    RoutingEngine engine;
    RoutingRunResult result;
    for (SolverType type : cfg.solvers ) {
        auto res = engine.solve(*graph, cfg, type);
        if (res) {
            result = std::move(res.value());
        } else {
            std::cerr << "Failed to solve for solver type: " << getSolverName(type) << std::endl;
            return -1;
        }
        if (!parser.write(result, "result_"+getSolverName(type)+".json", OutputFormat::JSON)) {
            throw std::runtime_error("Failed to write results to file for solver type: " + getSolverName(type));
        }
    }



    return 0;
}
