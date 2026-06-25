#include "../include/algorithms/oblivious/mwu/tree_mwu.h"
#include "../include/io/parse_argurment_io.h"
#include "../include/routing/routing_engine.h"
#include "../include/core/utils.h"

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
    IRoutingResult result;
    for (SolverType type : cfg.solvers ) {

        // Toy example of including try & catch to handle exceptions
        try {
            auto res = engine.solve(*graph, cfg, type);
            if (res) {
                result = std::move(res.value());
            } else {
                throw std::runtime_error("Failed to solve for solver type: " + getSolverName(type));
            }
            if (!parser.write(result, "result_"+getSolverName(type)+".json", OutputFormat::JSON)) {
                throw std::runtime_error("Failed to write results to file for solver type: " + getSolverName(type));
            }
        } catch(std::runtime_error& e) {
            std::cerr << e.what() << std::endl;
            return -1;
        } catch (std::invalid_argument& e) {
            std::cerr << e.what() << std::endl;
            return -1;
        }

    }



    return 0;
}
