#include "../include/io/cli.h"
#include "../include/routing/routing_engine.h"
#include "routing/storage/routing_table_adapter.h"

int main(int argc, char** argv) {

    // Parse command line arguments
    auto cmd = parse(argc, argv);
    if (!cmd) {
        std::cerr << "Error: " << getError(cmd).error().message << std::endl;
        return 1;
    }

    switch (cmd->action) {

    case Action::Help:
        std::cout << help(argv[0]);
        return 0;

    case Action::ListSolvers:
        std::cout << solverList();
        return 0;

    case Action::Run:
        break;
    }

    // Running engine
    auto cfg = cmd.value().config;
    RoutingEngine engine;
    auto res = engine.entry(cfg);
    if (!res) {
        std::cerr << "Error: " << getError(res).error().message << std::endl;
        return 1;
    }

    // Access the graph without transferring ownership.
    auto graphResult = engine.getGraph();
    if (!graphResult) {
        std::cerr << "Error: " << getError(graphResult).error().message << '\n';
        return 1;
    }

    auto& graph = *graphResult.value();
    constexpr std::uint64_t graphFingerprint = 0; // Temporary test placeholder.
    constexpr double tolerance = 1e-10;

    for (std::size_t i = 0; i < res.value().solver_results.size(); ++i) {
        auto& solverResult = res.value().solver_results[i];
        auto* originalScheme = solverResult.scheme.get();

        std::unique_ptr<RoutingScheme> restoredScheme;
        const std::string path = "routing_table_" + std::to_string(i) + ".bin";

        if (auto* linear = dynamic_cast<LinearRoutingScheme*>(originalScheme)) {
            auto saved = RoutingTableAdapter::saveLinear(linear->routing_table, graph, path, graphFingerprint);
            if (!saved) {
                std::cerr << getError(saved).error().message << '\n';
                return 1;
            }

            auto loaded = RoutingTableAdapter::loadLinear(graph, path, graphFingerprint);
            if (!loaded) {
                std::cerr << getError(loaded).error().message << '\n';
                return 1;
            }

            restoredScheme = std::make_unique<LinearRoutingScheme>(graph, linear->root_x, std::move(loaded.value()));
        } else if (auto* allPairs = dynamic_cast<AllPairRoutingScheme*>(originalScheme)) {
            auto saved = RoutingTableAdapter::saveAllPair(allPairs->routing_table, graph, path, graphFingerprint);
            if (!saved) {
                std::cerr << getError(saved).error().message << '\n';
                return 1;
            }

            auto loaded = RoutingTableAdapter::loadAllPair(graph, path, graphFingerprint);
            if (!loaded) {
                std::cerr << getError(loaded).error().message << '\n';
                return 1;
            }

            restoredScheme = std::make_unique<AllPairRoutingScheme>(graph, std::move(loaded.value()));
        } else {
            std::cout << "Skipping unsupported scheme at index " << i << '\n';
            continue;
        }

        double maxError = 0.0;

        for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
            for (int s = 0; s < graph.getNumNodes(); ++s) {
                for (int t = 0; t < graph.getNumNodes(); ++t) {
                    if (s == t)
                        continue;

                    const double original = originalScheme->getFlow(e, s, t);
                    const double restored = restoredScheme->getFlow(e, s, t);

                    if (!std::isfinite(original) || !std::isfinite(restored)) {
                        std::cerr << "Non-finite flow detected.\n";
                        return 1;
                    }

                    maxError = std::max(maxError, std::abs(original - restored));
                }
            }
        }

        std::cout << "Scheme " << i << " maximum round-trip flow error: " << maxError << '\n';

        if (maxError > tolerance) {
            std::cerr << "Round-trip verification failed.\n";
            return 1;
        }

        std::cout << "Scheme " << i << " round-trip verification passed.\n";
    }

    return 0; // something new
}
