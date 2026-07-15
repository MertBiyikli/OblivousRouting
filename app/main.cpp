#include "../include/io/parse_argument_io.h"
#include "../include/routing/routing_engine.h"

#include "algorithms/semi_oblivious/expander_hierarchy/tree_flow_electrical_embedder.h"
#include "algorithms/semi_oblivious/expander_hierarchy/tree_flow_router.h"
#include "algorithms/semi_oblivious/expander_hierarchy/tree_sparsifier.h"
#include "algorithms/semi_oblivious/expander_hierarchy/preprocessing/hierarchy_preprocessor.h"
#include "routing/routing_runner.h"

/*
LinearRoutingTable runExpanderHierarchySmokeTest(IGraph& graph) {
    //graph.print();
    if (graph.getNumNodes() < 2) {
        std::cerr
            << "[expander smoke] Graph requires at least two vertices.\n";
    }

    std::cout
        << "[expander smoke] preprocessing graph with "
        << graph.getNumNodes() << " vertices and "
        << graph.getNumUndirectedEdges() << " edges\n";

    XCutHierarchyPreprocessor preprocessor;

    auto hierarchy_result = preprocessor.build(graph);

    if (!hierarchy_result) {
        std::cerr
            << "[expander smoke] preprocessing failed\n";
        std::cerr << hierarchy_result.error().message << '\n';
    }

    const HierarchyResult& hierarchy = *hierarchy_result;

    TreeSparsifierBuilder tree_builder;

    auto tree_result =
        tree_builder.build(graph, hierarchy);

    if (!tree_result) {
        std::cerr
            << "[tree sparsifier] construction failed\n"
            << tree_result.error().message
            << '\n';

    }

    const TreeSparsifier& tree = *tree_result;
    //tree.print();


    const int root = 0;
    //const int target = graph.getNumNodes() - 3;
    TreeFlowRouter tree_router(tree);
    constexpr double demand_value = 1.0;
    LinearRoutingTable linear_routing_table;
    linear_routing_table.init(graph);

    for (int target = 0; target < graph.getNumNodes(); ++target) {
        if (target == root) {
            continue;
        }

        auto tree_flow_result = tree_router.routePair(
            root,
            target,
            demand_value
        );


        if (!tree_flow_result) {
            std::cerr
                << "[tree routing] failed\n"
                << tree_flow_result.error().message
                << '\n';

        }

        const TreeFlowResult& tree_flow = *tree_flow_result;

        TreeFlowElectricalEmbedder embedder(
        graph,
        hierarchy,
        tree
    );

        auto embedding_result = embedder.embed(tree_flow);

        if (!embedding_result) {
            std::cerr
                << "[electrical embedding] failed\n"
                << embedding_result.error().message
                << '\n';
        }

        const auto& embedding = *embedding_result;


        for (int e = 0; e < graph.getNumDirectedEdges();++e) {
            const double flow = embedding.signed_edge_flow[e];

            if (std::abs(flow) <= 1e-10) {
                continue;
            }

            const auto [u, v] = graph.getEdgeEndpoints(e);

            // get sign of the flow
            if (flow < 0) {
                int anti_e = graph.getAntiEdge(e);
                linear_routing_table.addFlow(anti_e, target, std::abs(flow));
            }else {
                linear_routing_table.addFlow(e, target, flow);
            }
        }
    }
    return linear_routing_table;
}
*/

int main(int argc, char **argv) {
    // Parse command line arguments

    RoutingEngine engine;
    auto result = engine.entry(argc, argv);
    if (!result) {
        std::cerr << "Error: " << result.error().message << std::endl;
        return 1;
    }
/*
    IRoutingResult res;
    if (auto graph = engine.getGraph()) {
        auto t0 = timeNow();
        auto table = runExpanderHierarchySmokeTest(*graph);
        std::cout << "Running time expander: " << (duration( timeNow()-t0)) << " microseconds" << std::endl;
        const std::unique_ptr<RoutingScheme> scheme = std::make_unique<LinearRoutingScheme>(*graph, 0, std::move(table));


        auto dem = DemandEvaluator::evaluate(*graph, scheme ,engine.getConfig(), res);

        for (std::size_t i = 0; i < res.demand_evaluations.size(); ++i) {
            const auto& eval = res.demand_evaluations[i];

            std::cout << "    {\n";
            std::cout << "      \"demand_model\": \""
                << (demandModelName(eval.demand_type)) << "\",\n";
            std::cout << "      \"congestion\": "
                << eval.congestion << ",\n";
            std::cout << "      \"runtime_microseconds\": "
                << eval.runtime_microseconds << "\n";
            std::cout << "    }";

            if (i + 1 < res.demand_evaluations.size()) {
                std::cout << ",";
            }

            std::cout << "\n";
        }

    }
    */
}
