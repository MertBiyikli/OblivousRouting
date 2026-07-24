#include "../include/routing/routing_engine.h"
#include "../include/algorithms/oblivious/flow_sparsifier/boundary_linked_exp_decomp/expander_decomp.h"
#include "algorithms/oblivious/flow_sparsifier/boundary_linked_exp_decomp/practical_adapter.h"
#include "algorithms/oblivious/flow_sparsifier/boundary_linked_exp_decomp/validation.h"
#include <cassert>

#include "algorithms/oblivious/mwu/flow_sparsifier_mwu.h"

int main(int argc, char **argv) {

    // Parse command line arguments
    RoutingEngine engine;
    auto res = engine.entry(argc, argv);
    if (!res) {
        std::cerr << "Error: " << getError(res).error().message << std::endl;
        return 1;
    }


/*
    auto graph_p = engine.getGraph();
    auto graph = Graph(*graph_p);
    std::vector<int> all_vertices(graph.node_count());
    std::iota(all_vertices.begin(), all_vertices.end(), 0);


    const Capacity root_volume = graph.volume(all_vertices);
    const double log_volume = std::log10(
        static_cast<double>(std::max<Capacity>(2, root_volume)));

    const double gamma_cmp =
        22.0 + std::ceil(5.0 * log_volume * log_volume);

    DecompositionConfig decomposition_config;
    const double log2_volume = std::max(
    1.0,
    std::log2(static_cast<double>(
        std::max<Capacity>(2, root_volume))));

    const double alpha_max =1.0 /(4.0 * gamma_cmp * log2_volume * log2_volume);

    decomposition_config.alpha = 0.9 * alpha_max;
    decomposition_config.phi = 0.2;
    decomposition_config.gamma_cmp = gamma_cmp;


    constexpr std::size_t exact_cluster_limit = 18;
    PracticalOracleConfig oracle_config;
    Capacity total_multiplicity = 0;
    for (const Edge& edge : graph.edges()) {
        if (total_multiplicity >
            std::numeric_limits<Capacity>::max() - edge.multiplicity) {
            throw std::overflow_error("terminal upper bound overflow");
            }
        total_multiplicity += edge.multiplicity;
    }

    const Capacity maximum_boundary_copies =
        std::max<Capacity>(
            1,
            static_cast<Capacity>(std::ceil(
                decomposition_config.alpha / decomposition_config.phi)));

    if (total_multiplicity >
        std::numeric_limits<Capacity>::max() / maximum_boundary_copies) {
        throw std::overflow_error("terminal upper bound overflow");
        }

    const Capacity terminal_upper_bound =
        total_multiplicity * maximum_boundary_copies;

    if (terminal_upper_bound >
        std::numeric_limits<std::size_t>::max()) {
        throw std::overflow_error("terminal upper bound exceeds size_t");
        }

    oracle_config.maximum_explicit_terminals =
        static_cast<std::size_t>(terminal_upper_bound);
    oracle_config.exact_fallback_max_vertices = exact_cluster_limit;
    oracle_config.gamma_cmp = gamma_cmp;
    oracle_config.allow_exact_fallback = true;
    oracle_config.selection_mode = OracleSelectionMode::Auto;


    const auto oracle = std::make_shared<PracticalCutMatchingOracle>(oracle_config);
    auto t0 = timeNow();
    const DecompositionResult result = BoundaryLinkedDecomposer(decomposition_config, oracle).decompose(graph);
    std::cout << "Running time boundary-linked: " << duration(timeNow()-t0) << " microseconds." << std::endl;

    assert(result.boundary_accounting_identity_verified && "end-to-end boundary accounting failed");

    ValidationOptions validation_options;
    validation_options.exact_max_vertices = exact_cluster_limit;
    validation_options.require_exact_for_all = false;
    const ValidationReport validation = validate_decomposition(graph, decomposition_config, result, {}, validation_options);
    for (const std::string& error : validation.errors) {
        std::cerr << "Boundary-linked validation: " << error << '\n';
    }
    assert(validation.valid() && "end-to-end result failed independent validation");

    std::vector<bool> seen(graph.node_count(), false);
    const ExactConductanceOracle exact(exact_cluster_limit);
    for (const OutputCluster& cluster : result.clusters) {
        for (int v : cluster.vertices) {
            assert(!seen[v] &&  "end-to-end clusters overlap");
            seen[v] = true;
        }
        if (cluster.vertices.size() <= exact_cluster_limit) {
            const OracleResult check = exact.analyze(
                graph,
                cluster.vertices,
                decomposition_config.alpha / cluster.phi,
                cluster.phi);
            assert(check.kind == OracleResult::Kind::ExpanderCertificate &&
                "end-to-end final cluster failed exact boundary-linked verification");
        }
    }
    for (bool vertex_was_seen : seen) {
        assert(vertex_was_seen && "end-to-end decomposition omitted a vertex");
    }

    std::cout
    << "Clusters: " << result.clusters.size() << '\n'
    << "Splits: " << result.splits.size() << '\n'
    << "Randomized certificates: "
    << result.randomized_expansion_proof_count << '\n'
    << "Boundary accounting: "
    << result.boundary_accounting_identity_verified << '\n'
    << "Output expansion certified: "
    << result.output_expansion_certified << '\n';

    for (std::size_t i = 0; i < result.clusters.size(); ++i) {
        std::cout << "Cluster " << i << ": ";
        for (int v : result.clusters[i].vertices) {
            std::cout << v << ' ';
        }
        std::cout << '\n';
    }*/
}
