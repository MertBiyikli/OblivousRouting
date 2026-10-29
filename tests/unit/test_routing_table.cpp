#include "../../include/routing/routing_table.h"
#include "data_structures/graph/graph.h"
#include "routing/storage/allpair_routing_table.h"
#include "routing/storage/linear_routing_table.h"
#include "routing/storage/routing_scheme_verifier.h"
#include "routing/storage/routing_table_adapter.h"
#include "utils/my_math.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <string>
#include <vector>

using Catch::Approx;

static optimized::Graph<EdgeData> createSimpleGraph() {
    std::vector<optimized::Graph<EdgeData>::InputEdge> input = {
        {0, 1, EdgeData{1.0, 1.0}},
        {1, 2, EdgeData{1.0, 1.0}},
        {2, 3, EdgeData{1.0, 1.0}},
        {1, 3, EdgeData{1.0, 1.0}},
    };
    return optimized::Graph<EdgeData>(4, input);
}

TEST_CASE("AllPairRoutingTable - Construction and Initialization", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    AllPairRoutingTable table;

    table.init(graph);

    REQUIRE(table.n == 4);
    REQUIRE(table.adj_ids.size() == 8);  // 4 undirected edges = 8 directed edges
    REQUIRE(table.adj_vals.size() == 8);
}

TEST_CASE("AllPairRoutingTable - Add Flow by Nodes", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    AllPairRoutingTable table;

    table.init(graph);

    // Add flow for commodity (0, 3) on edge 0
    table.addFlow(0, 0, 3, 0.5);

    REQUIRE(!table.adj_ids[0].empty());
    REQUIRE(!table.adj_vals[0].empty());
}

TEST_CASE("AllPairRoutingTable - Get Flow by Nodes", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    AllPairRoutingTable table;

    table.init(graph);

    // Add flow for commodity (0, 3) on edge 0
    table.addFlow(0, 0, 3, 0.5);

    double f = table.getFlow(0, 0, 3);
    REQUIRE(f == Approx(0.5));

}


TEST_CASE("AllPairRoutingTable - Multiple Flows on Same Edge", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    AllPairRoutingTable table;

    table.init(graph);

    // Add multiple flows to the same edge
    table.addFlow(0, 0, 1, 0.3);
    table.addFlow(0, 0, 2, 0.2);
    table.addFlow(0, 1, 3, 0.5);

    REQUIRE(table.adj_ids[0].size() == 3);
    REQUIRE(table.adj_vals[0].size() == 3);

    REQUIRE(table.getFlow(0, 0, 1) == Approx(0.3));
    REQUIRE(table.getFlow(0, 0, 2) == Approx(0.2));
    REQUIRE(table.getFlow(0, 1, 3) == Approx(0.5));
}



TEST_CASE("AllPairRoutingTable - Multiple Flows Different Edges", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    AllPairRoutingTable table;

    table.init(graph);

    // Add flows to different edges
    table.addFlow(0, 0, 1, 0.3);
    table.addFlow(1, 0, 1, 0.2);
    table.addFlow(2, 0, 1, 0.5);

    REQUIRE(table.getFlow(0, 0, 1) == Approx(0.3));
    REQUIRE(table.getFlow(1, 0, 1) == Approx(0.2));
    REQUIRE(table.getFlow(2, 0, 1) == Approx(0.5));
}

TEST_CASE("LinearRoutingTable - Construction and Initialization", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    LinearRoutingTable table;

    table.init(graph);

    REQUIRE(table.n == 4);
    REQUIRE(table.src_ids.size() == 8);
    REQUIRE(table.src_flows.size() == 8);
}

TEST_CASE("LinearRoutingTable - Add Flow", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    LinearRoutingTable table;

    table.init(graph);

    // Add source flow for commodity from source 0
    table.addFlow(0, 0, 0.5);

    REQUIRE(!table.src_ids[0].empty());
    REQUIRE(!table.src_flows[0].empty());
}

TEST_CASE("LinearRoutingTable - Get Flow", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    LinearRoutingTable table;

    table.init(graph);

    table.addFlow(0, 0, 0.5);
    double flow = table.getFlow(0, 0);

    REQUIRE(flow == Approx(0.5));
}

TEST_CASE("LinearRoutingTable - Multiple Sources Same Edge", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    LinearRoutingTable table;

    table.init(graph);

    // Add flows from multiple sources on same edge
    table.addFlow(0, 0, 0.3);
    table.addFlow(0, 1, 0.2);
    table.addFlow(0, 2, 0.5);

    REQUIRE(table.src_ids[0].size() == 3);
    REQUIRE(table.getFlow(0, 0) == Approx(0.3));
    REQUIRE(table.getFlow(0, 1) == Approx(0.2));
    REQUIRE(table.getFlow(0, 2) == Approx(0.5));
}

TEST_CASE("LinearRoutingTable - Erase Flow", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    LinearRoutingTable table;

    table.init(graph);

    table.addFlow(0, 0, 0.5);
    table.addFlow(0, 1, 0.3);

    int initial_size = table.src_ids[0].size();

    // Erase flow from source 0
    table.eraseFlow(0, 0);

    REQUIRE(table.src_ids[0].size() <= initial_size);
}

TEST_CASE("LinearRoutingTable - Non-existent Flow", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    LinearRoutingTable table;

    table.init(graph);

    // Try to get flow that wasn't added
    double flow = table.getFlow(0, 3);

    REQUIRE(flow == Approx(0.0));
}


TEST_CASE("LinearRoutingTable - Multiple Edges Different Sources", "[RoutingTable]") {
    auto graph = createSimpleGraph();
    LinearRoutingTable table;

    table.init(graph);

    // Different sources on different edges
    table.addFlow(0, 0, 0.3);
    table.addFlow(1, 1, 0.2);
    table.addFlow(2, 2, 0.5);

    REQUIRE(table.getFlow(0, 0) == Approx(0.3));
    REQUIRE(table.getFlow(1, 1) == Approx(0.2));
    REQUIRE(table.getFlow(2, 2) == Approx(0.5));
}




namespace fs = std::filesystem;

optimized::Graph<EdgeData> makeGraph(double capacity = 10.0) {
    std::vector<optimized::Graph<EdgeData>::InputEdge> edges = {
        {0, 1, EdgeData{capacity, 1.0}},
        {1, 2, EdgeData{10.0, 1.0}}
    };

    return optimized::Graph<EdgeData>(3, edges);
}

struct TemporaryRoutingFile {
    fs::path path;

    explicit TemporaryRoutingFile(std::string name)
        : path(fs::temp_directory_path() / std::move(name)) {}

    ~TemporaryRoutingFile() {
        std::error_code ec;
        fs::remove(path, ec);
    }
};

TEST_CASE("Routing table rejects changed graph capacity", "[routing][persistence]") {
    auto original = makeGraph(10.0);
    auto modified = makeGraph(20.0);

    const auto originalFingerprint = optimized::GraphFingerprint::compute(original);
    const auto modifiedFingerprint = optimized::GraphFingerprint::compute(modified);

    REQUIRE(originalFingerprint != modifiedFingerprint);

    LinearRoutingTable table;
    table.init(original);
    table.addFlow(original.edgeId(1, 0), 1, 1.0);
    table.addFlow(original.edgeId(2, 1), 2, 1.0);
    table.addFlow(original.edgeId(1, 0), 2, 1.0);

    TemporaryRoutingFile file("erouting_wrong_graph_test.bin");

    auto saved = RoutingTableAdapter::saveLinear(table, original, file.path.string(), originalFingerprint, 0);
    REQUIRE(saved.has_value());

    std::int32_t root = -1;
    auto loaded = RoutingTableAdapter::loadLinear(modified, file.path.string(), modifiedFingerprint, &root);

    REQUIRE_FALSE(loaded.has_value());
}

TEST_CASE("Routing table rejects truncated binary", "[routing][persistence]") {
    auto graph = makeGraph();
    const auto fingerprint = optimized::GraphFingerprint::compute(graph);

    LinearRoutingTable table;
    table.init(graph);
    table.addFlow(graph.edgeId(1, 0), 1, 1.0);

    TemporaryRoutingFile file("erouting_truncated_integration_test.bin");

    auto saved = RoutingTableAdapter::saveLinear(table, graph, file.path.string(), fingerprint, 0);
    REQUIRE(saved.has_value());

    const auto originalSize = fs::file_size(file.path);
    REQUIRE(originalSize > 44);

    fs::resize_file(file.path, originalSize - 1);

    std::int32_t root = -1;
    auto loaded = RoutingTableAdapter::loadLinear(graph, file.path.string(), fingerprint, &root);

    REQUIRE_FALSE(loaded.has_value());
}

TEST_CASE("All-pairs scheme survives binary round trip", "[routing][persistence]") {
    auto graph = makeGraph();
    const auto fingerprint = optimized::GraphFingerprint::compute(graph);

    AllPairRoutingTable table;
    table.init(graph);

    const int edge01 = graph.edgeId(0, 1);
    const int edge12 = graph.edgeId(1, 2);

    // Unit flow for 0 -> 1.
    table.addFlow(edge01, 0, 1, 1.0);

    // Unit flow for 0 -> 2.
    table.addFlow(edge01, 0, 2, 1.0);
    table.addFlow(edge12, 0, 2, 1.0);

    // Unit flow for 1 -> 2.
    table.addFlow(edge12, 1, 2, 1.0);

    AllPairRoutingScheme original(graph, std::move(table));

    TemporaryRoutingFile file("erouting_allpairs_integration_test.bin");

    auto saved = RoutingTableAdapter::saveAllPair(original.routing_table, graph, file.path.string(), fingerprint);
    REQUIRE(saved.has_value());

    auto loaded = RoutingTableAdapter::loadAllPair(graph, file.path.string(), fingerprint);
    REQUIRE(loaded.has_value());

    AllPairRoutingScheme restored(graph, std::move(loaded.value()));

    for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
        for (int s = 0; s < graph.getNumNodes(); ++s) {
            for (int t = 0; t < graph.getNumNodes(); ++t) {
                if (s == t) continue;

                CHECK(restored.getFlow(e, s, t) == Catch::Approx(original.getFlow(e, s, t)));
            }
        }
    }

    auto verified = RoutingSchemeVerifier::verify(graph, restored);
    REQUIRE(verified.has_value());
}


TEST_CASE(
    "LinearRoutingScheme - Optimized demand routing matches reference",
    "[RoutingTable][CombinedRouting]"
) {
    auto graph = createSimpleGraph();

    LinearRoutingTable table;
    table.init(graph);

    // Populate deterministic routing coefficients.
    // These are test coefficients, not a complete valid routing scheme.
    for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
        for (int s = 0; s < graph.getNumNodes(); ++s) {
            const double coefficient =
                0.01 * static_cast<double>((e + 1) * (s + 2));

            table.addFlow(e, s, coefficient);
        }
    }

    LinearRoutingScheme scheme(graph, 0, std::move(table));

    SECTION("Asymmetric demand matrix") {
        demands demandMap;

        demandMap.addDemand(0, 1, 100.0);
        demandMap.addDemand(1, 0, 300.0);
        demandMap.addDemand(0, 2, 75.0);
        demandMap.addDemand(2, 3, 125.0);
        demandMap.addDemand(3, 1, 50.0);

        std::vector<double> actual;
        scheme.routeDemands(actual, demandMap);

        REQUIRE(actual.size() ==
                static_cast<std::size_t>(graph.getNumDirectedEdges()));

        for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
            const int reverse = graph.reverse(e).id;

            if (e > reverse) {
                continue;
            }

            std::vector<double> basis;
            scheme.getEdgeBasis(e, basis);

            double referenceLoad = 0.0;

            for (std::size_t i = 0; i < demandMap.size(); ++i) {
                const int s = demandMap.source[i];
                const int t = demandMap.target[i];
                const double value = demandMap.demand_values[i];

                if (s == t || value <= 0.0) {
                    continue;
                }

                const double magnitude = std::abs(
                    basis[s] - basis[t]
                );

                if (magnitude < EPS) {
                    continue;
                }

                referenceLoad += value * magnitude;
            }

            const double expectedForward =
                referenceLoad / graph.edgeData(e).capacity;

            const double expectedReverse =
                referenceLoad / graph.edgeData(reverse).capacity;

            CHECK(actual[e] == Catch::Approx(expectedForward).epsilon(1e-10));
            CHECK(actual[reverse] == Catch::Approx(expectedReverse).epsilon(1e-10));
        }
    }

    SECTION("Sparse and zero-valued demands") {
        demands demandMap;

        demandMap.addDemand(0, 3, 50.0);
        demandMap.addDemand(3, 0, 0.0);
        demandMap.addDemand(1, 2, 0.0);
        demandMap.addDemand(2, 1, 25.0);

        std::vector<double> actual;
        scheme.routeDemands(actual, demandMap);

        REQUIRE(actual.size() ==
                static_cast<std::size_t>(graph.getNumDirectedEdges()));

        for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
            const int reverse = graph.reverse(e).id;

            if (e > reverse) {
                continue;
            }

            std::vector<double> basis;
            scheme.getEdgeBasis(e, basis);

            double referenceLoad = 0.0;

            for (std::size_t i = 0; i < demandMap.size(); ++i) {
                const int s = demandMap.source[i];
                const int t = demandMap.target[i];
                const double value = demandMap.demand_values[i];

                if (s == t || value <= 0.0) {
                    continue;
                }

                const double magnitude = std::abs(
                    basis[s] - basis[t]
                );

                if (magnitude < EPS) {
                    continue;
                }

                referenceLoad += value * magnitude;
            }

            CHECK(actual[e] ==
                  Catch::Approx(referenceLoad / graph.edgeData(e).capacity)
                      .epsilon(1e-10));

            CHECK(actual[reverse] ==
                  Catch::Approx(referenceLoad / graph.edgeData(reverse).capacity)
                      .epsilon(1e-10));
        }
    }
}
