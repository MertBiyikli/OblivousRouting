//
// Created by Mert on 29.10.26.
//
// Independent Catch2 translation unit for N-1 analysis.
// Uses an in-memory, 7-vertex/9-link topology with the same edge IDs and
// endpoints as Backbone 1221. Capacities are synthetic: no LGF fixture needed.
// Add this .cpp to the integration_tests target.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "analysis/resilience/resilience_analyzer.h"
#include "analysis/resilience/resilience_result_json.h"
#include "data_structures/graph/graph.h"

#include <cstdint>
#include <set>
#include <vector>

namespace {
using Graph = optimized::Graph<EdgeData>;

Graph makeN1TestGraph() {
    // Input order defines canonical stable physical IDs 0,2,...,16.
    const std::vector<Graph::InputEdge> edges{
        {0, 2, {10.0, 1.0}}, // edge  0
        {0, 4, {10.0, 1.0}}, // edge  2
        {0, 5, {10.0, 1.0}}, // edge  4
        {1, 3, {10.0, 1.0}}, // edge  6
        {1, 4, {10.0, 1.0}}, // edge  8
        {1, 6, {10.0, 1.0}}, // edge 10
        {2, 5, {10.0, 1.0}}, // edge 12
        {3, 4, {10.0, 1.0}}, // edge 14
        {4, 6, {10.0, 1.0}}  // edge 16
    };
    return Graph(7, edges);
}

std::set<int> physicalIds(const NetworkData& state) {
    std::set<int> ids;
    for (const auto& link : state.links) {
        REQUIRE(link.edge_id >= 0);
        REQUIRE(ids.insert(link.edge_id).second);
        REQUIRE(link.source >= 0);
        REQUIRE(link.target >= 0);
        REQUIRE(link.source != link.target);
    }

    return ids;
}
} // namespace

TEST_CASE("N-1 prepared execution records independent link failures", "[n1]") {
    Graph graph = makeN1TestGraph();
    const auto fingerprint = optimized::GraphFingerprint::compute(graph);
    const std::set<int> expected{0, 2, 4, 6, 8, 10, 12, 14, 16};

    Config cfg{};
    cfg.solvers = {SolverType::ELECTRICAL_SKETCHING};
    cfg.demand_models = {DemandModelType::UNIFORM};
    cfg.demand_scale = 1.0;

    REQUIRE(graph.getNumUndirectedEdges() == 9);
    auto computed = analyzeN1(graph, cfg, ExecutionMode::FAST);
    REQUIRE(computed);
    const auto& result = computed.value();

    REQUIRE(result.graph_fingerprint == fingerprint);
    REQUIRE(result.outcomes.size() == expected.size());
    REQUIRE(physicalIds(result.baseline) == expected);

    std::set<int> observed;
    std::size_t successful = 0;
    std::size_t disconnected = 0;
    std::size_t failed = 0;

    for (const auto& outcome : result.outcomes) {
        REQUIRE(observed.insert(outcome.edge_id).second);
        REQUIRE(expected.count(outcome.edge_id) == 1);
        REQUIRE(graph.edgeAlive(outcome.edge_id));
        const auto endpoints = graph.getEdgeEndpoints(outcome.edge_id);
        REQUIRE(outcome.source == endpoints.first);
        REQUIRE(outcome.target == endpoints.second);

        switch (outcome.status) {
        case FailureStatus::SUCCESS: {
            ++successful;
            REQUIRE(outcome.reoptimized.has_value());
            REQUIRE(outcome.comparison.has_value());
            const auto& after = *outcome.reoptimized;
            auto survivors = expected;
            survivors.erase(outcome.edge_id);
            REQUIRE(physicalIds(after) == survivors);
            REQUIRE(after.total_demand_volume == Catch::Approx(result.baseline.total_demand_volume));
            REQUIRE(outcome.comparison->congestion_before == Catch::Approx(result.baseline.max_congestion));
            REQUIRE(outcome.comparison->congestion_after == Catch::Approx(after.max_congestion));
            break;
        }
        case FailureStatus::DISCONNECTED:
            ++disconnected;
            REQUIRE_FALSE(outcome.reoptimized.has_value());
            break;
        case FailureStatus::EXECUTION_FAILED:
            ++failed;
            REQUIRE_FALSE(outcome.reoptimized.has_value());
            break;
        }

        auto written = ResilienceResultJsonWriter::write(result, "results/scenarios/n1-demo.json");

        REQUIRE(written);
    }

    REQUIRE(observed == expected);
    REQUIRE(successful + disconnected + failed == expected.size());
    REQUIRE(graph.getNumUndirectedEdges() == expected.size());
    REQUIRE(optimized::GraphFingerprint::compute(graph) == fingerprint);
    for (int id : expected)
        REQUIRE(graph.edgeAlive(id));
}
