//
// Created by Mert on 28.10.26.
//

#include "analysis/failures/link_failures.h"
#include "analysis/failures/network_connectivity.h"
#include "utils.h"

#include <catch2/catch_test_macros.hpp>

TEST_CASE("Connected graph remains connected after non-bridge failure", "[ScenarioAnalysis][Connectivity]") {
    auto graph = createScenarioTestGraph();

    REQUIRE(isConnected(graph));

    const auto checkpoint = graph.checkpoint();

    // In the existing four-node fixture:
    // 0-1, 1-2, 2-3, 0-3.
    // Removing 1-2 does not disconnect the graph.
    REQUIRE(graph.removeEdge(2));

    CHECK(isConnected(graph));

    REQUIRE(graph.rollback(checkpoint).has_value());

    CHECK(isConnected(graph));
}

TEST_CASE("Link failure rejects demand across disconnected components", "[ScenarioAnalysis][Connectivity]") {
    auto graph = createScenarioTestGraph();

    const auto originalFingerprint = optimized::GraphFingerprint::compute(graph);

    const auto checkpoint = graph.checkpoint();

    // Remove two physical links so node 1
    // becomes isolated: 0-1 and 1-2.
    REQUIRE(graph.removeEdge(0));
    REQUIRE(graph.removeEdge(2));

    CHECK_FALSE(isConnected(graph));

    demands demand;
    demand.addDemand(1, 3, 10.0);

    const auto result = verifyDemandConnectivity(graph, demand);

    CHECK_FALSE(result.has_value());

    REQUIRE(graph.rollback(checkpoint).has_value());

    CHECK(optimized::GraphFingerprint::compute(graph) == originalFingerprint);

    CHECK(graph.edgeAlive(0));
    CHECK(graph.edgeAlive(2));
}

TEST_CASE("Disconnected components without crossing demand are allowed", "[ScenarioAnalysis][Connectivity]") {
    auto graph = createScenarioTestGraph();

    const auto checkpoint = graph.checkpoint();

    REQUIRE(graph.removeEdge(0));
    REQUIRE(graph.removeEdge(2));

    demands demand;

    // Nodes 0 and 3 remain connected.
    demand.addDemand(0, 3, 10.0);

    const auto result = verifyDemandConnectivity(graph, demand);

    CHECK(result.has_value());

    REQUIRE(graph.rollback(checkpoint).has_value());
}

TEST_CASE("LinkFailure rejects out-of-range edge ID", "[ScenarioAnalysis][Validation]") {
    auto graph = createScenarioTestGraph();

    const auto before = optimized::GraphFingerprint::compute(graph);

    LinkFailure scenario;
    scenario.linkFailure(999999);

    const auto result = scenario.verify(graph);

    REQUIRE_FALSE(result.has_value());

    CHECK(optimized::GraphFingerprint::compute(graph) == before);
}