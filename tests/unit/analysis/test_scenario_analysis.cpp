//
// Created by Mert on 22.10.26.
//

#include "routing/routing_result.h"
#include "analysis/scenario_analysis.h"

#include <absl/strings/str_format.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <map>

static optimized::Graph<EdgeData> createScenarioTestGraph() {
    std::vector<optimized::Graph<EdgeData>::InputEdge> edges = {{0, 1, EdgeData{10.0, 1.0}}, {1, 2, EdgeData{20.0, 1.0}}, {2, 3, EdgeData{30.0, 1.0}}, {0, 3, EdgeData{40.0, 1.0}}};

    return optimized::Graph<EdgeData>(4, edges);
}

TEST_CASE("CapacityUpgrade accepts valid upgrades", "[ScenarioAnalysis][CapacityUpgrade]") {

    auto graph = createScenarioTestGraph();
    CapacityUpgrade scenario(graph);

    REQUIRE(scenario.upgrade(0, 1.5).has_value());
    REQUIRE(scenario.upgrade(2, 2.0).has_value());
    REQUIRE(scenario.verify().has_value());

    std::map<int, double> upgrades(scenario.upgrades.begin(), scenario.upgrades.end());

    REQUIRE(upgrades.size() == 2);
    CHECK(upgrades.at(0) == Catch::Approx(1.5));
    CHECK(upgrades.at(2) == Catch::Approx(2.0));
}

TEST_CASE("CapacityUpgrade rejects invalid edge IDs", "[ScenarioAnalysis][CapacityUpgrade]") {

    auto graph = createScenarioTestGraph();
    CapacityUpgrade scenario(graph);

    const auto originalSize = scenario.upgrades.size();

    CHECK_FALSE(scenario.upgrade(-1, 1.5).has_value());
    CHECK_FALSE(scenario.upgrade(graph.getNumDirectedEdges(), 1.5).has_value());

    CHECK(scenario.upgrades.size() == originalSize);
}

TEST_CASE("CapacityUpgrade rejects invalid factors", "[ScenarioAnalysis][CapacityUpgrade]") {

    auto graph = createScenarioTestGraph();
    CapacityUpgrade scenario(graph);

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();

    CHECK_FALSE(scenario.upgrade(0, -1.0).has_value());
    CHECK_FALSE(scenario.upgrade(0, nan).has_value());
    CHECK_FALSE(scenario.upgrade(0, inf).has_value());

    CHECK(scenario.upgrades.empty());
}

TEST_CASE("CapacityUpgrade rejects duplicate physical links", "[ScenarioAnalysis][CapacityUpgrade]") {

    auto graph = createScenarioTestGraph();
    CapacityUpgrade scenario(graph);

    const int edge = 0;
    const int reverse = graph.reverse(edge).id;

    REQUIRE(scenario.upgrade(edge, 1.5).has_value());

    CHECK_FALSE(scenario.upgrade(edge, 2.0).has_value());
    CHECK_FALSE(scenario.upgrade(reverse, 2.0).has_value());

    CHECK(scenario.upgrades.size() == 1);
    CHECK(scenario.verify().has_value());
}

TEST_CASE("CapacityUpgrade does not modify baseline graph", "[ScenarioAnalysis][CapacityUpgrade]") {

    auto graph = createScenarioTestGraph();

    std::vector<double> originalCapacities;

    for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
        originalCapacities.push_back(graph.edgeData(e).capacity);
    }

    CapacityUpgrade scenario(graph);

    REQUIRE(scenario.upgrade(0, 1.5).has_value());
    REQUIRE(scenario.upgrade(2, 2.0).has_value());
    REQUIRE(scenario.verify().has_value());

    for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
        CHECK(graph.edgeData(e).capacity == Catch::Approx(originalCapacities[e]));
    }
}

TEST_CASE("CapacityUpgrade preserves valid state after rejected update", "[ScenarioAnalysis][CapacityUpgrade]") {

    auto graph = createScenarioTestGraph();
    CapacityUpgrade scenario(graph);

    REQUIRE(scenario.upgrade(0, 1.5).has_value());

    const auto before = scenario.upgrades;


    CHECK_FALSE(scenario.upgrade(-1, 2.0).has_value());
    CHECK_FALSE(scenario.upgrade(0, 3.0).has_value());

    REQUIRE(scenario.upgrades.size() == before.size());

    for (const auto& entry : before) {
        CHECK(scenario.upgrades.count(entry) == 1);
    }

    CHECK(scenario.verify().has_value());
}

TEST_CASE("LinkFailure accepts a single failure", "[ScenarioAnalysis][LinkFailure]") {
    auto graph = createScenarioTestGraph();
    LinkFailure scenario(graph);

    REQUIRE(scenario.linkFailure(0).has_value());
    REQUIRE(scenario.verify().has_value());
    REQUIRE(scenario.failed_links.size() == 1);
}

TEST_CASE("LinkFailure accepts multiple physical links", "[ScenarioAnalysis][LinkFailure]") {
    auto graph = createScenarioTestGraph();
    LinkFailure scenario(graph);

    int edge = graph.edgeId(0, 1);

    REQUIRE(scenario.linkFailure(edge).has_value());
    REQUIRE(scenario.verify().has_value());
    CHECK(scenario.failed_links.size() == 1);
}

TEST_CASE("LinkFailure rejects invalid edges", "[ScenarioAnalysis][LinkFailure]") {
    auto graph = createScenarioTestGraph();
    LinkFailure scenario(graph);

    CHECK_FALSE(scenario.linkFailure(-1).has_value());
    CHECK_FALSE(scenario.linkFailure(graph.getNumDirectedEdges()).has_value());
    CHECK(scenario.failed_links.empty());
}

TEST_CASE("LinkFailure rejects reverse-edge duplicates", "[ScenarioAnalysis][LinkFailure]") {
    auto graph = createScenarioTestGraph();
    LinkFailure scenario(graph);

    const int edge = 0;
    const int reverse = graph.reverse(edge).id;

    REQUIRE(scenario.linkFailure(edge).has_value());
    REQUIRE(scenario.linkFailure(reverse).has_value());
    CHECK(scenario.failed_links.size() == 1);
}

TEST_CASE("LinkFailure rejects invalid batches atomically", "[ScenarioAnalysis][LinkFailure]") {
    auto graph = createScenarioTestGraph();
    LinkFailure scenario(graph);

    REQUIRE(scenario.linkFailure(0).has_value());
    const auto before = scenario.failed_links;

    CHECK_FALSE(scenario.linkFailure(graph.getNumDirectedEdges()).has_value());
    CHECK_FALSE(scenario.linkFailure(-1).has_value());

    CHECK(scenario.failed_links == before);
    CHECK(scenario.verify().has_value());
}

TEST_CASE("LinkFailure preserves the baseline graph", "[ScenarioAnalysis][LinkFailure]") {
    auto graph = createScenarioTestGraph();
    const auto fingerprint = optimized::GraphFingerprint::compute(graph);

    LinkFailure scenario(graph);

    REQUIRE(scenario.linkFailure(0).has_value());
    REQUIRE(scenario.verify().has_value());

    CHECK(optimized::GraphFingerprint::compute(graph) == fingerprint);
}

TEST_CASE("TrafficGrowth accepts global scaling", "[ScenarioAnalysis][TrafficGrowth]") {
    auto graph = createScenarioTestGraph();
    TrafficGrowth scenario(graph);

    REQUIRE(scenario.setGlobalScale(1.5).has_value());
    REQUIRE(scenario.verify().has_value());

    CHECK(scenario.globalScale == Catch::Approx(1.5));
}

TEST_CASE("TrafficGrowth accepts commodity scaling", "[ScenarioAnalysis][TrafficGrowth]") {
    auto graph = createScenarioTestGraph();
    TrafficGrowth scenario(graph);

    REQUIRE(scenario.growPerCommodity(0, 1, 1.5).has_value());
    REQUIRE(scenario.growPerCommodity(1, 0, 2.0).has_value());
    REQUIRE(scenario.verify().has_value());

    CHECK(scenario.commodities.size() == 2);
}

TEST_CASE("TrafficGrowth rejects invalid scales", "[ScenarioAnalysis][TrafficGrowth]") {
    auto graph = createScenarioTestGraph();
    TrafficGrowth scenario(graph);

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();

    CHECK_FALSE(scenario.setGlobalScale(0.0).has_value());
    CHECK_FALSE(scenario.setGlobalScale(-1.0).has_value());
    CHECK_FALSE(scenario.setGlobalScale(nan).has_value());
    CHECK_FALSE(scenario.setGlobalScale(inf).has_value());

    CHECK(scenario.globalScale == Catch::Approx(1.0));
}

TEST_CASE("TrafficGrowth rejects invalid commodities", "[ScenarioAnalysis][TrafficGrowth]") {
    auto graph = createScenarioTestGraph();
    TrafficGrowth scenario(graph);

    CHECK_FALSE(scenario.growPerCommodity(-1, 1, 2.0).has_value());
    CHECK_FALSE(scenario.growPerCommodity(0, 0, 2.0).has_value());
    CHECK_FALSE(scenario.growPerCommodity(0, 1, 0.0).has_value());

    CHECK(scenario.commodities.empty());
}

TEST_CASE("TrafficGrowth rejects duplicate commodities", "[ScenarioAnalysis][TrafficGrowth]") {
    auto graph = createScenarioTestGraph();
    TrafficGrowth scenario(graph);

    REQUIRE(scenario.growPerCommodity(0, 1, 1.5).has_value());
    CHECK_FALSE(scenario.growPerCommodity(0, 1, 2.0).has_value());

    CHECK(scenario.commodities.size() == 1);
    CHECK(scenario.verify().has_value());
}

TEST_CASE("TrafficGrowth preserves valid state", "[ScenarioAnalysis][TrafficGrowth]") {
    auto graph = createScenarioTestGraph();
    TrafficGrowth scenario(graph);

    REQUIRE(scenario.setGlobalScale(1.5).has_value());
    REQUIRE(scenario.growPerCommodity(0, 1, 2.0).has_value());

    CHECK_FALSE(scenario.setGlobalScale(-2.0).has_value());
    CHECK_FALSE(scenario.growPerCommodity(0, 1, 3.0).has_value());

    CHECK(scenario.globalScale == Catch::Approx(1.5));
    CHECK(scenario.commodities.size() == 1);
    CHECK(scenario.verify().has_value());
}
