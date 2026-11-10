//
// Created by Mert on 29.10.26.
//

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "analysis/scenario_engine.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "analysis/scenario_engine.h"
#include "analysis/traffic_change/traffic_growth.h"
#include "analysis/workspace/workspace_preparer.h"
#include "routing/routing_engine.h"

#include <vector>

namespace {

using Graph = optimized::Graph<EdgeData>;

Graph makeTrafficIntegrationGraph() {
    const std::vector<Graph::InputEdge> edges{{0, 1, {10.0, 1.0}}, {1, 2, {20.0, 1.0}}, {0, 2, {15.0, 1.0}}};

    return Graph(3, edges);
}

Config makeTrafficIntegrationConfig() {
    Config cfg{};
    cfg.solvers = {SolverType::ELECTRICAL_SKETCHING};
    cfg.demand_models = {DemandModelType::UNIFORM};
    cfg.demand_scale = 1.0;
    return cfg;
}

} // namespace

TEST_CASE("Traffic evaluation consumes prepared workspace demands", "[traffic][workspace][cross-functional]") {
    auto graph = makeTrafficIntegrationGraph();
    auto cfg = makeTrafficIntegrationConfig();

    // 1. Compute a routing scheme on the original graph.
    RoutingEngine engine;

    auto solved = engine.solve(graph, cfg, cfg.solvers.front());

    REQUIRE(solved);
    REQUIRE(solved->scheme);
    REQUIRE(solved->scheme->isValid());

    const auto& routing = solved.value();

    // 2. Prepare the effective workspace.
    WorkspaceConfiguration modifications;

    modifications.capacity_upgrades = {{0, 1.5}};
    modifications.demand_scale = 1.2;
    modifications.commodity_overrides = {{0, 2, 2.0}};

    auto workspace = WorkspacePreparer::prepare(graph, cfg, modifications);

    REQUIRE(workspace);

    const auto fingerprint = workspace->effective_fingerprint;

    CHECK(graph.edgeData(0).capacity == Catch::Approx(15.0));

    // 3. Evaluate already-prepared original and modified demands.
    auto prepared_result = ScenarioEngine::evaluateTrafficGrowth(graph, routing, workspace->original_demands, workspace->effective_demands);

    REQUIRE(prepared_result);

    CHECK(prepared_result->graph_fingerprint == fingerprint);

    CHECK(prepared_result->baseline.total_demand_volume == Catch::Approx(workspace->original_demand_volume));

    CHECK(prepared_result->modified.total_demand_volume == Catch::Approx(workspace->effective_demand_volume));

    // 4. Evaluate the equivalent legacy transformation.
    TrafficGrowth growth;

    REQUIRE(growth.setGlobalScale(1.2));
    REQUIRE(growth.growPerCommodity(0, 2, 2.0));

    auto legacy_result = ScenarioEngine::evaluateTrafficGrowth(graph, routing, workspace->original_demands, growth);

    REQUIRE(legacy_result);

    // 5. Both paths must produce equivalent numerical results.
    CHECK(prepared_result->baseline.max_congestion == Catch::Approx(legacy_result->baseline.max_congestion).epsilon(1e-8));

    CHECK(prepared_result->modified.max_congestion == Catch::Approx(legacy_result->modified.max_congestion).epsilon(1e-8));

    CHECK(prepared_result->modified.total_demand_volume == Catch::Approx(legacy_result->modified.total_demand_volume));

    // 6. Validate individual link metrics.
    const auto& prepared_links = prepared_result->modified.links;
    const auto& legacy_links = legacy_result->modified.links;

    REQUIRE(prepared_links.size() == legacy_links.size());

    for (std::size_t i = 0; i < prepared_links.size(); ++i) {
        const auto& a = prepared_links[i];
        const auto& b = legacy_links[i];

        CHECK(a.edge_id == b.edge_id);

        CHECK(a.capacity == Catch::Approx(b.capacity));

        CHECK(a.load == Catch::Approx(b.load).epsilon(1e-8));

        CHECK(a.utilization == Catch::Approx(b.utilization).epsilon(1e-8));
    }

    // 7. Both evaluations must preserve the effective graph.
    CHECK(optimized::GraphFingerprint::compute(graph) == fingerprint);

    CHECK(graph.edgeData(0).capacity == Catch::Approx(15.0));
}