//
// Created by Mert on 30.10.26.
//

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "analysis/workspace/workspace_preparer.h"

#include "analysis/traffic_change/demand_vector.h"
#include <cmath>
#include <vector>

namespace {

using Graph = optimized::Graph<EdgeData>;

Graph makeWorkspaceGraph() {
    const std::vector<Graph::InputEdge> edges{{0, 1, {10.0, 1.0}}, {1, 2, {20.0, 1.0}}, {0, 2, {15.0, 1.0}}};

    return Graph(3, edges);
}

Config makeWorkspaceConfig() {
    Config cfg{};
    cfg.solvers = {SolverType::ELECTRICAL_SKETCHING};
    cfg.demand_models = {DemandModelType::UNIFORM};
    cfg.demand_scale = 1.0;
    return cfg;
}

} // namespace

TEST_CASE("Workspace composes traffic growth and capacity upgrades", "[workspace][composition]") {
    auto graph = makeWorkspaceGraph();
    auto cfg = makeWorkspaceConfig();

    const auto originalFingerprint = optimized::GraphFingerprint::compute(graph);

    WorkspaceConfiguration modifications;
    modifications.demand_scale = 1.3;
    modifications.capacity_upgrades = {{0, 1.5}};

    auto prepared = WorkspacePreparer::prepare(graph, cfg, modifications);

    REQUIRE(prepared);

    CHECK(prepared->original_fingerprint == originalFingerprint);

    CHECK(graph.edgeData(0).capacity == Catch::Approx(15.0));

    CHECK(prepared->effective_demand_volume == Catch::Approx(prepared->original_demand_volume * 1.3));

    CHECK(prepared->effective_fingerprint == optimized::GraphFingerprint::compute(graph));

    auto analysis = ScenarioEngine::prepare(graph, cfg, ExecutionMode::FAST, prepared->effective_demands);

    REQUIRE(analysis);

    CHECK(analysis->total_demand_volume == Catch::Approx(prepared->effective_demand_volume));
}

TEST_CASE("Invalid workspace modifications preserve the graph", "[workspace][composition]") {
    auto graph = makeWorkspaceGraph();
    auto cfg = makeWorkspaceConfig();

    const auto originalFingerprint = optimized::GraphFingerprint::compute(graph);

    WorkspaceConfiguration modifications;
    modifications.capacity_upgrades = {{0, 1.5}, {999, 2.0}};

    auto prepared = WorkspacePreparer::prepare(graph, cfg, modifications);

    REQUIRE_FALSE(prepared);

    CHECK(graph.edgeData(0).capacity == Catch::Approx(10.0));

    CHECK(optimized::GraphFingerprint::compute(graph) == originalFingerprint);
}

TEST_CASE("Workspace reset reconstructs the original state", "[workspace][composition]") {
    auto graph = makeWorkspaceGraph();
    auto cfg = makeWorkspaceConfig();

    WorkspaceConfiguration modifications;
    modifications.capacity_upgrades = {{0, 1.5}};
    modifications.demand_scale = 1.3;

    auto modified = WorkspacePreparer::prepare(graph, cfg, modifications);

    REQUIRE(modified);
    CHECK(graph.edgeData(0).capacity == Catch::Approx(15.0));

    // Recreate from the registered original graph.
    // Reset does not mutate the existing modified graph.
    modifications.reset();
    graph = makeWorkspaceGraph();

    auto restored = WorkspacePreparer::prepare(graph, cfg, modifications);

    REQUIRE(restored);

    CHECK(graph.edgeData(0).capacity == Catch::Approx(10.0));

    CHECK(restored->original_fingerprint == restored->effective_fingerprint);

    CHECK(restored->original_demand_volume == Catch::Approx(restored->effective_demand_volume));
}

TEST_CASE("Workspace composes capacity global and commodity modifications", "[workspace][composition][cross-functional]") {
    auto graph = makeWorkspaceGraph();
    auto cfg = makeWorkspaceConfig();

    const auto original_fingerprint = optimized::GraphFingerprint::compute(graph);

    const int checkpoint = graph.checkpoint();

    WorkspaceConfiguration modifications;
    modifications.capacity_upgrades = {{0, 1.5}};
    modifications.demand_scale = 1.2;
    modifications.commodity_overrides = {{0, 2, 2.0}};

    auto workspace = WorkspacePreparer::prepare(graph, cfg, modifications);

    REQUIRE(workspace);

    // Capacity modifications affect the graph.
    CHECK(graph.edgeData(0).capacity == Catch::Approx(15.0));

    CHECK(workspace->original_fingerprint == original_fingerprint);

    CHECK(workspace->effective_fingerprint == optimized::GraphFingerprint::compute(graph));

    auto original = DemandVector::fromDemands(workspace->original_demands, graph.getNumNodes());

    auto effective = DemandVector::fromDemands(workspace->effective_demands, graph.getNumNodes());

    REQUIRE(original);
    REQUIRE(effective);

    // Every original commodity must have the expected
    // multiplicative transformation.
    REQUIRE(original->entries().size() == effective->entries().size());

    bool found_override = false;

    for (const auto& entry : original->entries()) {
        double expected_factor = 1.2;

        if (entry.source == 0 && entry.target == 2) {
            expected_factor *= 2.0;
            found_override = true;
        }

        bool found = false;

        for (const auto& modified : effective->entries()) {
            if (entry.source == modified.source && entry.target == modified.target) {

                CHECK(modified.volume == Catch::Approx(entry.volume * expected_factor));

                found = true;
                break;
            }
        }

        REQUIRE(found);
    }

    REQUIRE(found_override);

    CHECK(effective->totalVolume() == Catch::Approx(workspace->effective_demand_volume));

    // Restoring the original graph must be possible
    // after successful workspace preparation.
    auto restored = graph.rollback(checkpoint);

    REQUIRE(restored);

    CHECK(optimized::GraphFingerprint::compute(graph) == original_fingerprint);

    CHECK(graph.edgeData(0).capacity == Catch::Approx(10.0));
}