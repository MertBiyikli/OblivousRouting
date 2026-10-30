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
#include "analysis/workspace/workspace_configuration.h"
#include "analysis/workspace/workspace_executer.h"
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

TEST_CASE("Prepared analysis accepts fixed external demands", "[workspace]") {
    Graph graph = makeN1TestGraph();

    Config cfg{};
    cfg.solvers = {SolverType::ELECTRICAL_SKETCHING};
    cfg.demand_models = {DemandModelType::UNIFORM};
    cfg.demand_scale = 1.0;

    const auto originalFingerprint = optimized::GraphFingerprint::compute(graph);

    auto originalDemands = ScenarioEngine::computeTotalDemandVolume(graph, cfg);

    REQUIRE(originalDemands);

    auto demandVector = DemandVector::fromDemands(originalDemands->first, graph.getNumNodes());

    REQUIRE(demandVector);

    auto scaled = demandVector->scaled(1.3);

    REQUIRE(scaled);

    auto fixedDemands = scaled->toDemands();

    REQUIRE(fixedDemands);

    auto prepared = ScenarioEngine::prepare(graph, cfg, ExecutionMode::FAST, fixedDemands.value());

    REQUIRE(prepared);

    REQUIRE(prepared->total_demand_volume == Catch::Approx(originalDemands->second * 1.3));

    REQUIRE(prepared->graph_fingerprint == originalFingerprint);

    REQUIRE(optimized::GraphFingerprint::compute(graph) == originalFingerprint);

    REQUIRE(prepared->solvers.size() == 1);

    REQUIRE(prepared->solvers.front().baseline.total_demand_volume == Catch::Approx(prepared->total_demand_volume));
}

TEST_CASE("Workspace configuration reset removes modifications", "[workspace]") {
    WorkspaceConfiguration workspace;

    workspace.capacity_upgrades.push_back({4, 1.5});
    workspace.demand_scale = 1.3;

    workspace.commodity_overrides.push_back({0, 2, 2.0});

    REQUIRE(workspace.isModified());

    workspace.reset();

    REQUIRE_FALSE(workspace.isModified());
    REQUIRE(workspace.capacity_upgrades.empty());
    REQUIRE(workspace.commodity_overrides.empty());
    REQUIRE(workspace.demand_scale == 1.0);
}

#include "analysis/workspace/workspace_preparer.h"

TEST_CASE("N-1 evaluates capacity and traffic modified workspace", "[workspace][n1]") {
    Graph graph = makeN1TestGraph();

    Config cfg{};
    cfg.solvers = {SolverType::ELECTRICAL_SKETCHING};
    cfg.demand_models = {DemandModelType::UNIFORM};
    cfg.demand_scale = 1.0;

    WorkspaceConfiguration modifications;
    modifications.capacity_upgrades = {{4, 1.5}};
    modifications.demand_scale = 1.3;

    auto workspace = WorkspacePreparer::prepare(graph, cfg, modifications);

    REQUIRE(workspace);

    const auto effectiveFingerprint = workspace->effective_fingerprint;

    REQUIRE(graph.edgeData(4).capacity == Catch::Approx(15.0));

    auto prepared = ScenarioEngine::prepare(graph, cfg, ExecutionMode::FAST, workspace->effective_demands);

    REQUIRE(prepared);

    auto result = analyzeN1(graph, prepared.value());

    REQUIRE(result);

    CHECK(result->graph_fingerprint == effectiveFingerprint);

    CHECK(result->baseline.total_demand_volume == Catch::Approx(workspace->effective_demand_volume));

    CHECK(workspace->effective_demand_volume == Catch::Approx(workspace->original_demand_volume * 1.3));

    CHECK(graph.edgeData(4).capacity == Catch::Approx(15.0));

    CHECK(optimized::GraphFingerprint::compute(graph) == effectiveFingerprint);

    for (const auto& outcome : result->outcomes) {
        CHECK(graph.edgeAlive(outcome.edge_id));

        if (outcome.status == FailureStatus::SUCCESS) {
            REQUIRE(outcome.reoptimized.has_value());

            for (const auto& link : outcome.reoptimized->links) {
                if (link.edge_id == 4) {
                    CHECK(link.capacity == Catch::Approx(15.0));
                }
            }
        }
    }
}

TEST_CASE("Single link failure uses modified workspace", "[workspace][composition]") {
    Graph graph = makeN1TestGraph();

    Config cfg{};
    cfg.solvers = {SolverType::ELECTRICAL_SKETCHING};
    cfg.demand_models = {DemandModelType::UNIFORM};
    cfg.demand_scale = 1.0;

    WorkspaceConfiguration modifications;
    modifications.capacity_upgrades = {{4, 1.5}};
    modifications.demand_scale = 1.3;

    auto workspace = WorkspacePreparer::prepare(graph, cfg, modifications);
    REQUIRE(workspace);

    const auto fingerprint = workspace->effective_fingerprint;

    auto prepared = ScenarioEngine::prepare(graph, cfg, ExecutionMode::FAST, workspace->effective_demands);
    REQUIRE(prepared);

    auto result = WorkspaceExecutor::linkFailure(graph, prepared.value(), 14);

    REQUIRE(result);

    CHECK(optimized::GraphFingerprint::compute(graph) == fingerprint);

    CHECK(graph.edgeData(4).capacity == Catch::Approx(15.0));

    CHECK(graph.edgeAlive(14));

    REQUIRE(result->solver_results.size() == 1);
    CHECK(result->solver_results.front().modified_Data.has_value());
}

TEST_CASE("Capacity planning uses effective workspace capacities", "[workspace][composition]") {
    Graph graph = makeN1TestGraph();

    Config cfg{};
    cfg.solvers = {SolverType::ELECTRICAL_SKETCHING};
    cfg.demand_models = {DemandModelType::UNIFORM};
    cfg.demand_scale = 1.0;

    WorkspaceConfiguration modifications;
    modifications.capacity_upgrades = {{4, 1.5}};
    modifications.demand_scale = 1.3;

    auto workspace = WorkspacePreparer::prepare(graph, cfg, modifications);
    REQUIRE(workspace);

    auto prepared = ScenarioEngine::prepare(graph, cfg, ExecutionMode::FAST, workspace->effective_demands);
    REQUIRE(prepared);

    const auto fingerprint = optimized::GraphFingerprint::compute(graph);

    auto result = WorkspaceExecutor::capacityUpgrade(graph, prepared.value(), 4, 2.0);

    REQUIRE(result);

    CHECK(graph.edgeData(4).capacity == Catch::Approx(15.0));

    CHECK(optimized::GraphFingerprint::compute(graph) == fingerprint);

    REQUIRE(result->solver_results.size() == 1);

    const auto& solver = result->solver_results.front();

    REQUIRE(solver.modified_Data.has_value());

    bool found_edge = false;

    for (const auto& link : solver.modified_Data->links) {
        if (link.edge_id == 4) {
            found_edge = true;
            CHECK(link.capacity == Catch::Approx(30.0));
        }
    }

    CHECK(found_edge);
}

TEST_CASE("Different analyses preserve the shared workspace", "[workspace][composition]") {
    Graph graph = makeN1TestGraph();

    Config cfg{};
    cfg.solvers = {SolverType::ELECTRICAL_SKETCHING};
    cfg.demand_models = {DemandModelType::UNIFORM};
    cfg.demand_scale = 1.0;

    WorkspaceConfiguration modifications;
    modifications.capacity_upgrades = {{4, 1.5}};
    modifications.demand_scale = 1.3;

    auto workspace = WorkspacePreparer::prepare(graph, cfg, modifications);
    REQUIRE(workspace);

    auto prepared = ScenarioEngine::prepare(graph, cfg, ExecutionMode::FAST, workspace->effective_demands);
    REQUIRE(prepared);

    const auto expected_fingerprint = workspace->effective_fingerprint;

    SECTION("Link failure followed by capacity planning") {
        auto failure = WorkspaceExecutor::linkFailure(graph, prepared.value(), 14);
        REQUIRE(failure);

        CHECK(optimized::GraphFingerprint::compute(graph) == expected_fingerprint);

        auto capacity = WorkspaceExecutor::capacityUpgrade(graph, prepared.value(), 4, 2.0);
        REQUIRE(capacity);

        CHECK(optimized::GraphFingerprint::compute(graph) == expected_fingerprint);
    }

    SECTION("Capacity planning followed by N-1") {
        auto capacity = WorkspaceExecutor::capacityUpgrade(graph, prepared.value(), 4, 2.0);
        REQUIRE(capacity);

        auto resilience = analyzeN1(graph, prepared.value());
        REQUIRE(resilience);

        CHECK(optimized::GraphFingerprint::compute(graph) == expected_fingerprint);
    }

    CHECK(graph.edgeData(4).capacity == Catch::Approx(15.0));
}
