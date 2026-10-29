
//
// Scenario Analysis Tests
// Updated for checkpoint-based scenario execution.
//

#include "analysis/capacity_upgrade/capacity_upgrade.h"
#include "analysis/failures/link_failures.h"
#include "analysis/graph_rollback_guard.h"
#include "analysis/scenario_analysis.h"
#include "analysis/scenario_engine.h"
#include "utils.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <vector>



// ============================================================
// Capacity Upgrade
// ============================================================

TEST_CASE("CapacityUpgrade accepts valid upgrades", "[ScenarioAnalysis][CapacityUpgrade]") {

    auto graph = createScenarioTestGraph();
    CapacityUpgrade scenario;

    REQUIRE(scenario.upgrade(0, 1.5).has_value());
    REQUIRE(scenario.upgrade(2, 2.0).has_value());
    REQUIRE(scenario.verify(graph).has_value());

    CHECK(scenario.kind() == ScenarioKind::CAPACITY_UPGRADE);
}

TEST_CASE("CapacityUpgrade rejects invalid edge IDs", "[ScenarioAnalysis][CapacityUpgrade]") {

    auto graph = createScenarioTestGraph();

    SECTION("Negative edge ID") {
        CapacityUpgrade scenario;

        CHECK_FALSE(scenario.upgrade(-1, 1.5).has_value());
        REQUIRE(scenario.verify(graph).has_value());
    }

    SECTION("Out of range edge ID") {
        CapacityUpgrade scenario;

        // Registration does not know the graph.
        REQUIRE(scenario.upgrade(graph.globalEdgeCount() * 2, 1.5).has_value());

        CHECK_FALSE(scenario.verify(graph).has_value());
    }
}

TEST_CASE("CapacityUpgrade rejects invalid factors", "[ScenarioAnalysis][CapacityUpgrade]") {

    const double nan = std::numeric_limits<double>::quiet_NaN();

    const double inf = std::numeric_limits<double>::infinity();

    CapacityUpgrade scenario;

    CHECK_FALSE(scenario.upgrade(0, -1.0).has_value());
    CHECK_FALSE(scenario.upgrade(0, 0.0).has_value());
    CHECK_FALSE(scenario.upgrade(0, 0.5).has_value());
    CHECK_FALSE(scenario.upgrade(0, nan).has_value());
    CHECK_FALSE(scenario.upgrade(0, inf).has_value());

    // A rejected request must not reserve edge 0.
    REQUIRE(scenario.upgrade(0, 1.5).has_value());
}

TEST_CASE("CapacityUpgrade rejects duplicate physical links", "[ScenarioAnalysis][CapacityUpgrade]") {

    auto graph = createScenarioTestGraph();

    CapacityUpgrade scenario;

    const int edge = 0;
    const int reverse = graph.reverse(edge).id;

    REQUIRE(scenario.upgrade(edge, 1.5).has_value());

    CHECK_FALSE(scenario.upgrade(edge, 2.0).has_value());

    // The reverse arc can be registered because the scenario
    // does not yet know the graph topology.
    REQUIRE(scenario.upgrade(reverse, 2.0).has_value());

    // Both IDs describe the same physical link.
    CHECK_FALSE(scenario.verify(graph).has_value());
}

TEST_CASE("CapacityUpgrade does not modify graph during configuration", "[ScenarioAnalysis][CapacityUpgrade]") {

    auto graph = createScenarioTestGraph();

    const auto fingerprint = optimized::GraphFingerprint::compute(graph);

    const int checkpoint = graph.checkpoint();

    CapacityUpgrade scenario;

    REQUIRE(scenario.upgrade(0, 1.5).has_value());
    REQUIRE(scenario.upgrade(2, 2.0).has_value());
    REQUIRE(scenario.verify(graph).has_value());

    CHECK(graph.checkpoint() == checkpoint);

    CHECK(optimized::GraphFingerprint::compute(graph) == fingerprint);

    CHECK(graph.edgeData(0).capacity == Catch::Approx(10.0));
    CHECK(graph.edgeData(2).capacity == Catch::Approx(20.0));
}

TEST_CASE("CapacityUpgrade applies correct physical link capacities", "[ScenarioAnalysis][CapacityUpgrade]") {

    auto graph = createScenarioTestGraph();

    CapacityUpgrade scenario;

    REQUIRE(scenario.upgrade(0, 1.5).has_value());
    REQUIRE(scenario.upgrade(2, 2.0).has_value());
    REQUIRE(scenario.verify(graph).has_value());

    {
        GraphRollbackGuard guard(graph);

        REQUIRE(scenario.apply(graph).has_value());

        CHECK(graph.edgeData(0).capacity == Catch::Approx(15.0));

        CHECK(graph.edgeData(2).capacity == Catch::Approx(40.0));

        // Reverse arcs share physical edge data.
        CHECK(graph.edgeData(1).capacity == Catch::Approx(15.0));

        CHECK(graph.edgeData(3).capacity == Catch::Approx(40.0));

        // Two physical links => two mutations.
        CHECK(graph.checkpoint() == guard.checkpoint() + 2);
    }

    CHECK(graph.edgeData(0).capacity == Catch::Approx(10.0));

    CHECK(graph.edgeData(2).capacity == Catch::Approx(20.0));
}

TEST_CASE("CapacityUpgrade rejects capacity overflow", "[ScenarioAnalysis][CapacityUpgrade]") {

    auto graph = createScenarioTestGraph();

    CapacityUpgrade scenario;

    REQUIRE(scenario.upgrade(0, std::numeric_limits<double>::max()).has_value());

    const auto original = optimized::GraphFingerprint::compute(graph);

    CHECK_FALSE(scenario.verify(graph).has_value());

    CHECK(optimized::GraphFingerprint::compute(graph) == original);
}

TEST_CASE("CapacityUpgrade preserves valid state after rejected update", "[ScenarioAnalysis][CapacityUpgrade]") {

    auto graph = createScenarioTestGraph();

    CapacityUpgrade scenario;

    REQUIRE(scenario.upgrade(0, 1.5).has_value());

    CHECK_FALSE(scenario.upgrade(-1, 2.0).has_value());
    CHECK_FALSE(scenario.upgrade(0, 3.0).has_value());

    REQUIRE(scenario.verify(graph).has_value());

    {
        GraphRollbackGuard guard(graph);

        REQUIRE(scenario.apply(graph).has_value());

        CHECK(graph.edgeData(0).capacity == Catch::Approx(15.0));
    }

    CHECK(graph.edgeData(0).capacity == Catch::Approx(10.0));
}

// ============================================================
// Link Failure
// ============================================================

TEST_CASE("LinkFailure accepts a single failure", "[ScenarioAnalysis][LinkFailure]") {

    auto graph = createScenarioTestGraph();

    LinkFailure scenario;

    REQUIRE(scenario.linkFailure(0).has_value());
    REQUIRE(scenario.verify(graph).has_value());

    CHECK(scenario.kind() == ScenarioKind::LINK_FAILURE);
}

TEST_CASE("LinkFailure accepts multiple physical links", "[ScenarioAnalysis][LinkFailure]") {

    auto graph = createScenarioTestGraph();

    LinkFailure scenario;

    REQUIRE(scenario.linkFailure(0).has_value());
    REQUIRE(scenario.linkFailure(2).has_value());
    REQUIRE(scenario.linkFailure(4).has_value());

    REQUIRE(scenario.verify(graph).has_value());

    {
        GraphRollbackGuard guard(graph);

        REQUIRE(scenario.apply(graph).has_value());

        CHECK_FALSE(graph.edgeAlive(0));
        CHECK_FALSE(graph.edgeAlive(2));
        CHECK_FALSE(graph.edgeAlive(4));

        CHECK(graph.getNumUndirectedEdges() == 1);
    }

    CHECK(graph.getNumUndirectedEdges() == 4);
}

TEST_CASE("LinkFailure rejects invalid edges", "[ScenarioAnalysis][LinkFailure]") {

    auto graph = createScenarioTestGraph();

    SECTION("Negative edge ID") {
        LinkFailure scenario;

        CHECK_FALSE(scenario.linkFailure(-1).has_value());
        REQUIRE(scenario.verify(graph).has_value());
    }

    SECTION("Out of range edge ID") {
        LinkFailure scenario;

        REQUIRE(scenario.linkFailure(graph.globalEdgeCount() * 2).has_value());

        CHECK_FALSE(scenario.verify(graph).has_value());
    }
}

TEST_CASE("LinkFailure rejects reverse-edge duplicates", "[ScenarioAnalysis][LinkFailure]") {

    auto graph = createScenarioTestGraph();

    LinkFailure scenario;

    const int edge = 0;
    const int reverse = graph.reverse(edge).id;

    REQUIRE(scenario.linkFailure(edge).has_value());

    CHECK_FALSE(scenario.linkFailure(edge).has_value());

    REQUIRE(scenario.linkFailure(reverse).has_value());

    // Detect same physical link during verification.
    CHECK_FALSE(scenario.verify(graph).has_value());
}

TEST_CASE("LinkFailure rejects inactive links", "[ScenarioAnalysis][LinkFailure]") {

    auto graph = createScenarioTestGraph();

    REQUIRE(graph.removeEdge(0));

    LinkFailure scenario;

    REQUIRE(scenario.linkFailure(0).has_value());

    CHECK_FALSE(scenario.verify(graph).has_value());

    REQUIRE(graph.restoreAll().has_value());
}

TEST_CASE("LinkFailure preserves baseline during configuration", "[ScenarioAnalysis][LinkFailure]") {

    auto graph = createScenarioTestGraph();

    const auto fingerprint = optimized::GraphFingerprint::compute(graph);

    const int checkpoint = graph.checkpoint();

    LinkFailure scenario;

    REQUIRE(scenario.linkFailure(0).has_value());
    REQUIRE(scenario.linkFailure(2).has_value());
    REQUIRE(scenario.verify(graph).has_value());

    CHECK(graph.checkpoint() == checkpoint);

    CHECK(graph.edgeAlive(0));
    CHECK(graph.edgeAlive(2));

    CHECK(optimized::GraphFingerprint::compute(graph) == fingerprint);
}

// ============================================================
// Graph Checkpoint / Rollback
// ============================================================

TEST_CASE("CapacityUpgrade rollback restores original graph", "[ScenarioAnalysis][Rollback]") {

    auto graph = createScenarioTestGraph();

    const auto fingerprint = optimized::GraphFingerprint::compute(graph);

    const int checkpoint = graph.checkpoint();

    CapacityUpgrade scenario;

    REQUIRE(scenario.upgrade(0, 1.5).has_value());
    REQUIRE(scenario.verify(graph).has_value());

    {
        GraphRollbackGuard guard(graph);

        REQUIRE(scenario.apply(graph).has_value());

        CHECK(graph.edgeData(0).capacity == Catch::Approx(15.0));
    }

    CHECK(graph.checkpoint() == checkpoint);

    CHECK(graph.edgeData(0).capacity == Catch::Approx(10.0));

    CHECK(optimized::GraphFingerprint::compute(graph) == fingerprint);
}

TEST_CASE("LinkFailure rollback restores original graph", "[ScenarioAnalysis][Rollback]") {

    auto graph = createScenarioTestGraph();

    const auto fingerprint = optimized::GraphFingerprint::compute(graph);

    const int checkpoint = graph.checkpoint();

    LinkFailure scenario;

    REQUIRE(scenario.linkFailure(0).has_value());
    REQUIRE(scenario.verify(graph).has_value());

    {
        GraphRollbackGuard guard(graph);

        REQUIRE(scenario.apply(graph).has_value());

        CHECK_FALSE(graph.edgeAlive(0));
        CHECK_FALSE(graph.edgeAlive(1));

        CHECK(graph.getNumUndirectedEdges() == 3);
    }

    CHECK(graph.edgeAlive(0));
    CHECK(graph.edgeAlive(1));

    CHECK(graph.getNumUndirectedEdges() == 4);
    CHECK(graph.checkpoint() == checkpoint);

    CHECK(optimized::GraphFingerprint::compute(graph) == fingerprint);
}

TEST_CASE("Nested scenarios restore correct checkpoints", "[ScenarioAnalysis][Rollback]") {

    auto graph = createScenarioTestGraph();

    CapacityUpgrade upgrade;
    LinkFailure failure;

    REQUIRE(upgrade.upgrade(0, 1.5).has_value());
    REQUIRE(failure.linkFailure(2).has_value());

    REQUIRE(upgrade.verify(graph).has_value());
    REQUIRE(failure.verify(graph).has_value());

    const int baseline_checkpoint = graph.checkpoint();

    {
        GraphRollbackGuard outer(graph);

        REQUIRE(upgrade.apply(graph).has_value());

        CHECK(graph.edgeData(0).capacity == Catch::Approx(15.0));

        {
            GraphRollbackGuard inner(graph);

            REQUIRE(failure.apply(graph).has_value());

            CHECK_FALSE(graph.edgeAlive(2));

            CHECK(graph.edgeData(0).capacity == Catch::Approx(15.0));
        }

        // Inner rollback restores the upgraded graph.
        CHECK(graph.edgeAlive(2));

        CHECK(graph.edgeData(0).capacity == Catch::Approx(15.0));
    }

    // Outer rollback restores the original graph.
    CHECK(graph.edgeAlive(2));

    CHECK(graph.edgeData(0).capacity == Catch::Approx(10.0));

    CHECK(graph.checkpoint() == baseline_checkpoint);
}

TEST_CASE("Repeated scenario execution does not accumulate mutations", "[ScenarioAnalysis][Rollback]") {

    auto graph = createScenarioTestGraph();

    const auto fingerprint = optimized::GraphFingerprint::compute(graph);

    const int original_checkpoint = graph.checkpoint();

    LinkFailure failure;

    REQUIRE(failure.linkFailure(0).has_value());
    REQUIRE(failure.verify(graph).has_value());

    for (int iteration = 0; iteration < 100; ++iteration) {
        {
            GraphRollbackGuard guard(graph);

            REQUIRE(failure.apply(graph).has_value());

            CHECK_FALSE(graph.edgeAlive(0));
            CHECK(graph.getNumUndirectedEdges() == 3);
        }

        REQUIRE(graph.edgeAlive(0));

        REQUIRE(graph.checkpoint() == original_checkpoint);

        REQUIRE(optimized::GraphFingerprint::compute(graph) == fingerprint);
    }
}

TEST_CASE("Stable edge IDs survive temporary link removals", "[ScenarioAnalysis][Graph]") {

    auto graph = createScenarioTestGraph();

    const int last_edge = 6;

    const auto endpoints = graph.getEdgeEndpoints(last_edge);

    const auto capacity = graph.edgeData(last_edge).capacity;

    {
        GraphRollbackGuard guard(graph);

        REQUIRE(graph.removeEdge(0));

        // Active edge count decreases, but stable IDs remain.
        REQUIRE(graph.edgeAlive(last_edge));

        CHECK(graph.getEdgeEndpoints(last_edge) == endpoints);

        CHECK(graph.edgeData(last_edge).capacity == Catch::Approx(capacity));
    }

    CHECK(graph.edgeAlive(0));
    CHECK(graph.edgeAlive(last_edge));

    CHECK(graph.getNumUndirectedEdges() == 4);
}

TEST_CASE("LinkFailure preserves remaining stable edge IDs", "[ScenarioAnalysis][LinkFailure]") {
    auto graph = createScenarioTestGraph();

    LinkFailure scenario;

    REQUIRE(scenario.linkFailure(0).has_value());
    REQUIRE(scenario.verify(graph).has_value());

    const auto original = optimized::GraphFingerprint::compute(graph);

    const int checkpoint = graph.checkpoint();

    {
        GraphRollbackGuard guard(graph);

        REQUIRE(scenario.apply(graph).has_value());

        CHECK_FALSE(graph.edgeAlive(0));
        CHECK_FALSE(graph.edgeAlive(1));

        CHECK(graph.edgeAlive(2));
        CHECK(graph.edgeAlive(4));
        CHECK(graph.edgeAlive(6));

        CHECK(graph.getNumUndirectedEdges() == 3);

        const auto modified = optimized::GraphFingerprint::compute(graph);

        CHECK(modified != original);
    }

    CHECK(graph.checkpoint() == checkpoint);
    CHECK(graph.getNumUndirectedEdges() == 4);

    CHECK(optimized::GraphFingerprint::compute(graph) == original);
}

inline void verifyScenarioRollback(
    optimized::Graph<EdgeData>& graph,
    Scenario& scenario,
    const Config& cfg,
    ExecutionMode mode
) {
    const auto fingerprint_before =
        optimized::GraphFingerprint::compute(graph);

    const int active_edges_before =
        graph.getNumDirectedEdges();

    const auto checkpoint_before =
        graph.checkpoint();

    const int stable_edge_bound =
        graph.globalEdgeCount() * 2;

    std::vector<bool> active_before(stable_edge_bound);
    std::vector<double> capacity_before(stable_edge_bound);

    for (int e = 0; e < stable_edge_bound; ++e) {
        active_before[e] = graph.edgeAlive(e);

        if (active_before[e]) {
            capacity_before[e] = graph.edgeData(e).capacity;
        }
    }

    // First scenario execution.
    auto first = ScenarioEngine::runScenario(
        graph, scenario, cfg, mode
    );

    REQUIRE(first.has_value());

    // Graph must be restored after runScenario().
    CHECK(
        optimized::GraphFingerprint::compute(graph)
        == fingerprint_before
    );

    CHECK(graph.getNumDirectedEdges() == active_edges_before);
    CHECK(graph.checkpoint() == checkpoint_before);

    for (int e = 0; e < stable_edge_bound; ++e) {
        CHECK(graph.edgeAlive(e) == active_before[e]);

        if (active_before[e]) {
            CHECK(
                graph.edgeData(e).capacity ==
                Catch::Approx(capacity_before[e])
            );
        }
    }

    // Repeat on the same graph instance.
    auto second = ScenarioEngine::runScenario(
        graph, scenario, cfg, mode
    );

    REQUIRE(second.has_value());

    CHECK(
        second->original_fingerprint
        == first->original_fingerprint
    );

    CHECK(
        second->modified_fingerprint
        == first->modified_fingerprint
    );

    REQUIRE(
        first->solver_results.size()
        == second->solver_results.size()
    );

    for (std::size_t i = 0;
         i < first->solver_results.size(); ++i) {

        const auto& a = first->solver_results[i];
        const auto& b = second->solver_results[i];

        REQUIRE(a.baseline_Data.has_value());
        REQUIRE(b.baseline_Data.has_value());

        REQUIRE(a.modified_Data.has_value());
        REQUIRE(b.modified_Data.has_value());

        CHECK(
            a.baseline_Data->max_congestion ==
            Catch::Approx(
                b.baseline_Data->max_congestion
            )
        );

        CHECK(
            a.modified_Data->max_congestion ==
            Catch::Approx(
                b.modified_Data->max_congestion
            )
        );
    }

    // Rollback must also hold after the second execution.
    CHECK(
        optimized::GraphFingerprint::compute(graph)
        == fingerprint_before
    );

    CHECK(graph.checkpoint() == checkpoint_before);
}


TEST_CASE("Scenario rejects invalid edge IDs",
          "[ScenarioAnalysis][Validation]") {
    auto graph = createScenarioTestGraph();
    const auto before = optimized::GraphFingerprint::compute(graph);

    LinkFailure failure;
    failure.linkFailure(999999);

    auto verified = failure.verify(graph);

    REQUIRE_FALSE(verified.has_value());
    CHECK(optimized::GraphFingerprint::compute(graph) == before);
}