//
// Created by Mert on 22.10.26.
//

#include <catch2/catch_test_macros.hpp>

#include <catch2/catch_approx.hpp>
#include <cmath>
#include <limits>
#include <unordered_map>

#include "analysis/fixed_routing_scenario.h"

static RoutingVisualizationResult makeBaseline() {
    RoutingVisualizationResult result{};

    result.graph_name = "TEST_CASE-network";
    result.solver_name = "electrical";
    result.demand_model = "uniform";

    result.edges = {
        {.id = 0, .source = 0, .target = 1,
         .capacity = 200.0, .load = 100.0,
         .utilization = 0.5, .enabled = true},

        {.id = 2, .source = 1, .target = 2,
         .capacity = 100.0, .load = 150.0,
         .utilization = 1.5, .enabled = true}
    };

    result.summary.maximum_congestion = 1.5;
    result.summary.average_utilization = 1.0;
    result.summary.overloaded_edges = 1;
    result.summary.total_routed_demand = 250.0;

    return result;
}

TEST_CASE("Identity preserves metrics", "[fixed-routing-scenario]") {
    const auto baseline = makeBaseline();

    auto transformed = transformScenario(baseline, 1.0, {});

    REQUIRE(transformed.has_value());

    CHECK(transformed->summary.maximum_congestion == Catch::Approx(1.5));
    CHECK(transformed->summary.average_utilization == Catch::Approx(1.0));
    CHECK(transformed->edges[0].utilization == Catch::Approx(0.5));
    CHECK(transformed->edges[1].utilization == Catch::Approx(1.5));
    CHECK(transformed->summary.overloaded_edges == 1);
}

TEST_CASE("Demand scaling", "[fixed-routing-scenario]") {
    const auto baseline = makeBaseline();

    auto transformed = transformScenario(baseline, 2.0, {});

    REQUIRE(transformed.has_value());

    CHECK(transformed->edges[0].load == Catch::Approx(200.0));
    CHECK(transformed->edges[1].load == Catch::Approx(300.0));

    CHECK(transformed->edges[0].utilization == Catch::Approx(1.0));
    CHECK(transformed->edges[1].utilization == Catch::Approx(3.0));

    CHECK(transformed->summary.maximum_congestion == Catch::Approx(3.0));
    CHECK(transformed->summary.average_utilization == Catch::Approx(2.0));
    CHECK(transformed->summary.total_routed_demand == Catch::Approx(500.0));
    CHECK(transformed->summary.overloaded_edges == 1);
}

TEST_CASE("Capacity upgrade preserves loads", "[fixed-routing-scenario]") {
    const auto baseline = makeBaseline();

    auto transformed = transformScenario(baseline, 1.0, {{2, 2.0}});

    REQUIRE(transformed.has_value());

    CHECK(transformed->edges[1].load == Catch::Approx(150.0));
    CHECK(transformed->edges[1].capacity == Catch::Approx(200.0));
    CHECK(transformed->edges[1].utilization == Catch::Approx(0.75));

    CHECK(transformed->edges[0].utilization == Catch::Approx(0.5));
    CHECK(transformed->summary.maximum_congestion == Catch::Approx(0.75));
    CHECK(transformed->summary.overloaded_edges == 0);
}

TEST_CASE("Combined transformation", "[fixed-routing-scenario]") {
    const auto baseline = makeBaseline();

    auto transformed = transformScenario(baseline, 1.5, {{2, 2.0}});

    REQUIRE(transformed.has_value());

    CHECK(transformed->edges[0].load == Catch::Approx(150.0));
    CHECK(transformed->edges[0].utilization == Catch::Approx(0.75));

    CHECK(transformed->edges[1].load == Catch::Approx(225.0));
    CHECK(transformed->edges[1].capacity == Catch::Approx(200.0));
    CHECK(transformed->edges[1].utilization == Catch::Approx(1.125));

    CHECK(transformed->summary.maximum_congestion == Catch::Approx(1.125));
    CHECK(transformed->summary.overloaded_edges == 1);
}

TEST_CASE("Baseline remains unchanged", "[fixed-routing-scenario]") {
    const auto baseline = makeBaseline();

    auto transformed = transformScenario(baseline, 2.0, {{2, 2.0}});

    REQUIRE(transformed.has_value());

    CHECK(baseline.edges[1].load == Catch::Approx(150.0));
    CHECK(baseline.edges[1].capacity == Catch::Approx(100.0));
    CHECK(baseline.edges[1].utilization == Catch::Approx(1.5));
    CHECK(baseline.summary.maximum_congestion == Catch::Approx(1.5));
}

TEST_CASE("Rejects invalid inputs", "[fixed-routing-scenario]") {
    const auto baseline = makeBaseline();

    CHECK_FALSE(transformScenario(baseline, 0.0, {}).has_value());
    CHECK_FALSE(transformScenario(baseline, -1.0, {}).has_value());

    CHECK_FALSE(transformScenario(
        baseline, std::numeric_limits<double>::quiet_NaN(), {}
    ).has_value());

    CHECK_FALSE(transformScenario(
        baseline, 1.0, {{2, 0.5}}
    ).has_value());

    CHECK_FALSE(transformScenario(
        baseline, 1.0, {{999, 2.0}}
    ).has_value());
}

TEST_CASE("Rejects overflow", "[fixed-routing-scenario]") {
    auto baseline = makeBaseline();

    baseline.edges[0].load = std::numeric_limits<double>::max();

    auto transformed = transformScenario(baseline, 2.0, {});

    CHECK_FALSE(transformed.has_value());
}