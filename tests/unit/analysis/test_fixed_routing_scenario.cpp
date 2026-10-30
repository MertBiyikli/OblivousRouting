
//
// Fixed-routing scenario transformation tests.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>

#include "analysis/fixed_routing_scenario.h"

static NetworkData makeBaseline() {
    NetworkData result{};

    result.links = {
        {
            .edge_id = 0,
            .capacity = 200.0,
            .load = 100.0,
            .utilization = 0.5,
            .source = 0,
            .target = 1
        },
        {
            .edge_id = 2,
            .capacity = 100.0,
            .load = 150.0,
            .utilization = 1.5,
            .source = 1,
            .target = 2
        }
    };

    result.max_congestion = 1.5;
    result.average_congestion = 1.0;
    result.total_demand_volume = 250.0;

    return result;
}

static std::size_t countOverloaded(const NetworkData& network) {
    return static_cast<std::size_t>(
        std::count_if(
            network.links.begin(),
            network.links.end(),
            [](const LinkData& link) {
                return link.utilization > 1.0;
            }
        )
    );
}

TEST_CASE("Identity preserves metrics", "[fixed-routing-scenario]") {
    const auto baseline = makeBaseline();

    auto transformed = transformScenario(baseline, 1.0, {});

    REQUIRE(transformed.has_value());

    CHECK(transformed->max_congestion == Catch::Approx(1.5));
    CHECK(transformed->average_congestion == Catch::Approx(1.0));

    CHECK(transformed->links[0].utilization == Catch::Approx(0.5));
    CHECK(transformed->links[1].utilization == Catch::Approx(1.5));

    CHECK(countOverloaded(*transformed) == 1);
}

TEST_CASE("Demand scaling", "[fixed-routing-scenario]") {
    const auto baseline = makeBaseline();

    auto transformed = transformScenario(baseline, 2.0, {});

    REQUIRE(transformed.has_value());

    CHECK(transformed->links[0].load == Catch::Approx(200.0));
    CHECK(transformed->links[1].load == Catch::Approx(300.0));

    CHECK(transformed->links[0].utilization == Catch::Approx(1.0));
    CHECK(transformed->links[1].utilization == Catch::Approx(3.0));

    CHECK(transformed->max_congestion == Catch::Approx(3.0));
    CHECK(transformed->average_congestion == Catch::Approx(2.0));
    CHECK(transformed->total_demand_volume == Catch::Approx(500.0));

    CHECK(countOverloaded(*transformed) == 1);
}

TEST_CASE("Capacity upgrade preserves loads", "[fixed-routing-scenario]") {
    const auto baseline = makeBaseline();

    auto transformed = transformScenario(
        baseline,
        1.0,
        {{2, 2.0}}
    );

    REQUIRE(transformed.has_value());

    CHECK(transformed->links[1].load == Catch::Approx(150.0));
    CHECK(transformed->links[1].capacity == Catch::Approx(200.0));
    CHECK(transformed->links[1].utilization == Catch::Approx(0.75));

    CHECK(transformed->links[0].utilization == Catch::Approx(0.5));
    CHECK(transformed->max_congestion == Catch::Approx(0.75));

    CHECK(countOverloaded(*transformed) == 0);
}

TEST_CASE("Combined transformation", "[fixed-routing-scenario]") {
    const auto baseline = makeBaseline();

    auto transformed = transformScenario(
        baseline,
        1.5,
        {{2, 2.0}}
    );

    REQUIRE(transformed.has_value());

    CHECK(transformed->links[0].load == Catch::Approx(150.0));
    CHECK(transformed->links[0].utilization == Catch::Approx(0.75));

    CHECK(transformed->links[1].load == Catch::Approx(225.0));
    CHECK(transformed->links[1].capacity == Catch::Approx(200.0));
    CHECK(transformed->links[1].utilization == Catch::Approx(1.125));

    CHECK(transformed->max_congestion == Catch::Approx(1.125));

    CHECK(countOverloaded(*transformed) == 1);
}

TEST_CASE("Baseline remains unchanged", "[fixed-routing-scenario]") {
    const auto baseline = makeBaseline();

    auto transformed = transformScenario(
        baseline,
        2.0,
        {{2, 2.0}}
    );

    REQUIRE(transformed.has_value());

    CHECK(baseline.links[1].load == Catch::Approx(150.0));
    CHECK(baseline.links[1].capacity == Catch::Approx(100.0));
    CHECK(baseline.links[1].utilization == Catch::Approx(1.5));
    CHECK(baseline.max_congestion == Catch::Approx(1.5));
}

TEST_CASE("Rejects invalid inputs", "[fixed-routing-scenario]") {
    const auto baseline = makeBaseline();

    CHECK_FALSE(
        transformScenario(baseline, 0.0, {}).has_value()
    );

    CHECK_FALSE(
        transformScenario(baseline, -1.0, {}).has_value()
    );

    CHECK_FALSE(
        transformScenario(
            baseline,
            std::numeric_limits<double>::quiet_NaN(),
            {}
        ).has_value()
    );

    CHECK_FALSE(
        transformScenario(baseline, 1.0, {{2, 0.5}}).has_value()
    );

    CHECK_FALSE(
        transformScenario(baseline, 1.0, {{999, 2.0}}).has_value()
    );
}

TEST_CASE("Rejects overflow", "[fixed-routing-scenario]") {
    auto baseline = makeBaseline();

    baseline.links[0].load =
        std::numeric_limits<double>::max();

    auto transformed = transformScenario(baseline, 2.0, {});

    CHECK_FALSE(transformed.has_value());
}
