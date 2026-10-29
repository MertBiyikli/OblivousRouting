//
// Created by Mert on 29.10.26.
//

#include <catch2/catch_test_macros.hpp>

#include "analysis/traffic_change/traffic_growth.h"

#include <limits>

TEST_CASE("TrafficGrowth identity", "[traffic-growth]") {
    auto baseline = DemandVector::create(3, {
        {0, 1, 100.0},
        {1, 0, 25.0},
        {1, 2, 50.0}
    });

    REQUIRE(baseline);

    TrafficGrowth scenario;
    auto result = scenario.transform(baseline.value());

    REQUIRE(result);
    CHECK(result->totalVolume() == 175.0);
    CHECK(result->entries().size() == 3);
    CHECK(baseline->totalVolume() == 175.0);
}

TEST_CASE("TrafficGrowth global scaling", "[traffic-growth]") {
    auto baseline = DemandVector::create(3, {
        {0, 1, 100.0},
        {1, 2, 50.0}
    });

    REQUIRE(baseline);

    TrafficGrowth scenario;
    REQUIRE(scenario.setGlobalScale(1.5));

    auto modified = scenario.transform(baseline.value());

    REQUIRE(modified);
    CHECK(modified->totalVolume() == 225.0);
    CHECK(modified->entries()[0].volume == 150.0);
    CHECK(modified->entries()[1].volume == 75.0);

    CHECK(baseline->totalVolume() == 150.0);
}

TEST_CASE("TrafficGrowth commodity scaling", "[traffic-growth]") {
    auto baseline = DemandVector::create(3, {
        {0, 1, 100.0},
        {1, 0, 25.0},
        {1, 2, 50.0}
    });

    REQUIRE(baseline);

    TrafficGrowth scenario;
    REQUIRE(scenario.growPerCommodity(0, 1, 2.0));

    auto modified = scenario.transform(baseline.value());

    REQUIRE(modified);
    CHECK(modified->totalVolume() == 275.0);
    CHECK(modified->entries()[0].volume == 200.0);
    CHECK(modified->entries()[1].volume == 25.0);
    CHECK(modified->entries()[2].volume == 50.0);
}

TEST_CASE("TrafficGrowth combined transformations", "[traffic-growth]") {
    auto baseline = DemandVector::create(3, {
        {0, 1, 100.0},
        {1, 2, 50.0}
    });

    REQUIRE(baseline);

    TrafficGrowth scenario;
    REQUIRE(scenario.setGlobalScale(1.5));
    REQUIRE(scenario.growPerCommodity(0, 1, 2.0));

    auto modified = scenario.transform(baseline.value());

    REQUIRE(modified);

    // 100 * 1.5 * 2 = 300
    //  50 * 1.5     = 75
    CHECK(modified->entries()[0].volume == 300.0);
    CHECK(modified->entries()[1].volume == 75.0);
    CHECK(modified->totalVolume() == 375.0);
}

TEST_CASE("TrafficGrowth rejects invalid configuration", "[traffic-growth]") {
    TrafficGrowth scenario;

    CHECK_FALSE(scenario.setGlobalScale(0.0));
    CHECK_FALSE(scenario.setGlobalScale(-1.0));
    CHECK_FALSE(scenario.setGlobalScale(
        std::numeric_limits<double>::infinity()
    ));

    CHECK_FALSE(scenario.growPerCommodity(0, 0, 2.0));
    CHECK_FALSE(scenario.growPerCommodity(-1, 2, 2.0));
    CHECK_FALSE(scenario.growPerCommodity(0, 1, 0.0));
    CHECK_FALSE(scenario.growPerCommodity(0, 1, -1.0));

    auto baseline = DemandVector::create(
        3, {{0, 1, 100.0}}
    );
    REQUIRE(baseline);

    REQUIRE(scenario.growPerCommodity(1, 2, 2.0));

    auto result = scenario.transform(baseline.value());
    CHECK_FALSE(result);
}

TEST_CASE("TrafficGrowth replacing an override", "[traffic-growth]") {
    auto baseline = DemandVector::create(
        3, {{0, 1, 100.0}}
    );

    REQUIRE(baseline);

    TrafficGrowth scenario;
    REQUIRE(scenario.growPerCommodity(0, 1, 2.0));
    REQUIRE(scenario.growPerCommodity(0, 1, 1.25));

    auto modified = scenario.transform(baseline.value());

    REQUIRE(modified);
    CHECK(modified->totalVolume() == 125.0);
}

