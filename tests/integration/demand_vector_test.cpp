//
// Created by Mert on 29.10.26.
//

#include <catch2/catch_test_macros.hpp>

#include "analysis/traffic_change/demand_vector.h"

#include <cmath>
#include <limits>
#include <vector>

TEST_CASE("DemandVector normalization", "[traffic][demand-vector]") {
    const auto result = DemandVector::create(4, {
        {0, 1, 100.0},
        {0, 1, 50.0},
        {1, 0, 25.0},
        {2, 3, 0.0}
    });

    REQUIRE(result);

    const auto& vector = result.value();

    REQUIRE(vector.nodeCount() == 4);
    REQUIRE(vector.entries().size() == 2);
    REQUIRE(vector.totalVolume() == 175.0);

    CHECK(vector.entries()[0].source == 0);
    CHECK(vector.entries()[0].target == 1);
    CHECK(vector.entries()[0].volume == 150.0);

    CHECK(vector.entries()[1].source == 1);
    CHECK(vector.entries()[1].target == 0);
    CHECK(vector.entries()[1].volume == 25.0);
}

TEST_CASE("DemandVector global scaling", "[traffic][demand-vector]") {
    auto baseline = DemandVector::create(3, {
        {0, 1, 100.0},
        {1, 2, 50.0}
    });

    REQUIRE(baseline);

    auto identity = baseline->scaled(1.0);
    auto growth = baseline->scaled(1.5);

    REQUIRE(identity);
    REQUIRE(growth);

    CHECK(identity->totalVolume() == 150.0);
    CHECK(growth->totalVolume() == 225.0);

    CHECK(growth->entries()[0].volume == 150.0);
    CHECK(growth->entries()[1].volume == 75.0);

    // The original data must not change.
    CHECK(baseline->totalVolume() == 150.0);
}

TEST_CASE("DemandVector commodity scaling", "[traffic][demand-vector]") {
    auto baseline = DemandVector::create(3, {
        {0, 1, 100.0},
        {1, 0, 25.0},
        {1, 2, 50.0}
    });

    REQUIRE(baseline);

    auto modified = baseline->scaledCommodity(0, 1, 2.0);
    REQUIRE(modified);

    CHECK(modified->totalVolume() == 275.0);
    CHECK(modified->entries()[0].volume == 200.0);
    CHECK(modified->entries()[1].volume == 25.0);
    CHECK(modified->entries()[2].volume == 50.0);

    CHECK_FALSE(baseline->scaledCommodity(2, 0, 2.0));
}

TEST_CASE("DemandVector engine round trip", "[traffic][demand-vector]") {
    demands raw;
    raw.addDemand(0, 1, 100.0);
    raw.addDemand(0, 1, 50.0);
    raw.addDemand(1, 0, 25.0);

    auto normalized = DemandVector::fromDemands(raw, 3);
    REQUIRE(normalized);

    auto converted = normalized->toDemands();
    REQUIRE(converted);

    CHECK(converted->size() == 2);
    REQUIRE(converted->getDemandValue(0, 1).has_value());
    REQUIRE(converted->getDemandValue(1, 0).has_value());

    CHECK(*converted->getDemandValue(0, 1) == 150.0);
    CHECK(*converted->getDemandValue(1, 0) == 25.0);

    auto roundtrip = DemandVector::fromDemands(
        converted.value(), 3
    );

    REQUIRE(roundtrip);
    CHECK(roundtrip->totalVolume() == 175.0);
}

TEST_CASE("DemandVector rejects invalid data", "[traffic][demand-vector]") {
    const double infinity =
        std::numeric_limits<double>::infinity();

    const double nan =
        std::numeric_limits<double>::quiet_NaN();

    CHECK_FALSE(DemandVector::create(0, {{0, 1, 1.0}}));
    CHECK_FALSE(DemandVector::create(3, {{0, 3, 1.0}}));
    CHECK_FALSE(DemandVector::create(3, {{1, 1, 1.0}}));
    CHECK_FALSE(DemandVector::create(3, {{0, 1, -1.0}}));
    CHECK_FALSE(DemandVector::create(3, {{0, 1, nan}}));
    CHECK_FALSE(DemandVector::create(3, {{0, 1, infinity}}));
    CHECK_FALSE(DemandVector::create(3, {{0, 1, 0.0}}));

    auto baseline = DemandVector::create(
        3, {{0, 1, 100.0}}
    );

    REQUIRE(baseline);

    CHECK_FALSE(baseline->scaled(0.0));
    CHECK_FALSE(baseline->scaled(-1.0));
    CHECK_FALSE(baseline->scaled(infinity));

    CHECK_FALSE(baseline->scaledCommodity(0, 1, 0.0));
    CHECK_FALSE(baseline->scaledCommodity(0, 1, -1.0));
    CHECK_FALSE(baseline->scaledCommodity(0, 1, infinity));

    demands malformed;
    malformed.source.push_back(0);
    malformed.target.push_back(1);
    CHECK_FALSE(DemandVector::fromDemands(malformed, 3));
}

TEST_CASE("DemandVector detects overflow", "[traffic][demand-vector]") {
    const double maximum =
        std::numeric_limits<double>::max();

    CHECK_FALSE(DemandVector::create(3, {
        {0, 1, maximum},
        {0, 1, maximum}
    }));

    CHECK_FALSE(DemandVector::create(3, {
        {0, 1, maximum},
        {1, 2, maximum}
    }));

    auto baseline = DemandVector::create(
        3, {{0, 1, maximum}}
    );

    REQUIRE(baseline);
    CHECK_FALSE(baseline->scaled(2.0));

    auto tiny = DemandVector::create(
    3, {{0, 1, std::numeric_limits<double>::denorm_min()}}
);

    REQUIRE(tiny);

    // Scaling the smallest positive representable value by 0.1
    // rounds to zero, which our implementation rejects.
    CHECK_FALSE(tiny->scaled(0.1));
}
