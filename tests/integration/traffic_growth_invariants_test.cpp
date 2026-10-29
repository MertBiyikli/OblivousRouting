//
// Created by Mert on 29.10.26.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "analysis/traffic_change/traffic_growth.h"
#include "analysis/traffic_change/demand_vector.h"

#include <cmath>
#include <map>
#include <utility>

namespace {

using Commodity = std::pair<int, int>;
using Volumes = std::map<Commodity, double>;

DemandVector makeBaseline() {
    auto result = DemandVector::create(4, {
        {0, 1, 10.0},
        {0, 2, 20.0},
        {1, 2, 30.0},
        {2, 3, 40.0}
    });

    REQUIRE(result.has_value());
    return std::move(result.value());
}

Volumes volumesOf(const DemandVector& vector) {
    Volumes volumes;

    for (const auto& entry : vector.entries()) {
        volumes[{entry.source, entry.target}] += entry.volume;
    }

    return volumes;
}

void requireVolumes(
    const DemandVector& actual,
    const Volumes& expected
) {
    const auto volumes = volumesOf(actual);

    REQUIRE(volumes.size() == expected.size());

    for (const auto& [commodity, expectedVolume] : expected) {
        const auto it = volumes.find(commodity);

        REQUIRE(it != volumes.end());
        CHECK(it->second == Catch::Approx(expectedVolume).epsilon(1e-10));
    }
}

} // namespace

TEST_CASE(
    "TrafficGrowth identity preserves every commodity",
    "[traffic-growth][invariants]"
) {
    const auto baseline = makeBaseline();

    TrafficGrowth scenario;
    REQUIRE(scenario.setGlobalScale(1.0));

    const auto result = scenario.transform(baseline);
    REQUIRE(result.has_value());

    CHECK(result->nodeCount() == baseline.nodeCount());

    CHECK(result->totalVolume() ==
          Catch::Approx(baseline.totalVolume()));

    requireVolumes(result.value(), {
        {{0, 1}, 10.0},
        {{0, 2}, 20.0},
        {{1, 2}, 30.0},
        {{2, 3}, 40.0}
    });
}

TEST_CASE(
    "TrafficGrowth uniform scaling is linear",
    "[traffic-growth][invariants]"
) {
    const auto baseline = makeBaseline();

    for (const double factor : {0.5, 1.0, 1.2, 2.0}) {
        DYNAMIC_SECTION("Global scale = " << factor) {
            TrafficGrowth scenario;

            REQUIRE(scenario.setGlobalScale(factor));

            const auto result = scenario.transform(baseline);
            REQUIRE(result.has_value());

            CHECK(result->totalVolume() ==
                Catch::Approx(100.0 * factor));

            const auto original = volumesOf(baseline);
            const auto modified = volumesOf(result.value());

            REQUIRE(original.size() == modified.size());

            for (const auto& [commodity, volume] : original) {
                REQUIRE(modified.contains(commodity));

                CHECK(modified.at(commodity) ==
                    Catch::Approx(volume * factor));
            }
        }
    }
}

TEST_CASE(
    "TrafficGrowth commodity scaling isolates the selected OD pair",
    "[traffic-growth][invariants]"
) {
    const auto baseline = makeBaseline();

    TrafficGrowth scenario;

    REQUIRE(scenario.setGlobalScale(1.0));
    REQUIRE(scenario.growPerCommodity(0, 2, 3.0));

    const auto result = scenario.transform(baseline);
    REQUIRE(result.has_value());

    requireVolumes(result.value(), {
        {{0, 1}, 10.0},
        {{0, 2}, 60.0},
        {{1, 2}, 30.0},
        {{2, 3}, 40.0}
    });

    CHECK(result->totalVolume() == Catch::Approx(140.0));
}

TEST_CASE(
    "TrafficGrowth combined scaling is deterministic",
    "[traffic-growth][invariants]"
) {
    const auto baseline = makeBaseline();

    TrafficGrowth first;
    TrafficGrowth second;

    REQUIRE(first.setGlobalScale(1.2));
    REQUIRE(first.growPerCommodity(0, 2, 2.0));

    REQUIRE(second.setGlobalScale(1.2));
    REQUIRE(second.growPerCommodity(0, 2, 2.0));

    const auto a = first.transform(baseline);
    const auto b = second.transform(baseline);

    REQUIRE(a.has_value());
    REQUIRE(b.has_value());

    const auto aVolumes = volumesOf(a.value());
    const auto bVolumes = volumesOf(b.value());

    REQUIRE(aVolumes.size() == bVolumes.size());

    for (const auto& [commodity, volume] : aVolumes) {
        REQUIRE(bVolumes.contains(commodity));
        CHECK(bVolumes.at(commodity) ==
            Catch::Approx(volume));
    }

    CHECK(a->totalVolume() ==
        Catch::Approx(b->totalVolume()));
}
