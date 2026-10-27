//
// Created by Mert on 22.10.26.
//

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <random>
#include <utility>
#include <vector>

#include "utils/demands.h"

namespace {

constexpr double kTestEps = 1e-9;

// Original algorithm: process every directed demand.
double referenceLoad(const std::vector<double>& basis, const demands& demandMap, double eps) {
    double load = 0.0;

    for (std::size_t i = 0; i < demandMap.size(); ++i) {
        const int s = demandMap.source[i];
        const int t = demandMap.target[i];
        const double demand = demandMap.demand_values[i];

        if (s == t || demand <= 0.0) {
            continue;
        }

        const double magnitude = std::abs(basis[s] - basis[t]);

        if (magnitude < eps) {
            continue;
        }

        load += demand * magnitude;
    }

    return load;
}

template <typename CombinedDemands>
double combinedLoad(const std::vector<double>& basis, const CombinedDemands& combined, double eps) {
    double load = 0.0;

    for (const auto& demand : combined) {
        const double magnitude = std::abs(basis[demand.source] - basis[demand.target]);

        if (magnitude < eps) {
            continue;
        }

        load += demand.value * magnitude;
    }

    return load;
}

void compareLoads(const std::vector<double>& basis, const demands& demandMap, double eps = kTestEps) {
    const double original = referenceLoad(basis, demandMap, eps);
    auto combined = combineDemandPairs(demandMap, static_cast<int>(basis.size()));
    const double optimized = combinedLoad(basis, combined, eps);

    REQUIRE(std::isfinite(original));
    REQUIRE(std::isfinite(optimized));

    CHECK(optimized == Catch::Approx(original).epsilon(1e-10));
}

} // namespace

TEST_CASE("Combined demands: symmetric traffic", "[routing][combined-demands]") {
    demands d;
    d.addDemand(0, 1, 100.0);
    d.addDemand(1, 0, 100.0);
    d.addDemand(0, 2, 200.0);
    d.addDemand(2, 0, 200.0);

    compareLoads({0.1, 0.4, -0.2}, d);
}

TEST_CASE("Combined demands: asymmetric traffic", "[routing][combined-demands]") {
    demands d;
    d.addDemand(0, 1, 100.0);
    d.addDemand(1, 0, 350.0);
    d.addDemand(0, 2, 20.0);
    d.addDemand(2, 0, 900.0);
    d.addDemand(1, 2, 55.0);

    compareLoads({-2.0, 0.5, 3.0}, d);
}

TEST_CASE("Combined demands: sparse and zero traffic", "[routing][combined-demands]") {
    demands d;
    d.addDemand(0, 3, 12.0);
    d.addDemand(3, 0, 0.0);
    d.addDemand(1, 2, 0.0);
    d.addDemand(2, 4, 7.0);
    d.addDemand(3, 3, 100.0);

    compareLoads({0.0, 1.5, -0.3, 2.0, 3.2}, d);
}

TEST_CASE("Combined demands: EPS threshold", "[routing][combined-demands]") {
    demands d;
    d.addDemand(0, 1, 100.0);
    d.addDemand(1, 0, 200.0);
    d.addDemand(0, 2, 50.0);
    d.addDemand(2, 1, 75.0);

    compareLoads({0.0, 0.25e-9, 2.0e-9}, d, kTestEps);
}

TEST_CASE("Combined demands: randomized dense matrix", "[routing][combined-demands]") {
    constexpr int n = 100;
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> demandDist(0.0, 500.0);
    std::uniform_real_distribution<double> basisDist(-1.0, 1.0);

    demands d;

    for (int s = 0; s < n; ++s) {
        for (int t = 0; t < n; ++t) {
            if (s != t) {
                d.addDemand(s, t, demandDist(rng));
            }
        }
    }

    for (int trial = 0; trial < 10; ++trial) {
        std::vector<double> basis(n);

        for (double& value : basis) {
            value = basisDist(rng);
        }

        compareLoads(basis, d);
    }
}

TEST_CASE("Combined demands: identical basis values", "[routing][combined-demands]") {
    demands d;
    d.addDemand(0, 1, 100.0);
    d.addDemand(1, 0, 250.0);
    d.addDemand(1, 2, 500.0);

    const std::vector<double> basis{2.0, 2.0, 2.0};

    CHECK(referenceLoad(basis, d, kTestEps) == 0.0);
    auto combinePairs = combineDemandPairs(d, 3);
    CHECK(combinedLoad(basis, combinePairs, kTestEps) == 0.0);
}
