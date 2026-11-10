//
// Created by Mert on 30.10.26.
//

#include <catch2/catch_test_macros.hpp>

#include "io/cli/validator.h"

TEST_CASE("CLI validator accepts valid N-1 configuration", "[cli][validation]") {
    CliValidationContext ctx{};
    ctx.action = Action::AnalyzeN1;
    ctx.solver_count = 1;
    ctx.demand_count = 1;
    ctx.valid_json_output = true;

    REQUIRE(validateCommandCompatibility(ctx).has_value());
}

TEST_CASE("CLI validator rejects multiple N-1 solvers", "[cli][validation]") {
    CliValidationContext ctx{};
    ctx.action = Action::AnalyzeN1;
    ctx.solver_count = 2;
    ctx.demand_count = 1;
    ctx.valid_json_output = true;

    REQUIRE_FALSE(validateCommandCompatibility(ctx).has_value());
}

TEST_CASE("CLI validator rejects N-1 link failures", "[cli][validation]") {
    CliValidationContext ctx{};
    ctx.action = Action::AnalyzeN1;
    ctx.solver_count = 1;
    ctx.demand_count = 1;
    ctx.valid_json_output = true;
    ctx.has_link_failures = true;

    REQUIRE_FALSE(validateCommandCompatibility(ctx).has_value());
}

TEST_CASE("CLI validator requires traffic routing table", "[cli][validation]") {
    CliValidationContext ctx{};
    ctx.action = Action::AnalyzeTraffic;
    ctx.solver_count = 1;
    ctx.demand_count = 1;
    ctx.valid_json_output = true;

    REQUIRE_FALSE(validateCommandCompatibility(ctx).has_value());
}

TEST_CASE("CLI validator accepts traffic with routing table", "[cli][validation]") {
    CliValidationContext ctx{};
    ctx.action = Action::AnalyzeTraffic;
    ctx.solver_count = 1;
    ctx.demand_count = 1;
    ctx.valid_json_output = true;
    ctx.has_routing_table_input = true;
    ctx.has_load_routing_table = true;

    REQUIRE(validateCommandCompatibility(ctx).has_value());
}

TEST_CASE("CLI validator rejects legacy solve capacity upgrades", "[integration][cli]") {
    CliValidationContext ctx;
    ctx.action = Action::Run;
    ctx.has_capacity_upgrades = true;

    auto result = validateCommandCompatibility(ctx);

    REQUIRE_FALSE(result);
}

TEST_CASE("CLI validator rejects legacy solve link failures", "[integration][cli]") {
    CliValidationContext ctx;
    ctx.action = Action::Run;
    ctx.has_link_failures = true;

    auto result = validateCommandCompatibility(ctx);

    REQUIRE_FALSE(result);
}

TEST_CASE("CLI validator accepts ordinary solve", "[integration][cli]") {
    CliValidationContext ctx;
    ctx.action = Action::Run;

    auto result = validateCommandCompatibility(ctx);

    REQUIRE(result);
}
