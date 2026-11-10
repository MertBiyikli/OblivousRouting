//
// Created by Mert on 30.10.26.
//

#ifndef E_ROUTING_VALIDATOR_H
#define E_ROUTING_VALIDATOR_H

#include "core/errors.h"
#include "io/cli/command.h"

#include <cstddef>

struct CliValidationContext {
    Action action = Action::Run;

    std::size_t solver_count = 0;
    std::size_t demand_count = 0;

    bool valid_json_output = false;
    bool has_routing_table_input = false;

    bool has_capacity_upgrades = false;
    bool has_link_failures = false;

    bool has_load_routing_table = false;
    bool has_save_routing_table = false;
    bool has_visualization = false;
    bool has_execution_mode = false;

    bool failure_analysis = false;
    bool failure_recovery = false;

    bool has_analysis_option = false;
    bool has_exactly_one_failure_target = false;
    bool has_invalid_failure_target_count = false;

    std::optional<AnalysisKind> analysis_kind;
    std::size_t capacity_upgrade_count = 0;
    std::size_t link_failure_count = 0;
};

Result<void> validateCommandCompatibility(const CliValidationContext& context);

#endif // E_ROUTING_VALIDATOR_H
