//
// Created by Mert on 30.10.26.
//

#include "io/cli/validator.h"

Result<void> validateCommandCompatibility(const CliValidationContext& ctx) {
    switch (ctx.action) {
    case Action::AnalyzeN1: {
        if (ctx.solver_count != 1) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "analyze-n1 requires exactly one --solver.");
        }

        if (ctx.demand_count != 1) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "analyze-n1 requires exactly one --demand model.");
        }

        if (!ctx.valid_json_output) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "analyze-n1 requires --output <filename.json>.");
        }

        if (ctx.has_load_routing_table || ctx.has_save_routing_table || ctx.has_link_failures || ctx.has_visualization || ctx.failure_analysis || ctx.failure_recovery) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "analyze-n1 does not accept routing-table, "
                                                                "link-failure, visualization, "
                                                                "or legacy failure-analysis flags.");
        }

        return {};
    }

    case Action::AnalyzeTraffic: {
        if (ctx.solver_count != 1 || ctx.demand_count != 1) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "analyze-traffic requires exactly one "
                                                                "solver and demand model.");
        }

        if (!ctx.has_routing_table_input) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "analyze-traffic requires --load-routing-table.");
        }

        if (!ctx.valid_json_output) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "analyze-traffic requires --output <filename.json>.");
        }

        if (ctx.has_link_failures || ctx.failure_analysis || ctx.failure_recovery || ctx.has_save_routing_table || ctx.has_visualization || ctx.has_execution_mode) {

            return makeErrorMessage(ErrorCode::InvalidArgument, "analyze-traffic does not accept failure, "
                                                                "routing-export, visualization, or execution-mode options.");
        }

        return {};
    }

    case Action::AnalyzeScenario: {
        if (!ctx.analysis_kind) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "analyze-scenario requires a valid --analysis.");
        }

        switch (*ctx.analysis_kind) {
        case AnalysisKind::LinkFailure:
            if (ctx.link_failure_count != 1) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "link-failure analysis requires exactly "
                                                                    "one --link-failure edge ID.");
            }
            break;

        case AnalysisKind::CapacityComparison:
            if (ctx.capacity_upgrade_count == 0 || ctx.capacity_upgrade_count > 100) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "capacity analysis requires 1 to 100 capacity upgrades.");
            }

            if (ctx.link_failure_count != 0) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "capacity analysis does not accept --link-failure.");
            }

            break;
        }

        return {};
    }

    case Action::Run: {
        if (ctx.has_capacity_upgrades || ctx.has_link_failures) { // TODO: This code will die...
            return makeErrorMessage(ErrorCode::InvalidArgument, "Legacy scenario execution through 'solve' "
                                                                "is no longer supported. Use 'analyze-scenario' instead.");
        }

        return {};
    }

    case Action::Help:
    case Action::ListSolvers:
        return {};
    }

    return makeErrorMessage(ErrorCode::InvalidArgument, "Unsupported CLI action.");
}
