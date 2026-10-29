//
// Created by Mert on 27.10.26.
//

#ifndef E_ROUTING_SCENARIO_RESULT_JSON_H
#define E_ROUTING_SCENARIO_RESULT_JSON_H

#include "analysis/scenario_results.h"
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <ostream>
#include <stdexcept>
#include <string>
#include "analysis/scenario_results.h"
#include "core/errors.h"



#include <boost/multiprecision/number.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>
#include <system_error>

// Dedicated scenario analysis JSON; intentionally does NOT change RoutingResultWriter.
namespace scenario_json_detail {

inline void quoted(std::ostream& os, const std::string& value) {
    constexpr char hex[] = "0123456789abcdef";
    os << '"';
    for (unsigned char c : value) {
        switch (c) {
        case '"':
            os << "\\\"";
            break;
        case '\\':
            os << "\\\\";
            break;
        case '\b':
            os << "\\b";
            break;
        case '\f':
            os << "\\f";
            break;
        case '\n':
            os << "\\n";
            break;
        case '\r':
            os << "\\r";
            break;
        case '\t':
            os << "\\t";
            break;
        default:
            if (c < 0x20)
                os << "\\u00" << hex[c >> 4] << hex[c & 15];
            else
                os << static_cast<char>(c);
        }
    }
    os << '"';
}
inline void number(std::ostream& os, double x) {
    if (!std::isfinite(x))
        throw std::runtime_error("Nonfinite scenario numeric field");
    os << x;
}
inline void optionalNumber(std::ostream& os, const std::optional<double>& value) {
    if (!value)
        os << "null";
    else
        number(os, *value);
}
inline void network(std::ostream& os, const std::optional<NetworkData>& data) {
    if (!data) {
        os << "null";
        return;
    }
    os << "{\"max_congestion\":";
    number(os, data->max_congestion);
    os << ",\"average_congestion\":";
    number(os, data->average_congestion);
    os << ",\"total_demand_volume\":";
    number(os, data->total_demand_volume);
    os << ",\"links\":[";
    for (std::size_t i = 0; i < data->links.size(); ++i) {
        if (i)
            os << ',';
        const auto& l = data->links[i];
        os << "{\"edge_id\":" << l.edge_id << ",\"capacity\":";
        number(os, l.capacity);
        os << ",\"load\":";
        number(os, l.load);
        os << ",\"utilization\":";
        number(os, l.utilization);
        os << '}';
    }
    os << "]}";
}
inline void performance(std::ostream& os, const std::optional<ScenarioPerformanceMetrics>& data) {
    if (!data) {
        os << "null";
        return;
    }
    const auto& m = *data;
    os << "{\"baseline_congestion\":";
    number(os, m.baseline_congestion);
    os << ",\"fixed_routing_congestion\":";
    number(os, m.fixed_routing_congestion);
    os << ",\"reoptimized_congestion\":";
    number(os, m.reoptimized_congestion);
    os << ",\"capacity_only_improvement\":";
    optionalNumber(os, m.capacity_only_improvement);
    os << ",\"total_improvement\":";
    optionalNumber(os, m.total_improvement);
    os << ",\"reoptimization_gain\":";
    optionalNumber(os, m.reoptimization_gain);
    os << ",\"total_absolute_load_redistribution\":";
    optionalNumber(os, m.total_absolute_load_redistribution);
    os << ",\"maximum_absolute_link_load_change\":";
    optionalNumber(os, m.maximum_absolute_link_load_change);
    os << ",\"changed_physical_links\":";
    if (m.changed_physical_links)
        os << *m.changed_physical_links;
    else
        os << "null";
    os << '}';
}
inline const char* kind(ScenarioKind k) {
    switch (k) {
    case ScenarioKind::GENERIC:
        return "generic";
    case ScenarioKind::CAPACITY_UPGRADE:
        return "capacity_upgrade";
    case ScenarioKind::LINK_FAILURE:
        return "link_failure";
    case ScenarioKind::TRAFFIC_GROWTH:
        return "traffic_growth";
    }
    return "unknown";
}
inline const char* mode(ExecutionMode m) {
    switch (m) {
    case ExecutionMode::FAST:
        return "fast";
    case ExecutionMode::STANDARD:
        return "standard";
    case ExecutionMode::RESEARCH:
        return "research";
    }
    return "unknown";
}

inline void comparison(
    std::ostream& os,
    const std::optional<ComparisonData>& data
) {
    if (!data) {
        os << "null";
        return;
    }

    os << "{\"congestion_before\":";
    number(os, data->congestion_before);

    os << ",\"congestion_after\":";
    number(os, data->congestion_after);

    os << ",\"absolute_reduction\":";
    number(os, data->absolute_reduction);

    os << ",\"relative_reduction\":";
    optionalNumber(os, data->relative_reduction);

    os << '}';
}
} // namespace scenario_json_detail



class ScenarioResultJsonWriter {
  public:
    // May throw for an invalid/non-finite scalar; use write() to get Result<void>.
    static void writeTo(std::ostream& os, const ScenarioResult& result) {
        using namespace scenario_json_detail;
        const auto old_flags = os.flags();
        const auto old_precision = os.precision();
        const auto old_locale = os.getloc();
        os.imbue(std::locale::classic());
        os << std::setprecision(std::numeric_limits<double>::max_digits10);
        os << "{\"schema_version\":\"1.1\",\"scenario_kind\":";
        quoted(os, kind(result.kind));
        os << ",\"execution_mode\":";
        quoted(os, mode(result.execution_mode));
        os << ",\"original_fingerprint\":";
        quoted(os, std::to_string(result.original_fingerprint));

        os << ",\"modified_fingerprint\":";
        quoted(os, std::to_string(result.modified_fingerprint));
        os << ",\"capacity_upgrades\":[";
        for (std::size_t i = 0; i < result.capacity_upgrades.size(); ++i) {
            if (i)
                os << ',';
            const auto& u = result.capacity_upgrades[i];
            os << "{\"edge_id\":" << u.edge_id << ",\"original_capacity\":";
            number(os, u.original_capacity);
            os << ",\"upgraded_capacity\":";
            number(os, u.upgraded_capacity);
            os << ",\"capacity_factor\":";
            number(os, u.capacity_factor);
            os << '}';
        }


        os << "],\"link_failures\":[";

        for (std::size_t i = 0;
             i < result.link_failures.size(); ++i) {

            if (i)
                os << ',';

            os << result.link_failures[i];
             }

        os << "],\"solver_results\":[";
        for (std::size_t i = 0; i < result.solver_results.size(); ++i) {
            if (i)
                os << ',';
            const auto& s = result.solver_results[i];
            os << "{\"solver_algorithms\":[";
            for (std::size_t j = 0; j < s.solver_algorithms.size(); ++j) {
                if (j)
                    os << ',';
                quoted(os, s.solver_algorithms[j]);
            }
            os << "],\"baseline_solver_seconds\":";
            number(os, s.baseline_solver_seconds);
            os << ",\"modified_solver_seconds\":";
            number(os, s.modified_solver_seconds);
            os << ",\"baseline_evaluation_seconds\":";
            number(os, s.baseline_evaluation_seconds);
            os << ",\"modified_evaluation_seconds\":";
            number(os, s.modified_evaluation_seconds);
            os << ",\"baseline\":";
            network(os, s.baseline_Data);
            os << ",\"fixed_routing\":";
            network(os, s.fixed_routing_Data);
            os << ",\"reoptimized\":";
            network(os, s.modified_Data);

            os << ",\"comparison\":";
            comparison(os, s.comparison);

            os << ",\"fixed_routing_comparison\":";
            comparison(os, s.fixed_routing_comparison);

            os << ",\"performance_metrics\":";
            performance(os, s.performance_metrics);

            os << '}';
        }
        os << "]}";
        os.flags(old_flags);
        os.precision(old_precision);
        os.imbue(old_locale);
    }

    static Result<void> write(const ScenarioResult& result, const std::filesystem::path& destination) {
        try {
            if (!destination.parent_path().empty()) {
                std::filesystem::create_directories(destination.parent_path());
            }
            std::ofstream out(destination, std::ios::trunc);
            if (!out) {
                return makeErrorMessage(ErrorCode::RuntimeError, "Cannot open scenario JSON output: " + destination.string());
            }
            writeTo(out, result);
            out << '\n';
            out.flush();
            if (!out) {
                return makeErrorMessage(ErrorCode::RuntimeError, "Failed writing scenario JSON output: " + destination.string());
            }
            return {};
        } catch (const std::exception& e) {
            return makeErrorMessage(ErrorCode::RuntimeError, "Scenario JSON serialization failed: " + std::string(e.what()));
        }
    }
};



#endif // E_ROUTING_SCENARIO_RESULT_JSON_H
