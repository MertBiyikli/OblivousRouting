//
// Created by Mert on 29.10.26.
//

#ifndef E_ROUTING_TRAFFIC_GROWTH_RESULT_H
#define E_ROUTING_TRAFFIC_GROWTH_RESULT_H

#include "analysis/scenario_result_json.h"
#include "analysis/scenario_results.h"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <optional>
#include <string>

struct TrafficCommodityOverride {
    int source = 0;
    int target = 0;
    double factor = 1.0;
};

struct TrafficGrowthResult {
    std::uint64_t graph_fingerprint = 0;

    // Both evaluations use the identical routing scheme and graph.
    NetworkData baseline;
    NetworkData modified;

    double demand_scale = 1.0;
    double absolute_demand_change = 0.0;
    std::optional<double> relative_demand_change;

    double absolute_congestion_change = 0.0;
    std::optional<double> relative_congestion_change;

    // User-requested transformation, distinct from the
    // observed total-demand ratio.
    double requested_global_scale = 1.0;

    std::vector<TrafficCommodityOverride> commodity_overrides;
};

class TrafficGrowthResultJsonWriter {
  public:
    static Result<void> write(const TrafficGrowthResult& result, const std::filesystem::path& path) {
        try {
            if (!path.parent_path().empty()) {
                std::filesystem::create_directories(path.parent_path());
            }

            std::ofstream out(path, std::ios::trunc);

            if (!out) {
                return makeErrorMessage(ErrorCode::RuntimeError, "Cannot open traffic analysis output.");
            }

            out.imbue(std::locale::classic());
            out << std::setprecision(std::numeric_limits<double>::max_digits10);

            out << "{\"schema_version\":\"1.0\"";
            out << ",\"scenario_kind\":\"traffic_growth\"";
            out << ",\"execution_mode\":\"fixed_routing\"";

            out << ",\"graph_fingerprint\":";
            quoted(out, std::to_string(result.graph_fingerprint));

            out << ",\"baseline\":";
            network(out, std::optional<NetworkData>{result.baseline});

            out << ",\"modified\":";
            network(out, std::optional<NetworkData>{result.modified});

            out << ",\"demand_scale\":";
            number(out, result.demand_scale);

            out << ",\"requested_global_scale\":";
            number(out, result.requested_global_scale);

            out << ",\"commodity_overrides\":[";

            for (std::size_t i = 0; i < result.commodity_overrides.size(); ++i) {

                if (i != 0) {
                    out << ',';
                }

                const auto& item = result.commodity_overrides[i];

                out << "{\"source\":" << item.source << ",\"target\":" << item.target << ",\"factor\":";

                number(out, item.factor);

                out << '}';
            }

            out << ']';

            out << ",\"absolute_demand_change\":";
            number(out, result.absolute_demand_change);

            out << ",\"relative_demand_change\":";
            optionalNumber(out, result.relative_demand_change);

            out << ",\"absolute_congestion_change\":";
            number(out, result.absolute_congestion_change);

            out << ",\"relative_congestion_change\":";
            optionalNumber(out, result.relative_congestion_change);

            out << "}\n";
            out.flush();

            if (!out) {
                return makeErrorMessage(ErrorCode::RuntimeError, "Failed writing traffic analysis output.");
            }

            return {};
        } catch (const std::exception& error) {
            return makeErrorMessage(ErrorCode::RuntimeError, std::string("Traffic JSON serialization failed: ") + error.what());
        }
    }
};

#endif // E_ROUTING_TRAFFIC_GROWTH_RESULT_H
