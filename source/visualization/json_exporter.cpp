//
// Created by Mert Biyikli on 24.07.26.
//

#include "visualization/json_exporter.h"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>


std::string escapeJson(const std::string& value) {
    std::ostringstream output;

    for (const unsigned char character : value) {
        switch (character) {
            case '"':
                output << "\\\"";
                break;

            case '\\':
                output << "\\\\";
                break;

            case '\b':
                output << "\\b";
                break;

            case '\f':
                output << "\\f";
                break;

            case '\n':
                output << "\\n";
                break;

            case '\r':
                output << "\\r";
                break;

            case '\t':
                output << "\\t";
                break;

            default:
                if (character < 0x20) {
                    output
                        << "\\u"
                        << std::hex
                        << std::setw(4)
                        << std::setfill('0')
                        << static_cast<int>(character)
                        << std::dec
                        << std::setfill(' ');
                } else {
                    output << static_cast<char>(character);
                }
        }
    }

    return output.str();
}


std::string RoutingVisualizationJsonExporter::serialize(const RoutingVisualizationResult& result) {
    std::ostringstream output;

    /*
     * Keep enough precision to round-trip double values.
     */
    output << std::setprecision(17);

    output << "{\n";

    output << "  \"schemaVersion\": 1,\n";

    output
        << "  \"graphName\": \""
        << escapeJson(result.graph_name)
        << "\",\n";

    output
        << "  \"solver\": \""
        << escapeJson(result.solver_name)
        << "\",\n";

    output
        << "  \"demandModel\": \""
        << escapeJson(result.demand_model)
        << "\",\n";

    output << "  \"summary\": {\n";

    output
        << "    \"maximumCongestion\": "
        << result.summary.maximum_congestion
        << ",\n";

    output
        << "    \"averageUtilization\": "
        << result.summary.average_utilization
        << ",\n";

    output
        << "    \"overloadedEdges\": "
        << result.summary.overloaded_edges
        << ",\n";

    output
        << "    \"totalRoutedDemand\": "
        << result.summary.total_routed_demand
        << "\n";

    output << "  },\n";

    output << "  \"nodes\": [\n";

    for (
        std::size_t index = 0;
        index < result.nodes.size();
        ++index
    ) {
        const auto& node = result.nodes[index];

        output
            << "    {"
            << "\"id\": " << node.id
            << ", \"label\": \""
            << escapeJson(node.label)
            << "\""
            << "}";

        output
            << (
                index + 1 == result.nodes.size()
                    ? "\n"
                    : ",\n"
            );
    }

    output << "  ],\n";

    output << "  \"edges\": [\n";

    for (std::size_t index = 0;index < result.edges.size();++index) {
        const auto& edge = result.edges[index];

        output
            << "    {"
            << "\"id\": " << edge.id
            << ", \"source\": " << edge.source
            << ", \"target\": " << edge.target
            << ", \"capacity\": " << edge.capacity
            << ", \"load\": " << edge.load
            << ", \"utilization\": " << edge.utilization
            << ", \"enabled\": "
            << (edge.enabled ? "true" : "false")
            << "}";

        output
            << (
                index + 1 == result.edges.size()
                    ? "\n"
                    : ",\n"
            );
    }

    output << "  ],\n";
    output << "  \"linkFailures\": [\n";

    for (
        std::size_t failure_index = 0;
        failure_index < result.link_failures.size();
        ++failure_index
    ) {
        const auto& failure =
            result.link_failures[failure_index];

        output
            << "    {"
            << "\"failedEdgeId\": "
            << failure.failed_edge_id

            << ", \"antiEdgeId\": "
            << failure.anti_edge_id

            << ", \"source\": "
            << failure.source

            << ", \"target\": "
            << failure.target

            << ", \"capacity\": "
            << failure.capacity

            << ", \"affectedCommodities\": "
            << failure.affected_commodities

            << ", \"affectedDemandVolume\": "
            << failure.affected_demand_volume

            << ", \"lostTraffic\": "
            << failure.lost_traffic

            << ", \"totalDemand\": "
            << failure.total_demand

            << ", \"affectedDemandFraction\": "
            << failure.affected_demand_fraction

            << ", \"lostTrafficFraction\": "
            << failure.lost_traffic_fraction

            << "}";

        output
            << (
                failure_index + 1 ==
                        result.link_failures.size()
                    ? "\n"
                    : ",\n"
            );
    }

    output
    << "  ]\n"
    << "}\n";
    return output.str();
}

Result<void> RoutingVisualizationJsonExporter::write(const RoutingVisualizationResult& result,const std::filesystem::path& output_path) {
    std::error_code error;

    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path(),error);

        if (error) {
            return makeErrorMessage(ErrorCode::RuntimeError,"Could not create visualization output directory: " +error.message());
        }
    }

    std::ofstream output(output_path);

    if (!output) {
        return makeErrorMessage(ErrorCode::RuntimeError,"Could not open visualization output file: " +output_path.string());
    }

    output << serialize(result);

    if (!output) {
        return makeErrorMessage(ErrorCode::RuntimeError,"Could not write visualization output file: " +output_path.string());
    }

    return {};
}