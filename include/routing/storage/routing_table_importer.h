//
// Created by Mert on 15.10.26.
//

#ifndef E_ROUTING_ROUTING_TABLE_IMPORTER_H
#define E_ROUTING_ROUTING_TABLE_IMPORTER_H
#include "core/errors.h"
#include "routing/storage/routing_table_adapter.h"
#include "data_structures/graph/graph.h"

#include <cstdint>
#include <memory>
#include <string>

class RoutingTableImporter {
public:
    static Result<std::unique_ptr<RoutingScheme>> load(const optimized::Graph<EdgeData>& graph, const std::string& path) {
        try {
            const std::uint64_t fingerprint = optimized::GraphFingerprint::compute(graph);

            // Inspect the version-2 header and determine the table representation.
            // Do not guess the type from the solver name.
            auto header = RoutingTableSerializer::inspect(path);

            if (!header) {
                return getError(header);
            }

            if (header->type == RoutingTableType::Linear) {
                std::int32_t root = -1;
                auto table = RoutingTableAdapter::loadLinear(graph, path, fingerprint, &root);

                if (!table) {
                    return getError(table);
                }

                if (root < 0 || root >= graph.getNumNodes()) {
                    return makeErrorMessage(ErrorCode::InvalidRouting, "Invalid stored linear routing root");
                }

                std::unique_ptr<RoutingScheme> scheme = std::make_unique<LinearRoutingScheme>(graph, root, std::move(table.value()));
                return scheme;
            }

            if (header->type == RoutingTableType::AllPair) {
                auto table = RoutingTableAdapter::loadAllPair(graph, path, fingerprint);

                if (!table) {
                    return getError(table);
                }

                std::unique_ptr<RoutingScheme> scheme = std::make_unique<AllPairRoutingScheme>(graph, std::move(table.value()));
                return scheme;
            }

            return makeErrorMessage(ErrorCode::InvalidRouting, "Unsupported stored routing-table type");
        } catch (const std::exception& error) {
            return fromStdException(error, ErrorCode::RuntimeError);
        }
    }
};

#endif // E_ROUTING_ROUTING_TABLE_IMPORTER_H
