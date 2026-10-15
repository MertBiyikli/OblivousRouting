//
// Created by Mert on 15.10.26.
//

#ifndef E_ROUTING_ROUTING_TABLE_EXPORTER_H
#define E_ROUTING_ROUTING_TABLE_EXPORTER_H
#include "io/solver_io.h"
#include "routing/routing_result.h"
#include "routing_table_adapter.h"

class RoutingTableExporter {
  public:
    static Result<void> save(const RoutingExperimentResult& experiment, const optimized::Graph<EdgeData>& graph, const std::string& directory) {
        namespace fs = std::filesystem;

        if (directory.empty()) {
            return {};
        }

        std::error_code ec;
        fs::create_directories(directory, ec);

        if (ec) {
            return makeErrorMessage(ErrorCode::RuntimeError, "Cannot create routing table directory: " + ec.message());
        }

        // Temporary placeholder until deterministic graph hashing is implemented.
        const std::uint64_t graphFingerprint = optimized::GraphFingerprint::compute(graph);

        for (std::size_t i = 0; i < experiment.solver_results.size(); ++i) {
            const auto& result = experiment.solver_results[i];

            const std::string filename = std::to_string(i) + "_" + getSolverTypeName(result.type) + ".bin";
            const fs::path path = fs::path(directory) / filename;

            if (result.scheme) {
                auto saved = saveScheme(*result.scheme, graph, path.string(), graphFingerprint);

                if (!saved) {
                    return getError(saved);
                }

                continue;
            }

            // Semi-oblivious solvers can produce one routing scheme per demand model.
            for (const auto& [demandName, scheme] : result.demand_schemes) {
                if (!scheme) {
                    continue;
                }

                const fs::path demandPath = fs::path(directory) / (std::to_string(i) + "_" + getSolverTypeName(result.type) + "_" + demandName + ".bin");

                auto saved = saveScheme(*scheme, graph, demandPath.string(), graphFingerprint);

                if (!saved) {
                    return getError(saved);
                }
            }
        }

        return {};
    }

  private:
    static Result<void> saveScheme(const RoutingScheme& scheme, const optimized::Graph<EdgeData>& graph, const std::string& path, std::uint64_t fingerprint) {
        if (const auto* linear = dynamic_cast<const LinearRoutingScheme*>(&scheme)) {
            return RoutingTableAdapter::saveLinear(linear->routing_table, graph, path, fingerprint, linear->root_x);
        }

        if (const auto* allPairs = dynamic_cast<const AllPairRoutingScheme*>(&scheme)) {
            return RoutingTableAdapter::saveAllPair(allPairs->routing_table, graph, path, fingerprint);
        }

        return makeErrorMessage(ErrorCode::InvalidRouting, "Unsupported routing scheme for serialization: " + path);
    }
};

#endif // E_ROUTING_ROUTING_TABLE_EXPORTER_H
