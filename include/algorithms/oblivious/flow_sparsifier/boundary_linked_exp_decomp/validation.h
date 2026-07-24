//
// Created by Mert Biyikli on 21.07.26.
//

#ifndef OBLIVIOUSROUTING_VALIDATION_H
#define OBLIVIOUSROUTING_VALIDATION_H
#include "expander_decomp.h"

#include <cstddef>
#include <span>
#include <string>
#include <vector>

struct ValidationOptions {
    // Exact enumeration is exponential. Clusters larger than this bound are
    // still checked structurally and against their recorded certificates, but
    // their expansion is not independently recomputed. Zero disables it.
    std::size_t exact_max_vertices{20};
    bool require_exact_for_all{false};
    double tolerance{1e-12};
};

struct ValidationReport {
    bool parameters_valid{true};
    bool partition_valid{true};
    bool cluster_metadata_valid{true};
    bool split_audit_valid{true};
    bool proof_model_valid{true};
    bool certificate_metadata_valid{true};
    bool property3_valid{true};
    bool boundary_accounting_valid{true};
    bool exact_expansion_valid{true};
    bool exact_expansion_complete{true};
    std::size_t exact_clusters_checked{};
    std::size_t exact_clusters_skipped{};
    std::vector<std::string> errors;

    [[nodiscard]] bool valid() const noexcept { return errors.empty(); }
};

// Independently validates a decomposition result. An empty initial_cluster has
// the same meaning as in BoundaryLinkedDecomposer::decompose: all vertices.
// Exact enumeration is intended for tests and small proof-of-concept inputs.
[[nodiscard]] ValidationReport validate_decomposition(
    const Graph& graph,
    const DecompositionConfig& config,
    const DecompositionResult& result,
    std::span<const int> initial_cluster = {},
    ValidationOptions options = {});

#endif //OBLIVIOUSROUTING_VALIDATION_H
