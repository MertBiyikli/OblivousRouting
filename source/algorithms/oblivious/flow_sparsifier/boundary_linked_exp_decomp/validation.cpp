//
// Created by Mert Biyikli on 21.07.26.
//

#include "algorithms/oblivious/flow_sparsifier/boundary_linked_exp_decomp/validation.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <numeric>
#include <sstream>
#include <utility>

void add_error(
    ValidationReport& report, bool& category, std::string message) {
    category = false;
    report.errors.push_back(std::move(message));
}

std::string at_cluster(std::size_t index, const char* message) {
    std::ostringstream output;
    output << "cluster " << index << ": " << message;
    return output.str();
}

std::string at_split(std::size_t index, const char* message) {
    std::ostringstream output;
    output << "split " << index << ": " << message;
    return output.str();
}

bool close(double left, double right, double tolerance) {
    if (left == right) {
        return true;
    }
    if (!std::isfinite(left) || !std::isfinite(right)) {
        return false;
    }
    const double scale = std::max({1.0, std::abs(left), std::abs(right)});
    return std::abs(left - right) <= tolerance * scale;
}



ValidationReport validate_decomposition(
    const Graph& graph,
    const DecompositionConfig& config,
    const DecompositionResult& result,
    std::span<const int> initial_cluster,
    ValidationOptions options) {
    ValidationReport report;
    if (!(config.alpha > 0.0) || !std::isfinite(config.alpha) ||
        !(config.phi > 0.0) || !std::isfinite(config.phi) ||
        !(config.gamma_cmp >= 1.0) || !std::isfinite(config.gamma_cmp) ||
        !(config.property3_constant > 0.0) ||
        !std::isfinite(config.property3_constant) ||
        !(options.tolerance >= 0.0) || !std::isfinite(options.tolerance)) {
        add_error(report, report.parameters_valid,
            "invalid decomposition or validation parameters");
        return report;
    }

    std::vector<int> root;
    if (initial_cluster.empty()) {
        root.resize(graph.node_count());
        std::iota(root.begin(), root.end(), int{0});
    } else {
        root.assign(initial_cluster.begin(), initial_cluster.end());
    }

    std::vector<bool> expected(graph.node_count(), false);
    bool root_valid = true;
    for (int vertex : root) {
        if (vertex >= graph.node_count()) {
            add_error(report, report.partition_valid,
                "initial cluster contains an out-of-range vertex");
            root_valid = false;
            continue;
        }
        if (expected[vertex]) {
            add_error(report, report.partition_valid,
                "initial cluster contains a duplicate vertex");
            root_valid = false;
            continue;
        }
        expected[vertex] = true;
    }
    if (!root_valid) {
        return report;
    }

    ClusterStatistics root_stats;
    try {
        root_stats = cluster_statistics(graph, root);
    } catch (const std::exception& error) {
        add_error(report, report.cluster_metadata_valid,
            std::string("could not recompute initial-cluster statistics: ") +
                error.what());
        return report;
    }
    const double log_m = safe_log2(root_stats.volume);
    if (config.enforce_theorem_parameter_range) {
        const double alpha_limit =
            1.0 / (4.0 * config.gamma_cmp * log_m * log_m);
        if (config.alpha > alpha_limit + options.tolerance) {
            add_error(report, report.parameters_valid,
                "alpha violates the theorem parameter range");
        }
    }

    if (!result.output_expansion_certified) {
        add_error(report, report.certificate_metadata_valid,
            "result does not claim certified output expansion");
    }
    if (result.input_boundary != root_stats.boundary) {
        add_error(report, report.boundary_accounting_valid,
            "recorded input boundary differs from the graph");
    }
    if (result.input_volume != root_stats.volume) {
        add_error(report, report.cluster_metadata_valid,
            "recorded input volume differs from the graph");
    }

    std::vector<bool> seen(graph.node_count(), false);
    Capacity recomputed_output_boundary = 0;
    const ExactConductanceOracle exact(options.exact_max_vertices);
    for (std::size_t index = 0; index < result.clusters.size(); ++index) {
        const OutputCluster& cluster = result.clusters[index];
        bool cluster_vertices_valid = true;
        if (cluster.vertices.empty() && !root.empty()) {
            add_error(report, report.partition_valid,
                at_cluster(index, "is empty"));
            cluster_vertices_valid = false;
        }
        for (int vertex : cluster.vertices) {
            if (vertex >= graph.node_count()) {
                add_error(report, report.partition_valid,
                    at_cluster(index, "contains an out-of-range vertex"));
                cluster_vertices_valid = false;
                continue;
            }
            if (!expected[vertex]) {
                add_error(report, report.partition_valid,
                    at_cluster(index, "contains a vertex outside the initial cluster"));
                cluster_vertices_valid = false;
            }
            if (seen[vertex]) {
                add_error(report, report.partition_valid,
                    at_cluster(index, "overlaps an earlier output cluster"));
                cluster_vertices_valid = false;
            } else {
                seen[vertex] = true;
            }
        }

        if (!(cluster.phi > 0.0) || !std::isfinite(cluster.phi) ||
            cluster.phi + options.tolerance < config.phi) {
            add_error(report, report.cluster_metadata_valid,
                at_cluster(index, "has an invalid phi"));
        }
        const double expected_weight = config.alpha / cluster.phi;
        if (!close(cluster.augmentation_weight, expected_weight,
                options.tolerance)) {
            add_error(report, report.cluster_metadata_valid,
                at_cluster(index, "has the wrong augmentation weight"));
        }
        if (std::isnan(cluster.certified_conductance_lower_bound) ||
            cluster.certified_conductance_lower_bound + options.tolerance <
                cluster.phi) {
            add_error(report, report.certificate_metadata_valid,
                at_cluster(index, "recorded expansion certificate is below phi"));
        }

        if (!cluster_vertices_valid) {
            report.exact_expansion_complete = false;
            ++report.exact_clusters_skipped;
            continue;
        }

        try {
            const ClusterStatistics stats =
                cluster_statistics(graph, cluster.vertices);
            if (stats.volume != cluster.volume) {
                add_error(report, report.cluster_metadata_valid,
                    at_cluster(index, "recorded volume differs from the graph"));
            }
            if (stats.boundary != cluster.boundary) {
                add_error(report, report.cluster_metadata_valid,
                    at_cluster(index, "recorded boundary differs from the graph"));
            }
            auto check = _checked_add(recomputed_output_boundary, stats.boundary);
            if (!check) {
                add_error(report, report.boundary_accounting_valid,
                    "recomputed output-boundary sum overflows Capacity");
            }
            const long double property3_bound =
                static_cast<long double>(config.property3_constant) *
                static_cast<long double>(config.gamma_cmp) *
                std::pow(static_cast<long double>(log_m), 4.0L) *
                static_cast<long double>(cluster.phi) *
                static_cast<long double>(stats.volume);
            if (static_cast<long double>(stats.boundary) >
                property3_bound + options.tolerance) {
                add_error(report, report.property3_valid,
                    at_cluster(index, "violates Property 3"));
            }
        } catch (const std::exception& error) {
            add_error(report, report.cluster_metadata_valid,
                at_cluster(index, error.what()));
        }

        if (options.exact_max_vertices != 0 &&
            cluster.vertices.size() <= options.exact_max_vertices) {
            try {
                const OracleResult check = exact.analyze(
                    graph,
                    cluster.vertices,
                    cluster.augmentation_weight,
                    cluster.phi);
                ++report.exact_clusters_checked;
                if (check.kind != OracleResult::Kind::ExpanderCertificate) {
                    add_error(report, report.exact_expansion_valid,
                        at_cluster(index,
                            "fails independent exact augmented-conductance verification"));
                }
            } catch (const std::exception& error) {
                add_error(report, report.exact_expansion_valid,
                    at_cluster(index, error.what()));
            }
        } else {
            report.exact_expansion_complete = false;
            ++report.exact_clusters_skipped;
        }
    }

    for (int vertex = 0; vertex < graph.node_count(); ++vertex) {
        if (expected[vertex] && !seen[vertex]) {
            add_error(report, report.partition_valid,
                "an initial-cluster vertex is missing from the output partition");
            break;
        }
    }

    Capacity recomputed_split_capacity = 0;
    std::size_t previous_round = 0;
    for (std::size_t index = 0; index < result.splits.size(); ++index) {
        const SplitAudit& split = result.splits[index];
        if (split.round == 0 || split.round < previous_round) {
            add_error(report, report.split_audit_valid,
                at_split(index, "has an invalid round number"));
        }
        previous_round = split.round;
        if (!(split.target_phi > 0.0) ||
            split.target_phi + options.tolerance < config.phi ||
            !std::isfinite(split.target_phi)) {
            add_error(report, report.split_audit_valid,
                at_split(index, "has an invalid target phi"));
        } else if (!close(
                split.augmentation_weight,
                config.alpha / split.target_phi,
                options.tolerance)) {
            add_error(report, report.split_audit_valid,
                at_split(index, "has the wrong augmentation weight"));
        }
        if (split.side_size == 0 || split.side_size >= split.parent_size) {
            add_error(report, report.split_audit_valid,
                at_split(index, "records a trivial cut"));
        }
        if (split.parent_vertices.size() != split.parent_size ||
            split.side_vertices.size() != split.side_size) {
            add_error(report, report.proof_model_valid,
                at_split(index,
                    "does not retain the vertices required for independent verification"));
        } else {
            try {
                const ConductanceResult recomputed = augmented_conductance(
                    graph,
                    split.parent_vertices,
                    split.side_vertices,
                    split.augmentation_weight);
                if (recomputed.cut_capacity != split.cut_capacity) {
                    add_error(report, report.proof_model_valid,
                        at_split(index,
                            "cut capacity differs from exact fractional recomputation"));
                }
                if (!close(
                        recomputed.conductance,
                        split.cut_conductance,
                        options.tolerance)) {
                    add_error(report, report.proof_model_valid,
                        at_split(index,
                            "conductance differs from exact fractional recomputation"));
                }
            } catch (const std::exception& error) {
                const std::string message =
                    std::string("fractional cut recomputation failed: ") +
                    error.what();
                add_error(report, report.proof_model_valid,
                    at_split(index, message.c_str()));
            }
        }
        if (!std::isfinite(split.cut_conductance) ||
            split.cut_conductance < 0.0 ||
            split.cut_conductance >
                config.gamma_cmp * split.target_phi + options.tolerance) {
            add_error(report, report.split_audit_valid,
                at_split(index, "exceeds gamma_cmp times target phi"));
        }
        auto check = _checked_add(recomputed_split_capacity, split.cut_capacity);
        if (!check) {
            add_error(report, report.boundary_accounting_valid,
                "recomputed split-capacity sum overflows Capacity");
        }
    }

    if (recomputed_split_capacity != result.total_split_capacity) {
        add_error(report, report.boundary_accounting_valid,
            "recorded total split capacity is inconsistent");
    }
    if (recomputed_output_boundary != result.output_boundary_sum) {
        add_error(report, report.boundary_accounting_valid,
            "recorded output-boundary sum is inconsistent");
    }
    Capacity expected_output_boundary = root_stats.boundary;
    if (recomputed_split_capacity >
        (std::numeric_limits<Capacity>::max() - expected_output_boundary) / 2) {
        add_error(report, report.boundary_accounting_valid,
            "boundary-accounting identity overflows Capacity");
    } else {
        expected_output_boundary += 2 * recomputed_split_capacity;
        if (recomputed_output_boundary != expected_output_boundary) {
            add_error(report, report.boundary_accounting_valid,
                "boundary-accounting identity does not hold");
        }
    }
    if (!result.boundary_accounting_identity_verified) {
        add_error(report, report.boundary_accounting_valid,
            "result does not claim verified boundary accounting");
    }

    if (options.require_exact_for_all && !report.exact_expansion_complete) {
        add_error(report, report.exact_expansion_valid,
            "complete exact expansion verification was required but skipped");
    }
    return report;
}
