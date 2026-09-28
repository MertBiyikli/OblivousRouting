//
// Created by Mert Biyikli on 21.07.26.
//

#ifndef OBLIVIOUSROUTING_EXPANDER_DECOMP_H
#define OBLIVIOUSROUTING_EXPANDER_DECOMP_H
#include <cstdint>
#include <cmath>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include "data_structures/graph/Igraph.h"
#include "core/errors.h"

using Capacity = std::uint64_t;

constexpr double kTolerance = 1e-12;

struct Edge {
    int u{};
    int v{};
    Capacity multiplicity{1};
};

// An undirected rational-capacity graph represented as a normalized
// multigraph. IGraph inputs are multiplied by one global decimal scale and
// divided by the global GCD. Thus capacity_unit() original-capacity units
// correspond to one stored parallel edge. A self-loop contributes its
// multiplicity once to the degree, matching the augmented-graph convention in
// Goranci et al.
class Graph {
public:
    Graph(const IGraph& igraph);
    Graph(int int_count, std::vector<Edge> edges);

    int node_count() const noexcept { return node_count_; }
    std::span<const Edge> edges() const noexcept { return edges_; }
    Capacity degree(int v) const;
    Capacity volume(std::span<const int> vertices) const;
    long double capacity_unit() const noexcept { return capacity_unit_; }

private:
    int node_count_{};
    std::vector<Edge> edges_;
    std::vector<Capacity> degrees_;
    long double capacity_unit_{1.0L};
};



struct ClusterStatistics {
    std::vector<Capacity> internal_degree;
    std::vector<Capacity> boundary_degree;
    Capacity volume{};
    Capacity boundary{};
};

// The degree arrays are indexed by the position in `cluster`.

// Number of copies of original edges crossing (side, cluster \ side).
[[nodiscard]] Capacity internal_cut_capacity(
    const Graph& graph,
    std::span<const int> cluster,
    std::span<const int> side);


struct ConductanceResult {
    double conductance{1.0};
    Capacity cut_capacity{};
    long double side_volume{};
    long double complement_volume{};
};

enum class ExpansionProofBasis {
    Unspecified,
    ExactEnumeration,
    DeterministicTrivial,
    RandomizedCutMatching,
    RandomizedCutMatchingAndTrimming,
};

struct ExpansionProofAudit {
    ExpansionProofBasis basis{ExpansionProofBasis::Unspecified};
    bool deterministic{false};
    bool paper_preconditions_verified{false};
    // False means the high-level proof applies only after the remaining
    // hidden-constant correspondence audit of the concrete backend.
    bool implementation_constants_locked{false};
    // `10` records the paper's o(m^-10) per-call failure statement. It is an
    // exponent, not a claimed numeric probability bound.
    double little_o_failure_exponent{};
    std::uint32_t random_seed{};
    int cut_matching_iterations{};
    int cut_matching_rounds{};
    long long per_round_capacity{};
    long long flow_scale{1};
    long long realized_congestion{};
    std::size_t terminals{};
};

// Conductance in G[cluster]^augmentation_weight. Boundary self-loops change
// volumes but never contribute to the cut numerator.
[[nodiscard]] ConductanceResult augmented_conductance(
    const Graph& graph,
    std::span<const int> cluster,
    std::span<const int> side,
    double augmentation_weight);

struct OracleResult {
    // CertifiedSide is Lemma 4.7 case 2: `side` is certified expanding and its
    // nonempty complement remains active in the current Algorithm 1 round.
    enum class Kind { SparseCut, CertifiedSide, ExpanderCertificate };
    Kind kind{Kind::ExpanderCertificate};
    std::vector<int> side;
    // Exact for the reference oracle; otherwise a certified lower bound.
    double conductance{1.0};
    ExpansionProofAudit proof;
};

// Proof-facing boundary: an implementation may return ExpanderCertificate
// only if it certifies every nontrivial cut in the requested augmented graph.
class CertifiedSparseCutOracle {
public:
    virtual ~CertifiedSparseCutOracle() = default;
    [[nodiscard]] virtual OracleResult analyze(
        const Graph& graph,
        std::span<const int> cluster,
        double augmentation_weight,
        double target_phi) const = 0;
    [[nodiscard]] virtual std::string name() const = 0;
    [[nodiscard]] virtual bool has_goranci_running_time() const noexcept = 0;
};

// Exponential reference oracle. It is suitable for correctness tests and small
// clusters, and deliberately throws above max_vertices instead of guessing.
class ExactConductanceOracle final : public CertifiedSparseCutOracle {
public:
    explicit ExactConductanceOracle(std::size_t max_vertices = 24)
        : max_vertices_(max_vertices) {}

    [[nodiscard]] OracleResult analyze(
        const Graph& graph,
        std::span<const int> cluster,
        double augmentation_weight,
        double target_phi) const override;
    [[nodiscard]] std::string name() const override;
    [[nodiscard]] bool has_goranci_running_time() const noexcept override {
        return false;
    }

private:
    std::size_t max_vertices_;
};

struct DecompositionConfig {
    double alpha{};
    double phi{};
    double gamma_cmp{1.0};
    double property3_constant{80.0};
    bool enforce_theorem_parameter_range{true};
};

struct OutputCluster {
    std::vector<int> vertices;
    double phi{};
    double augmentation_weight{};
    Capacity boundary{};
    Capacity volume{};
    double certified_conductance_lower_bound{};
    ExpansionProofAudit proof;
};

struct SplitAudit {
    std::size_t round{};
    double target_phi{};
    double augmentation_weight{};
    double cut_conductance{};
    Capacity cut_capacity{};
    std::size_t parent_size{};
    std::size_t side_size{};
    bool side_was_certified_expanding{false};
    // Retained so the independent validator can recompute the exact
    // fractional augmented cut instead of trusting solver metadata.
    std::vector<int> parent_vertices;
    std::vector<int> side_vertices;
};

struct DecompositionResult {
    std::vector<OutputCluster> clusters;
    std::vector<SplitAudit> splits;
    std::string oracle_name;
    bool output_expansion_certified{false};
    bool goranci_running_time_certified{false};
    bool all_expansion_proofs_deterministic{true};
    bool all_expansion_proofs_have_locked_constants{true};
    std::size_t randomized_expansion_proof_count{};
    Capacity input_boundary{};
    Capacity input_volume{};
    Capacity output_boundary_sum{};
    Capacity total_split_capacity{};
    bool boundary_accounting_identity_verified{false};
    long double property1_upper_bound{};
    bool property1_bound_verified{false};
};

// Implements the adaptive outer loop of Algorithm 1 in Goranci et al. The
// quality of the oracle determines the running-time guarantee, but not whether
// an emitted expansion certificate is trusted.
class BoundaryLinkedDecomposer {
public:
    BoundaryLinkedDecomposer(
        DecompositionConfig config,
        std::shared_ptr<const CertifiedSparseCutOracle> oracle);

    DecompositionResult decompose(
        const Graph& graph,
        std::span<const int> initial_cluster = {}) const;

private:
    DecompositionConfig config_;
    std::shared_ptr<const CertifiedSparseCutOracle> oracle_;
};

inline Capacity checked_add(Capacity left, Capacity right, const char* message) {
    if (right > std::numeric_limits<Capacity>::max() - left) {
        throw std::overflow_error(message);
    }
    return left + right;
}

inline Result<void> _checked_add(Capacity& value, Capacity increment) {
    if (increment > std::numeric_limits<Capacity>::max() - value) {
        return makeErrorMessage(ErrorCode::LogicError ,"edge multiplicity overflow");
    }
    value += increment;
    return {};
}

inline double safe_log2(Capacity volume) {
    return std::max(1.0, std::log2(static_cast<double>(std::max<Capacity>(2, volume))));
}

inline std::unordered_map<int, std::size_t> positions(
    std::span<const int> vertices) {
    std::unordered_map<int, std::size_t> result;
    result.reserve(vertices.size());
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        if (!result.emplace(vertices[i], i).second) {
            throw std::invalid_argument("cluster contains a duplicate vertex");
        }
    }
    return result;
}

inline ClusterStatistics cluster_statistics(
    const Graph& graph,
    std::span<const int> cluster) {
    const auto position = positions(cluster);

    ClusterStatistics result;
    result.internal_degree.resize(cluster.size());
    result.boundary_degree.resize(cluster.size());

    for (const Edge& edge : graph.edges()) {
        const auto u = position.find(edge.u);
        const auto v = position.find(edge.v);

        const bool u_inside = u != position.end();
        const bool v_inside = v != position.end();

        if (u_inside && v_inside) {
            result.internal_degree[u->second] = checked_add(
                result.internal_degree[u->second],
                edge.multiplicity,
                "internal degree overflow");

            if (edge.u != edge.v) {
                result.internal_degree[v->second] = checked_add(
                    result.internal_degree[v->second],
                    edge.multiplicity,
                    "internal degree overflow");
            }
        } else if (u_inside) {
            result.boundary_degree[u->second] = checked_add(
                result.boundary_degree[u->second],
                edge.multiplicity,
                "boundary degree overflow");

            result.boundary = checked_add(
                result.boundary,
                edge.multiplicity,
                "cluster boundary overflow");
        } else if (v_inside) {
            result.boundary_degree[v->second] = checked_add(
                result.boundary_degree[v->second],
                edge.multiplicity,
                "boundary degree overflow");

            result.boundary = checked_add(
                result.boundary,
                edge.multiplicity,
                "cluster boundary overflow");
        }
    }

    result.volume = graph.volume(cluster);
    return result;
}

#endif //OBLIVIOUSROUTING_EXPANDER_DECOMP_H
