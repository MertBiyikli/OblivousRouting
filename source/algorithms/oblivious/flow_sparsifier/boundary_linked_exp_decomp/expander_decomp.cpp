//
// Created by Mert Biyikli on 21.07.26.
//

#include "algorithms/oblivious/flow_sparsifier/boundary_linked_exp_decomp/expander_decomp.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include "algorithms/oblivious/flow_sparsifier/boundary_linked_exp_decomp/augmented_cluster_view.h"
#include "core/errors.h"

namespace {

constexpr Capacity kMaximumDecimalScale = 1'000'000'000ULL;

struct InputEdge {
    int u{};
    int v{};
    double capacity{};
};

void validate_input_capacity(double value) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument(
            "boundary-linked decomposition requires finite positive edge capacities");
    }
}

bool close_to_integer(long double value) {
    const long double nearest = std::round(value);
    const long double tolerance =
        32.0L * static_cast<long double>(
            std::numeric_limits<double>::epsilon()) *
        std::max(1.0L, std::fabs(value));
    return std::fabs(value - nearest) <= tolerance;
}

Capacity choose_decimal_scale(std::span<const InputEdge> edges) {
    Capacity scale = 1;
    while (true) {
        bool represents_all = true;
        for (const InputEdge& edge : edges) {
            const long double scaled =
                static_cast<long double>(edge.capacity) *
                static_cast<long double>(scale);
            if (!close_to_integer(scaled)) {
                represents_all = false;
                break;
            }
        }
        if (represents_all) {
            return scale;
        }
        if (scale > kMaximumDecimalScale / 10) {
            throw std::invalid_argument(
                "edge capacities have no common decimal unit with at most "
                "nine fractional digits; use the weighted backend");
        }
        scale *= 10;
    }
}

Capacity scaled_capacity(double value, Capacity scale) {
    const long double scaled =
        static_cast<long double>(value) *
        static_cast<long double>(scale);
    if (!close_to_integer(scaled)) {
        throw std::logic_error(
            "capacity ceased to be integral under the selected global scale");
    }
    const long double nearest = std::round(scaled);
    if (nearest <= 0.0L ||
        nearest >
        static_cast<long double>(std::numeric_limits<Capacity>::max())) {
        throw std::overflow_error(
            "scaled edge capacity exceeds boundary-linked Capacity range");
    }
    return static_cast<Capacity>(nearest);
}

}  // namespace


Graph::Graph(const IGraph& igraph)
    : node_count_(igraph.getNumNodes()),
      degrees_(static_cast<std::size_t>(igraph.getNumNodes()), Capacity{0}) {
    std::vector<InputEdge> input_edges;
    input_edges.reserve(igraph.getNumUndirectedEdges());

    const int directed_edge_count = igraph.getNumDirectedEdges();
    for (int e = 0; e < directed_edge_count; ++e) {
        const auto [head, tail] = igraph.getEdgeEndpoints(e);

        // IGraph stores both orientations of each undirected edge.
        if (head >= tail) {
            continue;
        }

        if (head < 0 || tail < 0 ||
            head >= node_count_ || tail >= node_count_) {
            throw std::logic_error(
                "IGraph returned an invalid edge endpoint");
        }

        const int anti_edge = igraph.getAntiEdge(e);
        if (anti_edge == INVALID_EDGE_ID) {
            throw std::logic_error(
                "IGraph directed edge has no reverse edge");
        }
        const double forward_capacity = igraph.getEdgeCapacity(e);
        const double reverse_capacity = igraph.getEdgeCapacity(anti_edge);
        validate_input_capacity(forward_capacity);
        validate_input_capacity(reverse_capacity);
        if (forward_capacity != reverse_capacity) {
            throw std::invalid_argument(
                "boundary-linked decomposition requires symmetric capacities");
        }
        input_edges.push_back({head, tail, forward_capacity});
    }

    if (input_edges.size() !=
        static_cast<std::size_t>(igraph.getNumUndirectedEdges())) {
        throw std::logic_error(
            "IGraph directed-edge representation is not symmetric");
    }

    if (input_edges.empty()) {
        return;
    }

    const Capacity decimal_scale = choose_decimal_scale(input_edges);
    std::vector<Capacity> scaled;
    scaled.reserve(input_edges.size());
    Capacity common_divisor = 0;
    for (const InputEdge& edge : input_edges) {
        const Capacity value =
            scaled_capacity(edge.capacity, decimal_scale);
        scaled.push_back(value);
        common_divisor = std::gcd(common_divisor, value);
    }
    if (common_divisor == 0) {
        throw std::logic_error("positive capacities have zero global GCD");
    }
    capacity_unit_ =
        static_cast<long double>(common_divisor) /
        static_cast<long double>(decimal_scale);

    edges_.reserve(input_edges.size());
    for (std::size_t index = 0; index < input_edges.size(); ++index) {
        const InputEdge& input = input_edges[index];
        const Capacity multiplicity = scaled[index] / common_divisor;
        edges_.emplace_back(input.u, input.v, multiplicity);
        degrees_[input.u] = checked_add(
            degrees_[input.u], multiplicity, "graph degree overflow");
        degrees_[input.v] = checked_add(
            degrees_[input.v], multiplicity, "graph degree overflow");
    }
}

Graph::Graph(int node_count, std::vector<Edge> edges)
    : node_count_(node_count),
      edges_(std::move(edges)),
      degrees_(static_cast<std::size_t>(node_count), Capacity{0}) {
    for (const Edge& edge : edges_) {
        if (edge.u < 0 || edge.v < 0 ||
            edge.u >= node_count_ || edge.v >= node_count_) {
            throw std::invalid_argument(
                "edge endpoint is outside the graph");
        }
        if (edge.multiplicity == 0) {
            throw std::invalid_argument(
                "zero edge multiplicity is not allowed");
        }

        degrees_[edge.u] = checked_add(degrees_[edge.u],edge.multiplicity,"graph degree overflow");

        if (edge.u != edge.v) {
            degrees_[edge.v] = checked_add(degrees_[edge.v],edge.multiplicity,"graph degree overflow");
        }
    }
}

Capacity Graph::degree(int v) const {
    if (v >= node_count_) {
        throw std::out_of_range("node is outside the graph");
    }
    return degrees_[v];
}

Capacity Graph::volume(std::span<const int> vertices) const {
    Capacity result = 0;
    for (int v : vertices) {
        result = checked_add(
            result,
            degree(v),
            "graph volume overflow");
    }
    return result;
}

struct State {
    std::vector<int> vertices;
    bool active{true};
    bool expands{false};
    double phi{};
    double conductance_lower_bound{};
    ExpansionProofAudit proof;
};

ExpansionProofAudit exact_proof() {
    ExpansionProofAudit proof;
    proof.basis = ExpansionProofBasis::ExactEnumeration;
    proof.deterministic = true;
    proof.paper_preconditions_verified = true;
    proof.implementation_constants_locked = true;
    return proof;
}


std::vector<int> complement(
    std::span<const int> cluster, std::span<const int> side) {
    std::unordered_set<int> selected(side.begin(), side.end());
    std::vector<int> result;
    result.reserve(cluster.size() - side.size());
    for (int v : cluster) {
        if (!selected.contains(v)) {
            result.push_back(v);
        }
    }
    return result;
}

void validate_partition(
    std::span<const int> parent,
    std::span<const int> side,
    std::span<const int> other) {
    if (side.empty() || other.empty() || side.size() + other.size() != parent.size()) {
        throw std::logic_error("oracle returned a trivial cut");
    }
    std::unordered_set<int> union_set;
    for (int v : side) {
        union_set.insert(v);
    }
    for (int v : other) {
        if (!union_set.insert(v).second) {
            throw std::logic_error("oracle cut sides overlap");
        }
    }
    for (int v : parent) {
        if (!union_set.contains(v)) {
            throw std::logic_error("oracle cut does not partition its cluster");
        }
    }
}

ConductanceResult augmented_conductance(
    const Graph& graph,
    std::span<const int> cluster,
    std::span<const int> side,
    double augmentation_weight) {
    const ProofAugmentedCut cut = AugmentedClusterView(
        graph, cluster, augmentation_weight).evaluate_proof_cut(side);
    return {
        cut.conductance,
        cut.cut_capacity,
        cut.side_volume,
        cut.complement_volume};
}

OracleResult ExactConductanceOracle::analyze(
    const Graph& graph,
    std::span<const int> cluster,
    double augmentation_weight,
    double target_phi) const {
    if (cluster.size() > max_vertices_) {
        throw std::runtime_error(
            "exact conductance oracle limit exceeded; install a certified "
            "cut-matching/trimming oracle for larger clusters");
    }
    if (cluster.size() <= 1) {
        return {
            OracleResult::Kind::ExpanderCertificate, {},
            std::numeric_limits<double>::infinity(), exact_proof()};
    }
    if (cluster.size() >= 64) {
        throw std::runtime_error("exact oracle bit-mask representation exceeded");
    }

    const std::uint64_t subset_count = std::uint64_t{1} << cluster.size();
    const AugmentedClusterView view(graph, cluster, augmentation_weight);
    double minimum = std::numeric_limits<double>::infinity();
    std::vector<int> minimizer;

    // Fix vertex 0 in the side, thereby enumerating each unordered cut once.
    for (std::uint64_t mask = 1; mask + 1 < subset_count; ++mask) {
        if ((mask & 1U) == 0) {
            continue;
        }
        std::vector<int> side;
        side.reserve(std::popcount(mask));
        for (std::size_t i = 0; i < cluster.size(); ++i) {
            if ((mask & (std::uint64_t{1} << i)) != 0) {
                side.push_back(cluster[i]);
            }
        }
        const ProofAugmentedCut value = view.evaluate_proof_cut(side);
        if (value.conductance < minimum) {
            minimum = value.conductance;
            minimizer = std::move(side);
        }
    }

    if (minimum + kTolerance < target_phi) {
        return {
            OracleResult::Kind::SparseCut, std::move(minimizer), minimum, {}};
    }
    return {
        OracleResult::Kind::ExpanderCertificate, {}, minimum, exact_proof()};
}

std::string ExactConductanceOracle::name() const {
    return "exact-exponential-conductance";
}

BoundaryLinkedDecomposer::BoundaryLinkedDecomposer(
    DecompositionConfig config,
    std::shared_ptr<const CertifiedSparseCutOracle> oracle)
    : config_(config), oracle_(std::move(oracle)) {
    if (!oracle_) {
        throw std::invalid_argument("a certified sparse-cut oracle is required");
    }
    if (!(config_.alpha > 0.0) || !(config_.phi > 0.0) ||
        !(config_.gamma_cmp >= 1.0) || !(config_.property3_constant > 0.0)) {
        throw std::invalid_argument("invalid decomposition parameters");
    }
}

DecompositionResult BoundaryLinkedDecomposer::decompose(const Graph& graph, std::span<const int> initial_cluster) const {
    std::vector<int> root;
    if (initial_cluster.empty()) {
        root.resize(graph.node_count());
        std::iota(root.begin(), root.end(), int{0});
    } else {
        root.assign(initial_cluster.begin(), initial_cluster.end());
    }
    // Also validates uniqueness and membership.
    const ClusterStatistics root_stats = cluster_statistics(graph, root);
    const double log_m = safe_log2(root_stats.volume);
    const double alpha_limit = 1.0 /
        (4.0 * config_.gamma_cmp * log_m * log_m);
    if (config_.enforce_theorem_parameter_range &&
        config_.alpha > alpha_limit + kTolerance) {
        throw std::invalid_argument(
            "alpha violates alpha <= 1/(4 gamma_cmp log_2^2(m))");
    }

    std::vector<State> states;
    states.push_back(State{std::move(root), true, false, 0.0, 0.0, {}});
    DecompositionResult result;
    result.oracle_name = oracle_->name();
    result.goranci_running_time_certified = oracle_->has_goranci_running_time();
    result.input_boundary = root_stats.boundary;
    result.input_volume = root_stats.volume;

    std::size_t round = 0;
    while (std::any_of(states.begin(), states.end(), [](const State& state) {
        return state.active;
    })) {
        ++round;
        long double active_boundary = 0.0L;
        long double active_volume = 0.0L;
        for (State& state : states) {
            if (!state.active) {
                continue;
            }
            state.expands = false;
            const ClusterStatistics stats = cluster_statistics(graph, state.vertices);
            active_boundary += stats.boundary;
            active_volume += stats.volume;
        }
        const double adaptive = active_volume == 0.0L ? config_.phi :
            static_cast<double>(active_boundary /
                (8.0L * config_.gamma_cmp * log_m * log_m * active_volume));
        const double target_phi = std::max(config_.phi, adaptive);
        const double augmentation_weight = config_.alpha / target_phi;

        for (;;) {
            const auto it = std::find_if(states.begin(), states.end(), [](const State& state) {
                return state.active && !state.expands;
            });
            if (it == states.end()) {
                break;
            }
            const std::size_t index = static_cast<std::size_t>(it - states.begin());
            const OracleResult oracle_result = oracle_->analyze(
                graph, states[index].vertices, augmentation_weight, target_phi);
            if (oracle_result.kind == OracleResult::Kind::ExpanderCertificate) {
                if (oracle_result.conductance + kTolerance < target_phi) {
                    throw std::logic_error(
                        "oracle expansion lower bound is below requested phi");
                }
                states[index].expands = true;
                states[index].conductance_lower_bound = oracle_result.conductance;
                states[index].proof = oracle_result.proof;
                continue;
            }

            std::vector<int> other = complement(states[index].vertices, oracle_result.side);
            validate_partition(states[index].vertices, oracle_result.side, other);
            const ConductanceResult verified = augmented_conductance(
                graph, states[index].vertices, oracle_result.side, augmentation_weight);
            if (verified.conductance >
                config_.gamma_cmp * target_phi + 1e-12) {
                            throw std::runtime_error(
                                "practical balanced cut failed augmented conductance verification: "
                                "observed=" + std::to_string(verified.conductance) +
                                ", permitted=" +
                                std::to_string(config_.gamma_cmp * target_phi) +
                                ", target_phi=" + std::to_string(target_phi) +
                                ", gamma_cmp=" + std::to_string(config_.gamma_cmp));
                }
            result.splits.push_back({
                round,
                target_phi,
                augmentation_weight,
                verified.conductance,
                verified.cut_capacity,
                states[index].vertices.size(),
                oracle_result.side.size(),
                oracle_result.kind == OracleResult::Kind::CertifiedSide,
                states[index].vertices,
                oracle_result.side});
            if (oracle_result.kind == OracleResult::Kind::CertifiedSide) {
                if (oracle_result.conductance + kTolerance < target_phi) {
                    throw std::logic_error(
                        "oracle side-certificate lower bound is below requested phi");
                }
                states[index] = State{
                    oracle_result.side, true, true, 0.0,
                    oracle_result.conductance, oracle_result.proof};
            } else {
                states[index] = State{
                    oracle_result.side, true, false, 0.0, 0.0, {}};
            }
            states.push_back(State{
                std::move(other), true, false, 0.0, 0.0, {}});
        }

        bool deactivated = false;
        for (State& state : states) {
            if (!state.active) {
                continue;
            }
            const ClusterStatistics stats = cluster_statistics(graph, state.vertices);
            const long double threshold = config_.property3_constant *
                config_.gamma_cmp * std::pow(log_m, 4.0) * target_phi * stats.volume;
            if (static_cast<long double>(stats.boundary) <= threshold + kTolerance) {
                state.active = false;
                state.phi = target_phi;
                deactivated = true;
            }
        }
        if (!deactivated) {
            throw std::logic_error(
                "adaptive round made no progress; check constants and graph arithmetic");
        }
    }

    result.output_expansion_certified = true;
    for (const State& state : states) {
        const ClusterStatistics stats = cluster_statistics(graph, state.vertices);
        if (state.conductance_lower_bound + kTolerance < state.phi) {
            throw std::logic_error("final boundary-linked expansion certificate failed");
        }
        result.output_boundary_sum += stats.boundary;
        result.clusters.push_back({
            state.vertices,
            state.phi,
            config_.alpha / state.phi,
            stats.boundary,
            stats.volume,
            state.conductance_lower_bound,
            state.proof});
        result.all_expansion_proofs_deterministic =
            result.all_expansion_proofs_deterministic &&
            state.proof.deterministic;
        result.all_expansion_proofs_have_locked_constants =
            result.all_expansion_proofs_have_locked_constants &&
            state.proof.paper_preconditions_verified &&
            state.proof.implementation_constants_locked;
        if (state.proof.basis ==
                ExpansionProofBasis::RandomizedCutMatching ||
            state.proof.basis ==
                ExpansionProofBasis::RandomizedCutMatchingAndTrimming) {
            ++result.randomized_expansion_proof_count;
        }
    }
    for (const SplitAudit& split : result.splits) {
        if (split.cut_capacity >
            std::numeric_limits<Capacity>::max() - result.total_split_capacity) {
            throw std::overflow_error("total split capacity overflow");
        }
        result.total_split_capacity += split.cut_capacity;
    }
    if (result.total_split_capacity >
        (std::numeric_limits<Capacity>::max() - result.input_boundary) / 2) {
        throw std::overflow_error("boundary accounting identity overflow");
    }
    const Capacity expected_output_boundary =
        result.input_boundary + 2 * result.total_split_capacity;
    result.boundary_accounting_identity_verified =
        expected_output_boundary == result.output_boundary_sum;
    if (!result.boundary_accounting_identity_verified) {
        throw std::logic_error(
            "boundary accounting identity failed: input=" +
            std::to_string(result.input_boundary) +
            ", split_sum=" +
            std::to_string(result.total_split_capacity) +
            ", expected_output=" +
            std::to_string(expected_output_boundary) +
            ", actual_output=" +
            std::to_string(result.output_boundary_sum));
    }
    // Claim 4.12 in Goranci et al.: with Z=log_2(vol(U)), the explicit
    // charging bound is 4b + 8*gamma_cmp*Z*phi*m.
    result.property1_upper_bound =
        4.0L * static_cast<long double>(result.input_boundary) +
        8.0L * static_cast<long double>(config_.gamma_cmp) *
            static_cast<long double>(log_m) *
            static_cast<long double>(config_.phi) *
            static_cast<long double>(result.input_volume);
    result.property1_bound_verified =
        static_cast<long double>(result.output_boundary_sum) <=
        result.property1_upper_bound + kTolerance;
    if (!result.property1_bound_verified) {
        throw std::logic_error("Property 1 explicit boundary bound failed");
    }
    return result;
}
