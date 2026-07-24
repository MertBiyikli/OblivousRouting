//
// Created by Mert Biyikli on 21.07.26.
//
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "algorithms/oblivious/flow_sparsifier/boundary_linked_exp_decomp/practical_adapter.h"
#include "algorithms/oblivious/flow_sparsifier/boundary_linked_exp_decomp/augmented_cluster_view.h"
#include "cut_matching.hpp"
#include "datastructures/unit_flow.hpp"
#include "trimming.hpp"
#include "util.hpp"

int checked_size(std::size_t value, const char* message) {
    if (value > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error(message);
    }
    return static_cast<int>(value);
}

int checked_capacity(Capacity value, const char* message) {
    if (value > static_cast<Capacity>(std::numeric_limits<int>::max())) {
        throw std::overflow_error(message);
    }
    return static_cast<int>(value);
}

CutMatching::Parameters certified_candidate_parameters() {
    return {
        .tConst = 22,
        .tFactor = 5.0,
        .minIterations = 0,
        .minBalance = 0.0,
        .samplePotential = false,
        .balancedCutStrategy = false,
        .use_cut_heuristics = false,
        .use_balanced_partitions = false,
        .use_potential_based_dynamic_stopping_criterion = false,
        .stop_flow_at_fraction = false,
        .warm_start_unit_flow = false,
        .trim_with_max_flow_first = false,
        .krv_step_first = true,  // changed
        .kahan_error = false,
        .num_flow_vectors = 1,
        .tune_num_flow_vectors = false,
        .break_at_empty_terminals = false,
    };
}

int cut_matching_rounds(std::size_t terminals) {
    if (terminals == 0) {
        return 1;
    }
    const double logarithm = std::log10(
        static_cast<double>(std::max<std::size_t>(2, terminals)));
    return std::max(1, 22 + static_cast<int>(std::ceil(5.0 * logarithm * logarithm)));
}

ExpansionProofAudit deterministic_trivial_proof(std::size_t terminals) {
    ExpansionProofAudit proof;
    proof.basis = ExpansionProofBasis::DeterministicTrivial;
    proof.deterministic = true;
    proof.paper_preconditions_verified = true;
    proof.implementation_constants_locked = true;
    proof.terminals = terminals;
    return proof;
}

ExpansionProofAudit randomized_practical_proof(
    ExpansionProofBasis basis,
    const CutMatching::Result& practical,
    std::uint32_t seed,
    std::size_t terminals) {
    ExpansionProofAudit proof;
    proof.basis = basis;
    proof.deterministic = false;
    proof.paper_preconditions_verified = true;
    // Fixed-point scaling represents a capacity ratio no larger than
    // 1/(phi*T). The concrete T constants remain theorem-unlocked.
    proof.implementation_constants_locked = false;
    proof.little_o_failure_exponent = 10.0;
    proof.random_seed = seed;
    proof.cut_matching_iterations = practical.iterations;
    proof.cut_matching_rounds = practical.roundsPlanned;
    proof.per_round_capacity = practical.perRoundCapacity;
    proof.flow_scale = practical.flowScale;
    proof.realized_congestion = practical.congestion;
    proof.terminals = terminals;
    return proof;
}

std::optional<std::vector<int>> disconnected_positive_component(
    const AugmentedClusterView& view,
    const ExplicitAugmentedInstance& instance) {
    std::vector<std::size_t> parent(instance.original_vertices.size());
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    const auto find = [&](std::size_t start, auto&& self) -> std::size_t {
        if (parent[start] != start) {
            parent[start] = self(parent[start], self);
        }
        return parent[start];
    };
    for (const LocalRoutingEdge& edge : instance.routing_edges) {
        const std::size_t u = find(edge.u, find);
        const std::size_t v = find(edge.v, find);
        if (u != v) {
            parent[v] = u;
        }
    }

    std::vector<Capacity> component_volume(parent.size());
    for (std::size_t v = 0; v < parent.size(); ++v) {
        const std::size_t root = find(v, find);
        component_volume[root] += view.augmented_degree(
            instance.original_vertices[v]);
    }
    for (std::size_t root = 0; root < parent.size(); ++root) {
        if (component_volume[root] == 0 ||
            component_volume[root] == view.augmented_volume()) {
            continue;
        }
        std::vector<int> side;
        for (std::size_t v = 0; v < parent.size(); ++v) {
            if (find(v, find) == root) {
                side.push_back(instance.original_vertices[v]);
            }
        }
        return side;
    }
    return std::nullopt;
}


struct PracticalFlowInstance::Impl {
    std::unique_ptr<UnitFlow::Graph> flow;
    std::unique_ptr<UnitFlow::Graph> subdivision;
    std::vector<int> subdivision_index;
    std::vector<int> original_vertices;
    PracticalInstanceSummary summary;
};

PracticalFlowInstance::PracticalFlowInstance(std::unique_ptr<Impl> impl)
    : impl_(std::move(impl)) {}

PracticalFlowInstance::PracticalFlowInstance(PracticalFlowInstance&&) noexcept = default;
PracticalFlowInstance& PracticalFlowInstance::operator=(PracticalFlowInstance&&) noexcept = default;
PracticalFlowInstance::~PracticalFlowInstance() = default;

PracticalFlowInstance PracticalFlowInstance::build(
    const AugmentedClusterView& view,
    std::size_t maximum_explicit_terminals) {
    const ExplicitAugmentedInstance source =
        view.build_explicit_instance(maximum_explicit_terminals);
    const int node_count = checked_size(
        source.original_vertices.size(), "practical vertex limit exceeded");
    const int terminal_count = checked_size(
        source.terminal_count, "practical terminal limit exceeded");
    if (terminal_count > std::numeric_limits<int>::max() - node_count) {
        throw std::overflow_error("practical subdivision vertex limit exceeded");
    }
    (void)checked_capacity(
        view.augmented_volume(),
        "practical augmented volume exceeds int range");

    std::vector<UnitFlow::Edge> routing_edges;
    routing_edges.reserve(source.routing_edges.size());
    for (const LocalRoutingEdge& edge : source.routing_edges) {
        routing_edges.emplace_back(
            checked_size(edge.u, "practical routing endpoint overflow"),
            checked_size(edge.v, "practical routing endpoint overflow"), 0);
    }

    std::vector<int> augmented_degrees;
    augmented_degrees.reserve(source.augmented_degrees.size());
    for (Capacity degree : source.augmented_degrees) {
        augmented_degrees.push_back(checked_capacity(
            degree, "practical augmented degree exceeds int range"));
    }

    std::vector<UnitFlow::Edge> subdivision_edges;
    subdivision_edges.reserve(source.terminal_incidences.size());
    for (const TerminalIncidence& incidence : source.terminal_incidences) {
        subdivision_edges.emplace_back(
            checked_size(incidence.vertex,
                "practical subdivision endpoint overflow"),
            node_count + checked_size(incidence.terminal,
                "practical terminal id overflow"),
            0);
    }

    auto impl = std::make_unique<Impl>();
    impl->flow = std::make_unique<UnitFlow::Graph>(node_count, routing_edges, std::move(augmented_degrees));
    impl->subdivision = std::make_unique<UnitFlow::Graph>(
        node_count + terminal_count, subdivision_edges);
    impl->subdivision_index.assign(
        static_cast<std::size_t>(node_count + terminal_count), -1);
    for (int terminal = 0; terminal < terminal_count; ++terminal) {
        impl->subdivision_index[static_cast<std::size_t>(node_count + terminal)] = 0;
    }
    impl->original_vertices = source.original_vertices;
    impl->summary = {
        source.original_vertices.size(),
        source.routing_edges.size(),
        source.terminal_count,
        static_cast<std::size_t>(node_count + terminal_count),
        source.terminal_incidences.size(),
        view.augmented_volume()};

    for (int v = 0; v < node_count; ++v) {
        const int expected = checked_capacity(
            view.augmented_degree(source.original_vertices[static_cast<std::size_t>(v)]),
            "practical degree validation overflow");
        if (impl->flow->globalDegree(v) != expected ||
            impl->subdivision->degree(v) != expected) {
            throw std::logic_error(
                "practical augmented-degree construction invariant failed");
        }
    }
    return PracticalFlowInstance(std::move(impl));
}

const PracticalInstanceSummary& PracticalFlowInstance::summary() const noexcept {
    return impl_->summary;
}

UnitFlow::Graph& PracticalFlowInstance::flow_graph() noexcept {
    return *impl_->flow;
}

UnitFlow::Graph& PracticalFlowInstance::subdivision_graph() noexcept {
    return *impl_->subdivision;
}

PracticalCutMatchingOracle::PracticalCutMatchingOracle(
    PracticalOracleConfig config)
    : config_(config), seed_generator_(config.random_seed) {
    if (config_.maximum_explicit_terminals == 0 ||
        (config_.allow_exact_fallback &&
         config_.exact_fallback_max_vertices == 0) ||
        !std::isfinite(config_.gamma_cmp) || config_.gamma_cmp < 1.0) {
        throw std::invalid_argument("invalid practical oracle configuration");
    }
}

OracleResult PracticalCutMatchingOracle::analyze(
    const Graph& graph,
    std::span<const int> cluster,
    double augmentation_weight,
    double target_phi) const {
    if (!std::isfinite(target_phi) || target_phi <= 0.0) {
        throw std::invalid_argument("target phi must be finite and positive");
    }
    if (!(augmentation_weight < 1.0 / (8.0 * target_phi))) {
        throw std::invalid_argument(
            "Goranci Lemma 4.7 requires w < 1/(8 phi)");
    }
    const AugmentedClusterView view(graph, cluster, augmentation_weight);
    const auto fallback_or_throw = [&](std::string message) -> OracleResult {
        if (!config_.allow_exact_fallback ||
            cluster.size() > config_.exact_fallback_max_vertices) {
            throw std::runtime_error(message);
        }
        return ExactConductanceOracle(config_.exact_fallback_max_vertices)
            .analyze(graph, cluster, augmentation_weight, target_phi);
    };

    const ExplicitAugmentedInstance preview =
        view.build_explicit_instance(config_.maximum_explicit_terminals);
    if (cluster.size() <= 1 || view.augmented_volume() == 0) {
        return {
            OracleResult::Kind::ExpanderCertificate, {},
            std::numeric_limits<double>::infinity(),
            deterministic_trivial_proof(preview.terminal_count)};
    }
    if (const auto component = disconnected_positive_component(view, preview)) {
        const ProofAugmentedCut cut = view.evaluate_proof_cut(*component);
        if (cut.cut_capacity != 0 || cut.conductance != 0.0) {
            throw std::logic_error("disconnected-component cut verification failed");
        }
        return {
            OracleResult::Kind::SparseCut,
            std::move(*component),
            0.0,
            {}};
    }
    const bool exact_available =
    cluster.size() <= config_.exact_fallback_max_vertices &&
    cluster.size() < 64;

    const auto run_exact = [&]() -> OracleResult {
        return ExactConductanceOracle(
            config_.exact_fallback_max_vertices)
            .analyze(
                graph,
                cluster,
                augmentation_weight,
                target_phi);
    };

    if (config_.selection_mode == OracleSelectionMode::ExactOnly) {
        if (!exact_available) {
            throw std::runtime_error(
                "ExactOnly selected, but cluster exceeds exact limit");
        }
        return run_exact();
    }

    if (config_.selection_mode == OracleSelectionMode::Auto &&
        exact_available) {
        return run_exact();
        }
    constexpr std::size_t minimum_practical_terminals = 8;

    if (preview.terminal_count < minimum_practical_terminals) {
        if (exact_available) {
            return run_exact();
        }

        throw std::runtime_error(
            "practical cut matching requires at least 8 terminals, "
            "but exact enumeration is unavailable");
    }
    if (preview.routing_edges.size() == 1 && preview.terminal_count == 1) {
        if (target_phi <= 1.0) {
            return {
                OracleResult::Kind::ExpanderCertificate, {}, 1.0,
                deterministic_trivial_proof(preview.terminal_count)};
        }
        return {
            OracleResult::Kind::SparseCut,
            {preview.original_vertices.at(preview.routing_edges.front().u)},
            1.0,
            {}};
    }
    if (preview.routing_edges.empty() || preview.terminal_count <= 1) {
        return fallback_or_throw(
            "degenerate augmented graph is unsupported by practical cut matching");
    }

    PracticalFlowInstance instance = PracticalFlowInstance::build(
        view, config_.maximum_explicit_terminals);
    const std::uint32_t call_seed = seed_generator_();
    std::mt19937 generator(call_seed);
    CutMatching::Solver solver(
        instance.impl_->flow.get(),
        instance.impl_->subdivision.get(),
        &generator,
        &instance.impl_->subdivision_index,
        target_phi,
        certified_candidate_parameters());
    const CutMatching::Result practical = solver.compute(
        certified_candidate_parameters());
    const int expected_rounds = cut_matching_rounds(preview.terminal_count);
    const long long expected_scale = expected_rounds;
    const long double requested_scaled_capacity =
        static_cast<long double>(expected_scale) /
        static_cast<long double>(target_phi) /
        static_cast<long double>(expected_rounds);
    if (requested_scaled_capacity < 1.0L ||
        requested_scaled_capacity >
            static_cast<long double>(std::numeric_limits<long long>::max())) {
        return fallback_or_throw(
            "target phi cannot be represented by scaled integer cut matching");
    }
    const long long expected_capacity =
        static_cast<long long>(std::floor(requested_scaled_capacity));
    if (practical.roundsPlanned != expected_rounds ||
        practical.terminalCount != checked_size(
            preview.terminal_count, "terminal count audit overflow") ||
        practical.flowScale != expected_scale ||
        practical.perRoundCapacity != expected_capacity ||
        practical.iterations < 0 ||
        practical.iterations > practical.roundsPlanned) {
        return fallback_or_throw(
            "practical cut-matching execution metadata is inconsistent");
    }

    if (practical.type == CutMatching::Result::Balanced) {
        std::vector<int> side;
        side.reserve(static_cast<std::size_t>(instance.impl_->flow->size()));
        for (int local : *instance.impl_->flow) {
            side.push_back(instance.impl_->original_vertices.at(
                static_cast<std::size_t>(local)));
        }
        if (!side.empty() && side.size() < cluster.size()) {
            const ProofAugmentedCut verified = view.evaluate_proof_cut(side);
            if (verified.conductance >
                config_.gamma_cmp * target_phi + 1e-12) {
                return fallback_or_throw(
                    "practical balanced cut failed augmented conductance verification");
            }
            // The practical backend proves balance in its integral
            // ceil-rounded augmentation, whereas the proof-facing graph uses
            // the exact fractional augmentation weight. That balance claim
            // need not transfer between the two volume measures. Balance is
            // needed for the SW19/Goranci running-time analysis, but not for
            // the validity of a nontrivial sparse cut. We therefore accept
            // only the independently verified conductance statement here.
            // has_goranci_running_time() remains false for this adapter.
            return {
                OracleResult::Kind::SparseCut,
                std::move(side),
                verified.conductance,
                {}};
        }
        return fallback_or_throw(
            "practical balanced result returned a trivial cut");
    }

    if (practical.type == CutMatching::Result::Expander) {
        if (practical.iterations != practical.roundsPlanned ||
            !practical.allLogicalMatchingsComplete ||
            practical.congestion <= 0) {
            return fallback_or_throw(
                "practical expander did not complete the theoretical rounds");
        }
        // On small instances, turn the randomized conclusion into a
        // deterministic certificate whenever the exact fallback is enabled.
        if (config_.allow_exact_fallback &&
            cluster.size() <= config_.exact_fallback_max_vertices) {
            return ExactConductanceOracle(config_.exact_fallback_max_vertices)
                .analyze(graph, cluster, augmentation_weight, target_phi);
        }
        // Lemma 4.8(1), rather than 1/realized_congestion by itself, is the
        // expansion proof. Congestion is retained only as audit metadata: an
        // embedding certificate also needs the cut-matching graph's expansion.
        return {
            OracleResult::Kind::ExpanderCertificate, {}, target_phi,
            randomized_practical_proof(
                ExpansionProofBasis::RandomizedCutMatching,
                practical, call_seed, preview.terminal_count)};
    }

    if (practical.type == CutMatching::Result::NearExpander) {

        if (practical.iterations != practical.roundsPlanned ||
            practical.removedSubdivisionVolume > practical.balanceThreshold) {
            return fallback_or_throw(
                "near-expander did not reach the theoretical round/volume case");
        }

        std::vector<int> pretrim_side;
        pretrim_side.reserve(static_cast<std::size_t>(instance.impl_->flow->size()));
        for (int local : *instance.impl_->flow) {
            pretrim_side.push_back(instance.impl_->original_vertices.at(
                static_cast<std::size_t>(local)));
        }
        if (pretrim_side.empty() || pretrim_side.size() == cluster.size()) {
            return fallback_or_throw(
                "near-expander returned a trivial pre-trimming side; "
                "cluster_size=" + std::to_string(cluster.size()) +
                ", current_size=" +
                std::to_string(instance.impl_->flow->size()) +
                ", removed_size=" +
                std::to_string(instance.impl_->flow->removedSize()) +
                ", iterations=" + std::to_string(practical.iterations) +
                ", rounds=" + std::to_string(practical.roundsPlanned) +
                ", all_matchings_complete=" +
                std::to_string(practical.allLogicalMatchingsComplete));
        }
        const ProofAugmentedCut pretrim_cut =
            view.evaluate_proof_cut(pretrim_side);
        if (pretrim_cut.conductance >
            config_.gamma_cmp * target_phi + 1e-12) {
            return fallback_or_throw(
                "near-expander pre-trimming cut failed conductance verification");
        }
        const long double lemma_410_cut_limit =
            static_cast<long double>(target_phi) *
            view.proof_terminal_measure() / 16.0L;
        if (static_cast<long double>(pretrim_cut.cut_capacity) >
            lemma_410_cut_limit + 1e-12L) {
            return fallback_or_throw(
                "near-expander cut does not satisfy Lemma 4.10 Equation (1)");
        }

        const Trimming::Result trimming = Trimming::SaranurakWangTrimming(
            instance.impl_->flow.get(), target_phi, false);
        if (!trimming.feasible || trimming.remaining_vertices == 0) {
            return fallback_or_throw(
                "practical trimming did not certify a feasible nonempty remainder");
        }
        std::vector<int> side;
        side.reserve(trimming.remaining_vertices);
        for (int local : *instance.impl_->flow) {
            side.push_back(instance.impl_->original_vertices.at(
                static_cast<std::size_t>(local)));
        }
        if (side.empty() || side.size() == cluster.size()) {
            return fallback_or_throw(
                "practical trimming returned a trivial side");
        }
        const ProofAugmentedCut verified = view.evaluate_proof_cut(side);
        if (verified.conductance >
            config_.gamma_cmp * target_phi + 1e-12) {
            return fallback_or_throw(
                "trimmed cut failed augmented conductance verification");
        }
        const long double remainder_limit =
            view.proof_terminal_measure() / 10.0L;
        if (verified.complement_volume > remainder_limit) {
            return fallback_or_throw(
                "trimmed complement exceeds the Lemma 4.7 volume bound");
        }
        long double pruned_volume = 0.0L;
        std::vector<bool> retained(graph.node_count(), false);
        for (int vertex : side) {
            retained[vertex] = true;
        }
        for (int vertex : pretrim_side) {
            if (!retained[vertex]) {
                pruned_volume += view.proof_augmented_degree(vertex);
            }
        }
        const long double pruning_volume_limit =
            4.0L * static_cast<long double>(pretrim_cut.cut_capacity) /
            static_cast<long double>(target_phi);
        if (pruned_volume > pruning_volume_limit + 1e-12L) {
            return fallback_or_throw(
                "trimming violates the Lemma 4.9 pruned-volume bound");
        }
        if (static_cast<long double>(verified.cut_capacity) >
            2.0L * static_cast<long double>(pretrim_cut.cut_capacity)) {
            return fallback_or_throw(
                "trimming violates the Lemma 4.9 boundary-growth bound");
        }
        ExpansionProofAudit proof = randomized_practical_proof(
            ExpansionProofBasis::RandomizedCutMatchingAndTrimming,
            practical, call_seed, preview.terminal_count);
        double certified_lower_bound = target_phi;
        if (config_.allow_exact_fallback &&
            side.size() <= config_.exact_fallback_max_vertices) {
            const OracleResult exact_side =
                ExactConductanceOracle(config_.exact_fallback_max_vertices)
                    .analyze(
                        graph, side, augmentation_weight, target_phi);
            if (exact_side.kind != OracleResult::Kind::ExpanderCertificate) {
                return fallback_or_throw(
                    "trimmed side failed exact expansion verification");
            }
            proof = exact_side.proof;
            certified_lower_bound = exact_side.conductance;
        }
        return {
            OracleResult::Kind::CertifiedSide,
            std::move(side),
            certified_lower_bound,
            proof};
    }

    return fallback_or_throw(
        "fake-edge near-expander outcome is forbidden in certified mode");
}

std::string PracticalCutMatchingOracle::name() const {
    return config_.allow_exact_fallback
        ? "practical-cut-matching-trimming+exact-fallback"
        : "practical-cut-matching-trimming";
}
