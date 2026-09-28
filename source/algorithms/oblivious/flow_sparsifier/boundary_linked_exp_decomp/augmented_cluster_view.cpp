//
// Created by Mert Biyikli on 21.07.26.
//

#include "algorithms/oblivious/flow_sparsifier/boundary_linked_exp_decomp/augmented_cluster_view.h"
#include <cmath>
#include <limits>
#include <stdexcept>

constexpr std::size_t kNotInCluster = std::numeric_limits<std::size_t>::max();


Capacity checked_multiply(Capacity left, Capacity right, const char* message) {
    if (left != 0 && right > std::numeric_limits<Capacity>::max() / left) {
        throw std::overflow_error(message);
    }
    return left * right;
}


AugmentedClusterView::AugmentedClusterView(
    const Graph& graph,
    std::span<const int> cluster,
    double augmentation_weight)
    : graph_(&graph),
      cluster_(cluster.begin(), cluster.end()),
      position_by_vertex_(graph.node_count(), kNotInCluster),
      internal_degree_(cluster.size()),
      boundary_degree_(cluster.size()),
      augmented_degree_(cluster.size()),
      proof_augmented_degree_(cluster.size()),
      augmentation_weight_(augmentation_weight) {
    if (!std::isfinite(augmentation_weight_) || augmentation_weight_ < 0.0) {
        throw std::invalid_argument(
            "augmentation weight must be finite and nonnegative");
    }
    const long double rounded = std::ceil(
        static_cast<long double>(augmentation_weight_));
    if (rounded > static_cast<long double>(
            std::numeric_limits<Capacity>::max())) {
        throw std::overflow_error("boundary-loop multiplicity overflow");
    }
    boundary_loop_copies_ = static_cast<Capacity>(rounded);
    const long double backend_weight =
        static_cast<long double>(boundary_loop_copies_);
    const long double proof_weight =
        static_cast<long double>(augmentation_weight_);
    if (backend_weight < proof_weight ||
        backend_weight >= proof_weight + 1.0L) {
        throw std::logic_error(
            "rounded backend does not conservatively dominate proof model");
    }

    for (std::size_t i = 0; i < cluster_.size(); ++i) {
        const int v = cluster_[i];
        if (v < 0 || v >= graph.node_count()) {
            throw std::invalid_argument("cluster vertex is outside the graph");
        }
        if (position_by_vertex_[v] != kNotInCluster) {
            throw std::invalid_argument("cluster contains a duplicate vertex");
        }
        position_by_vertex_[v] = i;
    }

    for (const Edge& edge : graph.edges()) {
        const bool u_inside = position_by_vertex_[edge.u] != kNotInCluster;
        const bool v_inside = position_by_vertex_[edge.v] != kNotInCluster;
        if (u_inside && v_inside) {
            const std::size_t u = position_by_vertex_[edge.u];
            internal_degree_[u] = checked_add(
                internal_degree_[u], edge.multiplicity,
                "internal degree overflow");
            if (edge.u != edge.v) {
                const std::size_t v = position_by_vertex_[edge.v];
                internal_degree_[v] = checked_add(
                    internal_degree_[v], edge.multiplicity,
                    "internal degree overflow");
            }
            const auto kind = edge.u == edge.v
                ? TerminalBundle::Kind::Loop
                : TerminalBundle::Kind::InternalEdge;
            terminal_bundles_.push_back(
                {kind, edge.u, edge.v, edge.multiplicity});
            terminal_count_ = checked_add(
                terminal_count_, edge.multiplicity,
                "terminal count overflow");
            proof_terminal_measure_ +=
                static_cast<long double>(edge.multiplicity);
        } else if (u_inside || v_inside) {
            const int endpoint = u_inside ? edge.u : edge.v;
            const std::size_t index = position_by_vertex_[endpoint];
            boundary_degree_[index] = checked_add(
                boundary_degree_[index], edge.multiplicity,
                "boundary degree overflow");
            boundary_capacity_ = checked_add(
                boundary_capacity_, edge.multiplicity,
                "boundary capacity overflow");
            proof_terminal_measure_ +=
                static_cast<long double>(augmentation_weight_) *
                static_cast<long double>(edge.multiplicity);
            const Capacity loop_multiplicity = checked_multiply(
                boundary_loop_copies_, edge.multiplicity,
                "boundary terminal multiplicity overflow");
            if (loop_multiplicity != 0) {
                terminal_bundles_.push_back({
                    TerminalBundle::Kind::Loop,
                    endpoint,
                    endpoint,
                    loop_multiplicity});
                terminal_count_ = checked_add(
                    terminal_count_, loop_multiplicity,
                    "terminal count overflow");
            }
        }
    }

    for (std::size_t i = 0; i < cluster_.size(); ++i) {
        const Capacity reconstructed_ordinary_degree = checked_add(
            internal_degree_[i], boundary_degree_[i],
            "ordinary degree reconstruction overflow");
        if (reconstructed_ordinary_degree != graph.degree(cluster_[i])) {
            throw std::logic_error("cluster degree accounting invariant failed");
        }
        const Capacity boundary_contribution = checked_multiply(
            boundary_loop_copies_, boundary_degree_[i],
            "augmented degree overflow");
        augmented_degree_[i] = checked_add(
            internal_degree_[i], boundary_contribution,
            "augmented degree overflow");
        augmented_volume_ = checked_add(
            augmented_volume_, augmented_degree_[i],
            "augmented volume overflow");
        proof_augmented_degree_[i] =
            static_cast<long double>(internal_degree_[i]) +
            static_cast<long double>(augmentation_weight_) *
                static_cast<long double>(boundary_degree_[i]);
        proof_augmented_volume_ += proof_augmented_degree_[i];
        ordinary_volume_ = checked_add(
            ordinary_volume_, graph.degree(cluster_[i]),
            "ordinary volume overflow");
    }

    if (!std::isfinite(proof_augmented_volume_) ||
        !std::isfinite(proof_terminal_measure_)) {
        throw std::overflow_error("fractional augmented volume overflow");
    }

    Capacity routing_terminal_count = 0;
    Capacity loop_terminal_count = 0;
    for (const TerminalBundle& bundle : terminal_bundles_) {
        Capacity& total = bundle.kind == TerminalBundle::Kind::InternalEdge
            ? routing_terminal_count
            : loop_terminal_count;
        total = checked_add(total, bundle.multiplicity,
            "terminal-kind count overflow");
    }
    const Capacity reconstructed_augmented_volume = checked_add(
        checked_multiply(2, routing_terminal_count,
            "augmented volume reconstruction overflow"),
        loop_terminal_count,
        "augmented volume reconstruction overflow");
    if (reconstructed_augmented_volume != augmented_volume_) {
        throw std::logic_error("augmented volume accounting invariant failed");
    }
}

std::size_t AugmentedClusterView::position(int v) const {
    if (v < 0 || static_cast<std::size_t>(v) >= position_by_vertex_.size() ||
        position_by_vertex_[v] == kNotInCluster) {
        throw std::invalid_argument("vertex is not in the augmented cluster");
    }
    return position_by_vertex_[v];
}

Capacity AugmentedClusterView::internal_degree(int v) const {
    return internal_degree_[position(v)];
}

Capacity AugmentedClusterView::boundary_degree(int v) const {
    return boundary_degree_[position(v)];
}

Capacity AugmentedClusterView::augmented_degree(int v) const {
    return augmented_degree_[position(v)];
}

long double AugmentedClusterView::proof_augmented_degree(int v) const {
    return proof_augmented_degree_[position(v)];
}

std::vector<ExplicitTerminal> AugmentedClusterView::materialize_terminals(
    std::size_t maximum_terminals) const {
    if (terminal_count_ > maximum_terminals) {
        throw std::runtime_error("explicit terminal materialization limit exceeded");
    }
    std::vector<ExplicitTerminal> result;
    result.reserve(static_cast<std::size_t>(terminal_count_));
    for (const TerminalBundle& bundle : terminal_bundles_) {
        for (Capacity copy = 0; copy < bundle.multiplicity; ++copy) {
            result.push_back({bundle.kind, bundle.u, bundle.v});
        }
    }
    return result;
}

ExplicitAugmentedInstance AugmentedClusterView::build_explicit_instance(
    std::size_t maximum_terminals) const {
    const std::vector<ExplicitTerminal> terminals =
        materialize_terminals(maximum_terminals);
    ExplicitAugmentedInstance result;
    result.original_vertices = cluster_;
    result.augmented_degrees = augmented_degree_;
    result.terminal_count = terminals.size();
    result.routing_edges.reserve(terminals.size());
    if (terminals.size() > std::numeric_limits<std::size_t>::max() / 2) {
        throw std::overflow_error("subdivision incidence count overflow");
    }
    result.terminal_incidences.reserve(2 * terminals.size());

    for (std::size_t terminal = 0; terminal < terminals.size(); ++terminal) {
        const ExplicitTerminal& item = terminals[terminal];
        const std::size_t u = position(item.u);
        result.terminal_incidences.push_back({terminal, u});
        if (item.kind == TerminalBundle::Kind::InternalEdge) {
            const std::size_t v = position(item.v);
            if (u == v) {
                throw std::logic_error(
                    "internal routing terminal unexpectedly has one endpoint");
            }
            result.routing_edges.push_back({u, v});
            result.terminal_incidences.push_back({terminal, v});
        }
    }
    return result;
}

AugmentedCut AugmentedClusterView::evaluate_cut(
    std::span<const int> side) const {
    std::vector<bool> selected(cluster_.size(), false);
    Capacity side_volume = 0;
    for (int v : side) {
        const std::size_t index = position(v);
        if (selected[index]) {
            throw std::invalid_argument("cut side contains a duplicate vertex");
        }
        selected[index] = true;
        side_volume = checked_add(
            side_volume, augmented_degree_[index],
            "cut-side volume overflow");
    }

    Capacity cut_capacity = 0;
    for (const Edge& edge : graph_->edges()) {
        if (edge.u == edge.v ||
            position_by_vertex_[edge.u] == kNotInCluster ||
            position_by_vertex_[edge.v] == kNotInCluster) {
            continue;
        }
        if (selected[position_by_vertex_[edge.u]] !=
            selected[position_by_vertex_[edge.v]]) {
            cut_capacity = checked_add(
                cut_capacity, edge.multiplicity,
                "cut capacity overflow");
        }
    }
    const Capacity complement_volume = augmented_volume_ - side_volume;
    const Capacity denominator = std::min(side_volume, complement_volume);
    const double conductance = denominator == 0
        ? std::numeric_limits<double>::infinity()
        : static_cast<double>(static_cast<long double>(cut_capacity) /
              static_cast<long double>(denominator));
    return {conductance, cut_capacity, side_volume, complement_volume};
}

ProofAugmentedCut AugmentedClusterView::evaluate_proof_cut(
    std::span<const int> side) const {
    std::vector<bool> selected(cluster_.size(), false);
    long double side_volume = 0.0L;
    for (int v : side) {
        const std::size_t index = position(v);
        if (selected[index]) {
            throw std::invalid_argument("cut side contains a duplicate vertex");
        }
        selected[index] = true;
        side_volume += proof_augmented_degree_[index];
    }

    Capacity cut_capacity = 0;
    for (const Edge& edge : graph_->edges()) {
        if (edge.u == edge.v ||
            position_by_vertex_[edge.u] == kNotInCluster ||
            position_by_vertex_[edge.v] == kNotInCluster) {
            continue;
        }
        if (selected[position_by_vertex_[edge.u]] !=
            selected[position_by_vertex_[edge.v]]) {
            cut_capacity = checked_add(
                cut_capacity, edge.multiplicity,
                "cut capacity overflow");
        }
    }

    const long double complement_volume =
        proof_augmented_volume_ - side_volume;
    const long double denominator =
        std::min(side_volume, complement_volume);
    const double conductance = denominator <= 0.0L
        ? std::numeric_limits<double>::infinity()
        : static_cast<double>(
              static_cast<long double>(cut_capacity) / denominator);
    return {conductance, cut_capacity, side_volume, complement_volume};
}
