//
// Created by Mert Biyikli on 21.07.26.
//

#ifndef OBLIVIOUSROUTING_AUGMENTED_CLUSTER_VIEW_H
#define OBLIVIOUSROUTING_AUGMENTED_CLUSTER_VIEW_H
#include "expander_decomp.h"

#include <cstddef>
#include <span>
#include <vector>

struct AugmentedCut {
    double conductance{};
    Capacity cut_capacity{};
    Capacity side_volume{};
    Capacity complement_volume{};
};

// Proof-facing cut in the exact fractional graph G[U]^w.  The practical
// backend still receives an integral graph obtained by rounding w upward, but
// no cut is accepted using those rounded volumes.
struct ProofAugmentedCut {
    double conductance{};
    Capacity cut_capacity{};
    long double side_volume{};
    long double complement_volume{};
};

struct TerminalBundle {
    enum class Kind { InternalEdge, Loop };
    Kind kind{Kind::InternalEdge};
    int u{};
    int v{};  // Equal to u for a loop terminal.
    Capacity multiplicity{};
};

struct ExplicitTerminal {
    TerminalBundle::Kind kind{TerminalBundle::Kind::InternalEdge};
    int u{};
    int v{};
};

struct LocalRoutingEdge {
    std::size_t u{};
    std::size_t v{};
};

struct TerminalIncidence {
    std::size_t terminal{};
    std::size_t vertex{};
};

// Neutral input matching the practical repository's two-graph construction:
// routing_edges build flowGraph; terminal_incidences build the subdivision
// graph with terminal vertex local_vertex_count + terminal.
struct ExplicitAugmentedInstance {
    std::vector<int> original_vertices;
    std::vector<Capacity> augmented_degrees;
    std::vector<LocalRoutingEdge> routing_edges;
    std::vector<TerminalIncidence> terminal_incidences;
    std::size_t terminal_count{};
};

// A checked, immutable representation of H = G[U]^w.  It deliberately keeps
// two volume models: an integral, rounded-up model used only by the practical
// UnitFlow backend, and an exact fractional model used by every proof-facing
// conductance check.
class AugmentedClusterView {
public:
    AugmentedClusterView(
        const Graph& graph,
        std::span<const int> cluster,
        double augmentation_weight);

    [[nodiscard]] const Graph& graph() const noexcept { return *graph_; }
    [[nodiscard]] std::span<const int> vertices() const noexcept {
        return cluster_;
    }
    [[nodiscard]] double augmentation_weight() const noexcept {
        return augmentation_weight_;
    }
    [[nodiscard]] Capacity boundary_loop_copies() const noexcept {
        return boundary_loop_copies_;
    }

    [[nodiscard]] Capacity internal_degree(int v) const;
    [[nodiscard]] Capacity boundary_degree(int v) const;
    [[nodiscard]] Capacity augmented_degree(int v) const;
    [[nodiscard]] long double proof_augmented_degree(int v) const;

    [[nodiscard]] Capacity ordinary_volume() const noexcept {
        return ordinary_volume_;
    }
    [[nodiscard]] Capacity augmented_volume() const noexcept {
        return augmented_volume_;
    }
    [[nodiscard]] long double proof_augmented_volume() const noexcept {
        return proof_augmented_volume_;
    }
    [[nodiscard]] long double proof_terminal_measure() const noexcept {
        return proof_terminal_measure_;
    }
    [[nodiscard]] Capacity boundary_capacity() const noexcept {
        return boundary_capacity_;
    }
    [[nodiscard]] Capacity terminal_count() const noexcept {
        return terminal_count_;
    }

    [[nodiscard]] std::span<const TerminalBundle> terminal_bundles() const noexcept {
        return terminal_bundles_;
    }

    // Intended for the first, explicit practical cut-matching adapter. The
    // caller must choose a memory guard; exceeding it fails closed.
    [[nodiscard]] std::vector<ExplicitTerminal> materialize_terminals(
        std::size_t maximum_terminals) const;

    [[nodiscard]] ExplicitAugmentedInstance build_explicit_instance(
        std::size_t maximum_terminals) const;

    [[nodiscard]] AugmentedCut evaluate_cut(
        std::span<const int> side) const;

    [[nodiscard]] ProofAugmentedCut evaluate_proof_cut(
        std::span<const int> side) const;

private:
    [[nodiscard]] std::size_t position(int v) const;

    const Graph* graph_{};
    std::vector<int> cluster_;
    std::vector<std::size_t> position_by_vertex_;
    std::vector<Capacity> internal_degree_;
    std::vector<Capacity> boundary_degree_;
    std::vector<Capacity> augmented_degree_;
    std::vector<long double> proof_augmented_degree_;
    std::vector<TerminalBundle> terminal_bundles_;
    double augmentation_weight_{};
    Capacity boundary_loop_copies_{};
    Capacity ordinary_volume_{};
    Capacity augmented_volume_{};
    long double proof_augmented_volume_{};
    long double proof_terminal_measure_{};
    Capacity boundary_capacity_{};
    Capacity terminal_count_{};
};


#endif //OBLIVIOUSROUTING_AUGMENTED_CLUSTER_VIEW_H
