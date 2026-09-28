//
// Created by Mert Biyikli on 21.07.26.
//

#ifndef OBLIVIOUSROUTING_PRACTICAL_ADAPTER_H
#define OBLIVIOUSROUTING_PRACTICAL_ADAPTER_H

#include "augmented_cluster_view.h"
#include "expander_decomp.h"

namespace UnitFlow {
    class Graph;
}

struct PracticalInstanceSummary {
    std::size_t original_vertices{};
    std::size_t routing_edges{};
    std::size_t terminals{};
    std::size_t subdivision_vertices{};
    std::size_t subdivision_incidences{};
    Capacity augmented_volume{};
};

class PracticalFlowInstance {
public:
    PracticalFlowInstance(PracticalFlowInstance&&) noexcept;
    PracticalFlowInstance& operator=(PracticalFlowInstance&&) noexcept;
    ~PracticalFlowInstance();

    PracticalFlowInstance(const PracticalFlowInstance&) = delete;
    PracticalFlowInstance& operator=(const PracticalFlowInstance&) = delete;

    [[nodiscard]] static PracticalFlowInstance build(
        const AugmentedClusterView& view,
        std::size_t maximum_explicit_terminals);

    [[nodiscard]] const PracticalInstanceSummary& summary() const noexcept;
    [[nodiscard]] UnitFlow::Graph& flow_graph() noexcept;
    [[nodiscard]] UnitFlow::Graph& subdivision_graph() noexcept;

private:
    struct Impl;
    explicit PracticalFlowInstance(std::unique_ptr<Impl> impl);
    std::unique_ptr<Impl> impl_;

    friend class PracticalCutMatchingOracle;
};

enum class OracleSelectionMode {
    Auto,
    ExactOnly,
    PracticalOnly
};

struct PracticalOracleConfig {
    std::size_t maximum_explicit_terminals{2'000'000};
    std::size_t exact_fallback_max_vertices{20};
    std::uint32_t random_seed{555};
    double gamma_cmp{1.0};
    bool allow_exact_fallback{false};
    OracleSelectionMode selection_mode{OracleSelectionMode::Auto};
};

// Non-heuristic proof-of-concept implementation of Goranci Lemma 4.7 using the
// practical cut-matching and Saranurak-Wang trimming primitives. Candidate cuts
// and all deterministic side conditions are checked independently. The optional
// exact fallback is off by default and intended only for differential tests.
class PracticalCutMatchingOracle final : public CertifiedSparseCutOracle {
public:
    explicit PracticalCutMatchingOracle(PracticalOracleConfig config = {});

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
    PracticalOracleConfig config_;
    // A fresh derived seed is consumed for every adaptive oracle call. Reusing
    // the identical random stream from its beginning on every cluster would
    // not match the randomized cut-matching analysis.
    mutable std::mt19937 seed_generator_;
};
#endif //OBLIVIOUSROUTING_PRACTICAL_ADAPTER_H