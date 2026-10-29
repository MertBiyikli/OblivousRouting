//
// Created by Mert on 22.10.26.
//

#ifndef E_ROUTING_SCENARIO_ANALYSIS_H
#define E_ROUTING_SCENARIO_ANALYSIS_H

#include "core/errors.h"
#include "core/types.h"
#include "data_structures/graph/graph.h"
#include "io/solver_io.h"
#include "scenario_results.h"

struct ScenarioContext {
    ScenarioContext(const optimized::Graph<EdgeData>& g) : graph(g), fingerprint(optimized::GraphFingerprint::compute(g)) {}
    const optimized::Graph<EdgeData>& graph;
    std::uint64_t fingerprint;
};

// The callback must capture ONE fixed demand matrix and evaluate either graph
// using the corresponding computed routing result.
using RoutingMetricsEvaluator = std::function<Result<NetworkData>(const optimized::Graph<EdgeData>&, const IRoutingResult&)>;

class Scenario {
  public:
    virtual ~Scenario() = default;

    virtual ScenarioKind kind() const noexcept = 0;
    // Register modifications from parsed configuration.
    virtual Result<void> configure(const Config& cfg) = 0;

    virtual Result<void> verify(const optimized::Graph<EdgeData>& graph) const = 0;

    virtual Result<void> apply(optimized::Graph<EdgeData>& graph) const = 0;
};



inline std::vector<SolverType> getSolverBasedOnExecutionMode(ExecutionMode mode) {
    switch (mode) {
    case ExecutionMode::FAST:
        return {SolverType::EXPANDER_HIERARCHY};
    case ExecutionMode::STANDARD:
        return {SolverType::EXPANDER_HIERARCHY, SolverType::ELECTRICAL_SKETCHING, SolverType::RAECKE_CKR_FLAT};
    case ExecutionMode::RESEARCH:
        return {SolverType::EXPANDER_HIERARCHY, SolverType::ELECTRICAL_SKETCHING, SolverType::RAECKE_CKR_FLAT, SolverType::RAECKE_RANDOM_MST_FLAT, SolverType::LP_APPLEGATE_COHEN};
    default:
        return {SolverType::ELECTRICAL_SKETCHING};
    }
}

inline std::vector<SolverType> getSolverFromConfigAndMode(const Config& cfg, ExecutionMode mode) {

    switch (mode) {
        case ExecutionMode::FAST:
            if (cfg.solvers.empty()) {
                return getSolverBasedOnExecutionMode(mode);
            }else {
                return {cfg.solvers.front()};
            }
        case ExecutionMode::STANDARD:
            if (!cfg.solvers.empty()) {
                return cfg.solvers;
            }
            return getSolverBasedOnExecutionMode(mode);
        case ExecutionMode::RESEARCH:
            return getSolverBasedOnExecutionMode(mode);
    }
}
inline Result<void> solveInto(SolverType type, optimized::Graph<EdgeData>& working_graph, std::optional<IRoutingResult>& destination, double& seconds) {
    auto solver_opt = makeSolver(type, working_graph);
    if (!solver_opt || !(*solver_opt)) {
        return makeErrorMessage(ErrorCode::InvalidSolver, getSolverName(type) + " solver was not found.");
    }

    const auto started = std::chrono::steady_clock::now();
    auto scheme = (*solver_opt)->solve();
    seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();

    if (!scheme)
        return getError(scheme);
    if (!(*scheme)) {
        return makeErrorMessage(ErrorCode::SolverFailed, "The computed scheme is null.");
    }
    IRoutingResult routing;
    routing.scheme = std::move(scheme.value());
    destination.emplace(std::move(routing));
    return {};
}

inline Result<NetworkData> evaluateFixedRoutingUpgrade(const NetworkData& baseline, const optimized::Graph<EdgeData>& upgraded_graph) {
    NetworkData result;
    result.total_demand_volume = baseline.total_demand_volume;
    result.links.reserve(baseline.links.size());

    const int m = upgraded_graph.getNumDirectedEdges();
    double utilization_sum = 0.0;

    for (const auto& link : baseline.links) {
        const int e = link.edge_id;

        if (e < 0 || e >= m) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid edge ID in baseline Data.");
        }

        const int reverse = upgraded_graph.reverse(e).id;

        if (reverse < 0 || reverse >= m || upgraded_graph.reverse(reverse).id != e) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid reverse edge mapping.");
        }

        const double forward_capacity = upgraded_graph.edgeData(e).capacity;

        const double reverse_capacity = upgraded_graph.edgeData(reverse).capacity;

        if (!std::isfinite(forward_capacity) || forward_capacity <= 0.0 || !std::isfinite(reverse_capacity) || reverse_capacity <= 0.0 || !std::isfinite(link.load) || link.load < 0.0) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid capacity or physical-link load.");
        }

        // Same routing and demands => identical physical load.
        const double load = link.load;

        const double utilization = std::max(load / forward_capacity, load / reverse_capacity);

        if (!std::isfinite(utilization)) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid upgraded link utilization.");
        }

        result.links.push_back({e, forward_capacity, load, utilization});

        result.max_congestion = std::max(result.max_congestion, utilization);

        utilization_sum += utilization;
    }

    if (!result.links.empty()) {
        result.average_congestion = utilization_sum / result.links.size();
    }

    return result;
}

// TODO: below
/*
class MultiLinkFailure : public Scenario {
  public:
    MultiLinkFailure() : Scenario() {}
    Result<void> verify() override;
};

class CapacityPlanning : public Scenario {
  public:
    CapacityPlanning(const optimized::Graph<EdgeData>& g) : Scenario(g) {}
    Result<void> verify() override;
};
*/
#endif // E_ROUTING_SCENARIO_ANALYSIS_H