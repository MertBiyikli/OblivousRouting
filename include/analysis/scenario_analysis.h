//
// Created by Mert on 22.10.26.
//

#ifndef E_ROUTING_SCENARIO_ANALYSIS_H
#define E_ROUTING_SCENARIO_ANALYSIS_H

#include "core/errors.h"
#include "core/types.h"
#include "data_structures/graph/graph.h"
#include "scenario_results.h"

#include <unordered_set>

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

struct ScenarioContext {
    ScenarioContext(const optimized::Graph<EdgeData>& g) : graph(g), fingerprint(optimized::GraphFingerprint::compute(g)) {}
    const optimized::Graph<EdgeData>& graph;
    std::uint64_t fingerprint;
};

class Scenario {
  public:
    virtual ~Scenario() = default;
    optimized::Graph<EdgeData> graph;
    std::uint64_t fingerprint;
    Scenario(const optimized::Graph<EdgeData>& g) : graph(g) { fingerprint = optimized::GraphFingerprint::compute(graph); }

    virtual Result<ScenarioResult> execute(const ScenarioContext& context, ExecutionMode mode);
    virtual Result<void> verify() = 0;

    Result<void> validEdge(int edge_id) const {
        if (edge_id < 0 || edge_id >= graph.getNumDirectedEdges()) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Edge ID must be between 0 and m (for which each edge is represented by two directed arcs.");
        }
        return {};
    }

    Result<void> validNode(int node_id) {
        if (node_id < 0 || node_id >= graph.getNumNodes()) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Node id must be in between 0 and n.");
        }
        return {};
    }

    int canonicalEdge(int edge_id) const {
        if (!validEdge(edge_id)) {
            return -1;
        }

        const int canonical_edge_id = (graph.getEdgeEndpoints(edge_id).first < graph.getEdgeEndpoints(edge_id).second ? edge_id : graph.reverse(edge_id).id);
        return canonical_edge_id;
    }
};

class LinkFailure : public Scenario {
  public:
    LinkFailure(const optimized::Graph<EdgeData>& g) : Scenario(g) {}
    Result<void> verify() override;

    Result<void> linkFailure(int edge_id);

    std::unordered_set<int> failed_links;
};


class TrafficGrowth : public Scenario {
  public:
    TrafficGrowth(const optimized::Graph<EdgeData>& g) : Scenario(g) {}
    Result<void> verify() override;

    Result<void> setGlobalScale(double scale);
    Result<void> growPerCommodity(int s, int t, double scale);
    double globalScale{};
    std::unordered_map<std::pair<int, int>, double, PairHash> commodities;
};

// The callback must capture ONE fixed demand matrix and evaluate either graph
// using the corresponding computed routing result.
using RoutingMetricsEvaluator = std::function<Result<NetworkData>(const optimized::Graph<EdgeData>&, const IRoutingResult&)>;

class CapacityUpgrade : public Scenario {
  public:
    CapacityUpgrade(const optimized::Graph<EdgeData>& g) : Scenario(g) {}
    Result<void> verify() override;
    Result<ScenarioResult> execute(const ScenarioContext& context, ExecutionMode mode) override;

    Result<void> upgrade(const Config& cfg);
    Result<void> upgrade(int edge_id, double scale);
    Result<void> apply(optimized::Graph<EdgeData>& target) const;

    Result<ScenarioResult> executeEvaluated(const ScenarioContext& context,
                                            ExecutionMode mode,
                                            const RoutingMetricsEvaluator& evaluator);

    // Direct evaluation using one fixed demand instance.
    Result<ScenarioResult> executeWithDemands(const ScenarioContext& context,
                                              ExecutionMode mode,
                                              const demands& fixed_demands,
                                              double total_demand_volume);


    std::unordered_map<int, double> upgrades;
};

class MultiLinkFailure : public Scenario {
  public:
    MultiLinkFailure(const optimized::Graph<EdgeData>& g) : Scenario(g) {}
    Result<void> verify() override;
};

class CapacityPlanning : public Scenario {
  public:
    CapacityPlanning(const optimized::Graph<EdgeData>& g) : Scenario(g) {}
    Result<void> verify() override;
};

#endif // E_ROUTING_SCENARIO_ANALYSIS_H