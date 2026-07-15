//
// Created by Mert Biyikli on 13.07.26.
//

#ifndef OBLIVIOUSROUTING_TREE_FLOW_ELECTRICAL_EMBEDDER_H
#define OBLIVIOUSROUTING_TREE_FLOW_ELECTRICAL_EMBEDDER_H

#pragma once

#include "hierarchy_results.h"

#include "tree_flow_router.h"
#include "tree_sparsifier.h"
#include "routing_oracle/local_electrical_solver.h"
#include "core/errors.h"
#include "data_structures/graph/Igraph.h"

#include <unordered_map>
#include <utility>
#include <vector>

struct ClusterBoundaryProfile {
    int cluster_id = -1;

    double total_capacity = 0.0;

    /*
     * For each vertex v in the cluster:
     *
     *     sum of capacities of edges in delta(cluster)
     *     incident to v.
     */
    std::vector<std::pair<int, double>> vertex_capacity;
};

struct ClusterEmbeddingInfo {
    int cluster_id = -1;
    int level = -1;

    double outgoing_tree_flow = 0.0;
    double local_demand_l1 = 0.0;

    double local_max_congestion = 0.0;
    double local_conservation_error = 0.0;

    bool electrical_solve_performed = false;
};

struct ElectricalEmbeddingResult {
    /*
     * Final feasible signed flow on original graph edges.
     *
     * Indexed by directed edge ID. Only the canonical representative
     * used by cluster.induced_edges normally has a nonzero value.
     */
    std::vector<double> signed_edge_flow;

    /*
     * Sum of absolute flow contributed by every cluster solve.
     *
     * This gives the conservative concatenation congestion used by the
     * hierarchy analysis. Different-level flows may algebraically cancel
     * in signed_edge_flow but not in edge_usage.
     */
    std::vector<double> edge_usage;

    double max_congestion = 0.0;
    double max_net_congestion = 0.0;
    double max_conservation_error = 0.0;

    int bottleneck_edge = -1;

    int electrical_solves = 0;

    std::vector<ClusterEmbeddingInfo> cluster_info;
};

class TreeFlowElectricalEmbedder {
public:
    TreeFlowElectricalEmbedder(const IGraph& graph,const HierarchyResult& hierarchy,const TreeSparsifier& tree)
        : graph_(graph),
          hierarchy_(hierarchy),
          tree_(tree) {}

    Result<ElectricalEmbeddingResult> embed(const TreeFlowResult& tree_flow) const;

private:
    const IGraph& graph_;
    const HierarchyResult& hierarchy_;
    const TreeSparsifier& tree_;

    LocalElectricalSolver electrical_solver_;

    Result<std::unordered_map<int, ClusterBoundaryProfile>> buildBoundaryProfiles() const;

    Result<double> clusterTreeFlow(int cluster, const TreeFlowResult& tree_flow) const;

    Result<double> vertexTreeFlow(int vertex, const TreeFlowResult& tree_flow) const;
};

#endif //OBLIVIOUSROUTING_TREE_FLOW_ELECTRICAL_EMBEDDER_H