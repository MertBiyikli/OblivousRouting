//
// Created by Mert Biyikli on 13.07.26.
//

#ifndef OBLIVIOUSROUTING_TREE_FLOW_ROUTER_H
#define OBLIVIOUSROUTING_TREE_FLOW_ROUTER_H
#pragma once

#include "tree_sparsifier.h"

#include "core/errors.h"
#include "utils/demands.h"

#include <vector>

struct TreePairDemand {
    int source = -1;
    int target = -1;
    double value = 0.0;
};

/**
 * direction:
 *
 *   +1 = child -> parent
 *   -1 = parent -> child
 *
 * TreeSparsifierEdge itself is stored child -> parent.
 */
struct TreePathStep {
    int edge_id = -1;
    int direction = 0;
};

struct TreeCommodityPath {
    int source = -1;
    int target = -1;
    double value = 0.0;

    int lca = -1;

    // Ordered from source leaf to target leaf.
    std::vector<TreePathStep> steps;
};

struct TreeFlowResult {
    /*
     * Algebraic flow relative to the stored child -> parent orientation.
     *
     * Positive:
     *     child -> parent
     *
     * Negative:
     *     parent -> child
     */
    std::vector<double> signed_edge_flow;

    /*
     * Commodity load on each edge.
     *
     * This is the sum of the magnitudes of all commodities traversing
     * the edge. This, not abs(signed_edge_flow), determines congestion.
     */
    std::vector<double> edge_load;

    std::vector<double> edge_congestion;

    double max_congestion = 0.0;
    double max_conservation_error = 0.0;

    int bottleneck_edge = -1;

    // Optional, primarily useful for smoke tests and debugging.
    std::vector<TreeCommodityPath> commodity_paths;
};

/**
 * Computes the unique route in the tree sparsifier that connects two leaf nodes. The idea is to identify
 * the LCA and store the orientation of the paths from the source -> lca and lca -> target. This represent the
 * inter-cluster flow in the tree sparsifier.
 */
class TreeFlowRouter {
public:
    explicit TreeFlowRouter(const TreeSparsifier& tree)
        : tree_(tree) {}


    Result<TreeFlowResult> routePair(int source,int target,double value = 1.0) const;


    Result<TreeFlowResult> route(const std::vector<TreePairDemand>& pair_demands,bool record_paths = false) const;

    /*
     * Adapter for your existing demand representation:
     *
     *     map<pair<int, int>, double>
     */

    Result<TreeFlowResult> route(const demands& demand,bool record_paths = false) const;

private:
    const TreeSparsifier& tree_;


    Result<void> addPairDemand(const TreePairDemand& demand,TreeFlowResult& result,std::vector<double>& expected_divergence,bool record_path) const;


    Result<void> finalizeAndValidate(TreeFlowResult& result,const std::vector<double>& expected_divergence) const;
};
#endif //OBLIVIOUSROUTING_TREE_FLOW_ROUTER_H