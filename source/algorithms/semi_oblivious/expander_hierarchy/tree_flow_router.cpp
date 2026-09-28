//
// Created by Mert Biyikli on 13.07.26.
//

#include "algorithms/semi_oblivious/expander_hierarchy/tree_flow_router.h"
#include "core/errors.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>



constexpr double kFlowEpsilon = 1e-10;

bool validNode(const TreeSparsifier& tree,const int node) {
    return node >= 0 && node < static_cast<int>(tree.nodes.size());
}

bool validEdge(const TreeSparsifier& tree,const int edge) {
    return edge >= 0 &&
           edge < static_cast<int>(tree.edges.size());
}



Result<TreeFlowResult> TreeFlowRouter::routePair(const int source,const int target,const double value) const {
    std::vector<TreePairDemand> pair_demands{
        TreePairDemand{
            .source = source,
            .target = target,
            .value = value
        }
    };

    return route(pair_demands, true);
}

Result<TreeFlowResult> TreeFlowRouter::route(const demands& demand,const bool record_paths) const {
    std::vector<TreePairDemand> pair_demands;
    pair_demands.reserve(demand.size());

    for (int i = 0; i<demand.size(); i++) {
        double value = demand.getDemandValue(i);
        const auto [source, target] = demand.getDemandPair(i);

        pair_demands.push_back(TreePairDemand{.source = source,.target = target,.value = value});
    }


    return route(pair_demands, record_paths);
}

Result<TreeFlowResult> TreeFlowRouter::route(const std::vector<TreePairDemand>& pair_demands,const bool record_paths) const {
    if (tree_.empty()) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Cannot route on an empty tree sparsifier.");
    }

    if (!validNode(tree_, tree_.root)) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Tree sparsifier has an invalid root.");
    }

    if (tree_.vertex_to_leaf.empty()) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Tree sparsifier has no original-vertex leaf mapping.");
    }

    TreeFlowResult result;

    result.signed_edge_flow.assign(tree_.edges.size(),0.0);
    result.edge_load.assign(tree_.edges.size(),0.0);
    result.edge_congestion.assign(tree_.edges.size(),0.0);

    if (record_paths) {
        result.commodity_paths.reserve(pair_demands.size());
    }

    /*
     * Expected divergence is indexed by tree node.
     *
     * For a pair s -> t:
     *
     *     +d at the tree leaf representing s
     *     -d at the tree leaf representing t
     *      0 at every internal tree node
     */
    std::vector<double> expected_divergence(tree_.nodes.size(),0.0);

    for (const auto& pair_demand : pair_demands) {
        auto added = addPairDemand(pair_demand,result,expected_divergence,record_paths);

        if (!added) {
            return getError(added);
        }
    }

    auto finalized = finalizeAndValidate(result,expected_divergence);

    if (!finalized) {
        return getError(finalized);
    }

    return result;
}

Result<void> TreeFlowRouter::addPairDemand(const TreePairDemand& demand,TreeFlowResult& result,std::vector<double>& expected_divergence,const bool record_path) const {
    if (demand.source < 0 ||demand.target < 0 ||demand.source >=static_cast<int>(tree_.vertex_to_leaf.size()) ||demand.target >=static_cast<int>(tree_.vertex_to_leaf.size())) {
        return makeErrorMessage(ErrorCode::InvalidDemand,"Tree demand contains an invalid original vertex.");
    }

    if (!std::isfinite(demand.value) ||demand.value < 0.0) {
        return makeErrorMessage(ErrorCode::InvalidDemand,"Tree demand value must be finite and non-negative.");
    }

    if (demand.source == demand.target || demand.value <= kFlowEpsilon) {
        if (record_path) {
            result.commodity_paths.push_back(
                TreeCommodityPath{
                    .source = demand.source,
                    .target = demand.target,
                    .value = demand.value,
                    .lca = tree_.vertex_to_leaf[demand.source],
                    .steps = {}
                }
            );
        }

        return {};
    }

    int source_node = tree_.vertex_to_leaf[demand.source];
    int target_node = tree_.vertex_to_leaf[demand.target];

    if (!validNode(tree_, source_node) || !validNode(tree_, target_node)) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Tree vertex-to-leaf mapping contains an invalid node.");
    }

    const int source_leaf = source_node;
    const int target_leaf = target_node;

    expected_divergence[source_leaf] += demand.value;
    expected_divergence[target_leaf] -= demand.value;

    std::vector<TreePathStep> source_steps;
    std::vector<TreePathStep> target_steps;

    /*
     * Move one node upward.
     *
     * side_direction = +1:
     *     the commodity moves child -> parent.
     *
     * side_direction = -1:
     *     the eventual source-to-target path moves parent -> child.
     *     We are temporarily traversing that side upward only to find
     *     the LCA.
     */
    auto moveUp = [&](int& node,const int side_direction,std::vector<TreePathStep>& steps) -> Result<void> {
            if (!validNode(tree_, node)) {
                return makeErrorMessage(ErrorCode::InvalidGraph,"Invalid node encountered while finding tree path.");
            }

            const auto& tree_node = tree_.nodes[node];

            if (!tree_node.parent.has_value() ||!tree_node.parent_edge.has_value()) {
                return makeErrorMessage(ErrorCode::InvalidGraph,"Reached a non-root tree node without a parent.");
            }

            const int parent = *tree_node.parent;

            const int edge_id = *tree_node.parent_edge;

            if (!validNode(tree_, parent) || !validEdge(tree_, edge_id)) {
                return makeErrorMessage(ErrorCode::InvalidGraph,"Tree node contains an invalid parent reference.");
            }

            const auto& edge = tree_.edges[edge_id];

            if (edge.child != node || edge.parent != parent) {
                return makeErrorMessage(ErrorCode::InvalidGraph,"Tree parent edge is inconsistent with node metadata.");
            }

            result.signed_edge_flow[edge_id] += side_direction * demand.value;

            /*
             * Commodity load is unsigned. Regardless of which direction
             * the commodity traverses the edge, it contributes value d.
             */
            result.edge_load[edge_id] += demand.value;

            steps.push_back(
                TreePathStep{
                    .edge_id = edge_id,
                    .direction = side_direction
                }
            );

            node = parent;
            return {};
        };

    /*
     * First equalize depths.
     */
    while (tree_.nodes[source_node].depth > tree_.nodes[target_node].depth) {
        auto moved = moveUp(source_node,+1,source_steps);
        if (!moved) {
            return getError(moved);
        }
    }

    while (tree_.nodes[target_node].depth > tree_.nodes[source_node].depth) {
        auto moved = moveUp(target_node,-1,target_steps);

        if (!moved) {
            return getError(moved);
        }
    }

    /*
     * Move both nodes upward until their first common ancestor.
     */
    while (source_node != target_node) {
        auto source_moved = moveUp(source_node,+1,source_steps);

        if (!source_moved) {
            return getError(source_moved);
        }

        auto target_moved = moveUp(target_node,-1,target_steps);

        if (!target_moved) {
            return getError(target_moved);
        }
    }

    const int lca = source_node;

    if (record_path) {
        TreeCommodityPath path;
        path.source = demand.source;
        path.target = demand.target;
        path.value = demand.value;
        path.lca = lca;

        path.steps.reserve(
            source_steps.size() +
            target_steps.size()
        );

        /*
         * Source steps are already ordered:
         *
         *     source leaf -> LCA.
         */
        path.steps.insert(
            path.steps.end(),
            source_steps.begin(),
            source_steps.end()
        );

        /*
         * target_steps were discovered while climbing:
         *
         *     target leaf -> LCA.
         *
         * The actual commodity direction is:
         *
         *     LCA -> target leaf.
         *
         * Therefore reverse their order. Each step already has
         * direction -1 relative to child -> parent storage.
         */
        for (auto it = target_steps.rbegin();it != target_steps.rend();++it) {
            path.steps.push_back(*it);
        }

        result.commodity_paths.push_back(std::move(path));
    }

    return {};
}

Result<void> TreeFlowRouter::finalizeAndValidate(TreeFlowResult& result,const std::vector<double>& expected_divergence) const {
    if (result.signed_edge_flow.size() != tree_.edges.size() ||result.edge_load.size() != tree_.edges.size() ||result.edge_congestion.size() != tree_.edges.size()) {
        return makeErrorMessage(ErrorCode::SolverFailed,"Tree flow result arrays have inconsistent sizes.");
    }

    std::vector<double> actual_divergence(tree_.nodes.size(),0.0);

    result.max_congestion = 0.0;
    result.bottleneck_edge = -1;

    for (const auto& edge : tree_.edges) {
        if (!validNode(tree_, edge.child) || !validNode(tree_, edge.parent)) {
            return makeErrorMessage(ErrorCode::InvalidGraph,"Tree contains an edge with invalid endpoints.");
        }

        if (!std::isfinite(edge.capacity) ||
            edge.capacity <= 0.0) {
            return makeErrorMessage(ErrorCode::InvalidGraph,"Tree contains a non-positive edge capacity.");
        }

        const double signed_flow = result.signed_edge_flow[edge.id];
        const double load = result.edge_load[edge.id];

        if (!std::isfinite(signed_flow) || !std::isfinite(load) || load < -kFlowEpsilon) {
            return makeErrorMessage(ErrorCode::SolverFailed,"Tree routing produced an invalid edge flow.");
        }

        /*
         * Net signed flow can never exceed total commodity load.
         */
        if (std::abs(signed_flow) > load + kFlowEpsilon) {
            return makeErrorMessage(ErrorCode::SolverFailed,"Signed tree flow exceeds total edge load.");
        }

        /*
         * Positive signed flow is child -> parent.
         */
        actual_divergence[edge.child] += signed_flow;
        actual_divergence[edge.parent] -= signed_flow;

        const double congestion = load / edge.capacity;

        result.edge_congestion[edge.id] = congestion;

        if (congestion > result.max_congestion) {
            result.max_congestion = congestion;
            result.bottleneck_edge = edge.id;
        }
    }

    double max_error = 0.0;

    for (int node = 0;node < static_cast<int>(tree_.nodes.size());++node) {
        const double error = std::abs(actual_divergence[node] -expected_divergence[node]);

        max_error = std::max(max_error, error);
    }

    result.max_conservation_error = max_error;

    if (!std::isfinite(max_error) ||max_error > kFlowEpsilon) {
        return makeErrorMessage(ErrorCode::SolverFailed,"Tree-flow conservation failed with error " +std::to_string(max_error) + ".");
    }

    return {};
}