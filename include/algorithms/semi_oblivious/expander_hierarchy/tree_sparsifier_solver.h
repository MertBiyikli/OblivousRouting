//
// Created by Mert Biyikli on 14.07.26.
//

#ifndef OBLIVIOUSROUTING_TREE_SPARSIFIER_SOLVER_H
#define OBLIVIOUSROUTING_TREE_SPARSIFIER_SOLVER_H
#include "tree_flow_electrical_embedder.h"
#include "tree_flow_router.h"
#include "tree_sparsifier.h"
#include "algorithms/oblivious/oblivious_solver.h"
#include "preprocessing/hierarchy_preprocessor.h"
#include "utils/time_tracking.h"

class ElectrifiedExpanderHierarchySolver : public ILinearObliviousSolverBase {
    public:
    std::unique_ptr<TreeSparsifier> tree_;
    std::unique_ptr<HierarchyResult> hierarchy_;
    ExpanderMetrics metrics_;

    ElectrifiedExpanderHierarchySolver(optimized::Graph<EdgeData>& g, int root)
        : ILinearObliviousSolverBase(g, root){
    }

    Result<void> computeBasisFlows(LinearRoutingTable& table) override {
        auto res = init();
        if (!res) {
            return getError(res);
        }

        res = run(table);
        if (!res) {
            return getError(res);
        }

        return {};
    }

    Result<void> init(std::vector<double> edge_weights = {}) {
        auto start = timeNow();
        XCutHierarchyPreprocessor preprocessor;
        if (edge_weights.empty()) {
            for (int e = 0; e<graph.getNumDirectedEdges(); e++) {
                graph.edgeData(e).weight = graph.edgeData(e).capacity;
            }
        }else {
            for (int e = 0; e<graph.getNumDirectedEdges(); e++) {
                graph.edgeData(e).weight = edge_weights[e];
            }
        }

        auto hierarchy_result = preprocessor.build(graph);

        metrics_.hierarchy_runtime_microseconds = duration(timeNow() - start);
        if(!hierarchy_result) {
            return getError(hierarchy_result);
        }
        hierarchy_ = std::make_unique<HierarchyResult>(hierarchy_result.value());

        start = timeNow();
        auto tmp = collectHierarchyMetrics();
        if (!tmp) {
            return getError(tmp);
        }

        TreeSparsifierBuilder tree_builder;

        auto tree_result = tree_builder.build(graph, *hierarchy_);

        metrics_.tree_runtime_microseconds = duration(timeNow() - start);
        if (!tree_result) {
            return getError(tree_result);
        }

        tree_ = std::make_unique<TreeSparsifier>(tree_result.value());
        tmp = collectTreeMetrics();
        if (!tmp) {
            return getError(tmp);
        }
        return {};
    }

    Result<void> run(LinearRoutingTable& table) {
        if (!tree_  || !hierarchy_) {
            return makeErrorMessage(ErrorCode::InvalidGraph,"Tree sparsifier or hierarchy is not initialized.");
        }

        if (tree_->empty() || hierarchy_->levels.empty()) {
            return makeErrorMessage(ErrorCode::InvalidGraph,"Cannot route on an empty tree sparsifier.");
        }

        auto start= timeNow();
        TreeFlowRouter tree_router(*tree_);
        for (int target = 0; target < graph.getNumNodes(); ++target) {
            if (target == root) {
                continue;
            }

            auto tree_flow_result = tree_router.routePair(
                target,
                root,
                1.0 // -> should be an oblivious routing scheme
            );


            if (!tree_flow_result) {
                return getError(tree_flow_result);
            }

            const TreeFlowResult& tree_flow = *tree_flow_result;

            TreeFlowElectricalEmbedder embedder(graph,*hierarchy_,*tree_);

            auto embedding_result = embedder.embed(tree_flow);

            if (!embedding_result) {
                return getError(embedding_result);
            }

            const auto& embedding = *embedding_result;

            metrics_.total_electrical_solves  += embedding.electrical_solves;
            metrics_.basis_flows += 1;
            metrics_.max_basis_embedding_congestion = std::max(metrics_.max_basis_embedding_congestion,embedding.max_congestion);
            metrics_.max_conservation_error = std::max(metrics_.max_conservation_error,embedding.max_conservation_error);

            for (int e = 0; e < graph.getNumDirectedEdges();++e) {
                const double flow = embedding.signed_edge_flow[e];

                if (std::abs(flow) <= 1e-10) {
                    continue;
                }

                // get sign of the flow
                if (flow < 0) {
                    int anti_e = graph.reverse(e).id;
                    table.addFlow(anti_e, target, std::abs(flow));
                }else {
                    table.addFlow(e, target, flow);
                }
            }
        }
        metrics_.basis_flow_runtime_microseconds = duration(timeNow() - start);
        return {};
    }

    const ExpanderMetrics& getMetrics() const noexcept {
        return metrics_;
    }
    Result<void> collectHierarchyMetrics() {
        if (!hierarchy_) {
            return makeErrorMessage(ErrorCode::InvalidArgument,"Hierarchy is not initialized.");
        }
        metrics_.hierarchy_levels =
            hierarchy_->levels.size();

        metrics_.clusters_per_level.clear();
        metrics_.clusters_per_level.reserve(
            hierarchy_->levels.size()
        );

        std::size_t total_vertex_memberships = 0;

        for (const auto& level : hierarchy_->levels) {
            metrics_.clusters_per_level.push_back(
                level.clusters.size()
            );

            metrics_.hierarchy_clusters +=
                level.clusters.size();

            for (const auto& cluster : level.clusters) {
                const std::size_t cluster_size =
                    cluster.original_vertices.size();

                total_vertex_memberships +=
                    cluster_size;

                metrics_.max_cluster_vertices =
                    std::max(
                        metrics_.max_cluster_vertices,
                        cluster_size
                    );
            }
        }

        if (metrics_.hierarchy_clusters > 0) {
            metrics_.average_cluster_vertices =
                static_cast<double>(
                    total_vertex_memberships
                ) /
                static_cast<double>(
                    metrics_.hierarchy_clusters
                );
        }
        return {};
    }

    Result<void> collectTreeMetrics() {
        if (!tree_) {
            return makeErrorMessage(ErrorCode::InvalidArgument,"Tree sparsifier is not initialized.");
        }
        metrics_.tree_nodes =
            tree_->nodes.size();

        metrics_.tree_edges =
            tree_->edges.size();

        for (const auto& node : tree_->nodes) {
            metrics_.tree_depth =
                std::max(
                    metrics_.tree_depth,
                    node.depth
                );
        }
        return {};
    }
};



#endif //OBLIVIOUSROUTING_TREE_SPARSIFIER_SOLVER_H