//
// Created by Mert Biyikli on 13.07.26.
//

#include "algorithms/semi_oblivious/expander_hierarchy/tree_flow_electrical_embedder.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>



constexpr double kEmbeddingEpsilon = 1e-8;

template<typename Callback>
Result<void> forEachUndirectedEdge(const IGraph& graph, Callback&& callback) {
    for (int e = 0;e < graph.getNumDirectedEdges(); ++e) {
        const auto [u, v] = graph.getEdgeEndpoints(e);

        // Retain one canonical orientation.
        if (u >= v) {
            continue;
        }
        const double capacity = graph.getEdgeCapacity(e);

        if (!std::isfinite(capacity) ||capacity <= 0.0) {
            return makeErrorMessage(ErrorCode::InvalidGraph,"Original graph contains an invalid capacity.");
        }

        callback(e, u, v, capacity);
    }

    return {};
}

double l1Norm(const std::vector<double>& vector) {
    double result = 0.0;

    for (const double value : vector) {
        result += std::abs(value);
    }

    return result;
}

Result<std::unordered_map<int,ClusterBoundaryProfile>> TreeFlowElectricalEmbedder::buildBoundaryProfiles() const {
    std::unordered_map<int,ClusterBoundaryProfile> profiles;

    /*
     * Temporary sparse maps:
     *
     * cluster -> (original vertex -> incident boundary capacity)
     */
    std::unordered_map<int,std::unordered_map<int, double>> per_vertex_capacity;

    for (const auto& level : hierarchy_.levels) {
        std::vector<int> owner(graph_.getNumNodes(),-1);

        for (const auto& cluster : level.clusters) {
            profiles.emplace(cluster.id,ClusterBoundaryProfile{.cluster_id = cluster.id});

            for (const int v : cluster.original_vertices) {
                if (v < 0 || v >= graph_.getNumNodes()) {
                    return makeErrorMessage(ErrorCode::InvalidGraph,"Hierarchy contains an invalid vertex.");
                }

                if (owner[v] != -1) {
                    return makeErrorMessage(ErrorCode::InvalidGraph, "Hierarchy level is not a partition.");
                }

                owner[v] = cluster.id;
            }
        }

        for (int v = 0;v < graph_.getNumNodes();++v) {
            if (owner[v] < 0) {
                return makeErrorMessage(ErrorCode::InvalidGraph,"Hierarchy level does not cover every vertex.");
            }
        }

        auto processed = forEachUndirectedEdge(graph_,[&](const int,const int u,const int v,const double capacity) {
                const int cluster_u = owner[u];
                const int cluster_v = owner[v];

                if (cluster_u == cluster_v) {
                    return;
                }

                profiles[cluster_u].total_capacity +=capacity;
                profiles[cluster_v].total_capacity += capacity;

                per_vertex_capacity[cluster_u][u] += capacity;
                per_vertex_capacity[cluster_v][v] +=capacity;
            }
        );

        if (!processed) {
            return getError(processed);
        }
    }

    for (auto& [cluster_id, profile] : profiles) {
        const auto values_it = per_vertex_capacity.find(cluster_id);

        if (values_it == per_vertex_capacity.end()) {
            continue;
        }

        profile.vertex_capacity.reserve(
            values_it->second.size()
        );

        for (const auto& [vertex, capacity] : values_it->second) {
            profile.vertex_capacity.emplace_back(vertex,capacity);
        }

        std::sort(profile.vertex_capacity.begin(),profile.vertex_capacity.end(),[](const auto& lhs, const auto& rhs) {return lhs.first < rhs.first;});
    }

    /*
     * Check that profile capacities match the capacities used in the
     * tree sparsifier.
     */
    for (const auto& level : hierarchy_.levels) {
        for (const auto& cluster : level.clusters) {
            const auto tree_node_it = tree_.cluster_to_node.find(cluster.id);

            if (tree_node_it == tree_.cluster_to_node.end()) {
                return makeErrorMessage(ErrorCode::InvalidGraph,"Tree is missing a hierarchy cluster.");
            }

            const int tree_node_id = tree_node_it->second;
            const auto& tree_node = tree_.nodes[tree_node_id];

            if (!tree_node.parent_edge.has_value()) {
                continue;
            }

            const int tree_e = *tree_node.parent_edge;
            const double tree_capacity = tree_.edges[tree_e].capacity;
            const double boundary_capacity = profiles.at(cluster.id).total_capacity;
            const double tolerance = kEmbeddingEpsilon * std::max(1.0, tree_capacity);

            if (std::abs(tree_capacity -boundary_capacity) > tolerance) {
                return makeErrorMessage(ErrorCode::InvalidGraph,"Tree edge capacity does not match boundary capacity for cluster " +
                        std::to_string(cluster.id) +". Tree=" +std::to_string(tree_capacity) +", boundary=" +std::to_string(boundary_capacity));
            }
        }
    }

    return profiles;
}

Result<double> TreeFlowElectricalEmbedder::clusterTreeFlow(const int cluster,const TreeFlowResult& tree_flow) const {
    const auto node_it = tree_.cluster_to_node.find(cluster);

    if (node_it == tree_.cluster_to_node.end()) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Cluster has no corresponding tree node.");
    }

    const TreeSparsifierNode& node = tree_.nodes[node_it->second];

    if (!node.parent_edge.has_value()) {
        // Root has no outgoing tree edge.
        return 0.0;
    }

    const int e = *node.parent_edge;

    if (e < 0 || e >= static_cast<int>(tree_flow.signed_edge_flow.size())) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Tree flow is missing a cluster edge.");
    }

    return tree_flow.signed_edge_flow[e];
}

Result<double> TreeFlowElectricalEmbedder::vertexTreeFlow(const int vertex,const TreeFlowResult& tree_flow) const {
    if (vertex < 0 || vertex >=static_cast<int>(tree_.vertex_to_leaf.size())) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Original vertex has no tree leaf.");
    }
    const int leaf = tree_.vertex_to_leaf[vertex];

    if (leaf < 0 ||leaf >=static_cast<int>(tree_.nodes.size())) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Tree leaf mapping is invalid.");
    }

    const auto& leaf_node = tree_.nodes[leaf];

    if (!leaf_node.parent_edge.has_value()) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Original vertex leaf has no parent edge.");
    }

    const int e = *leaf_node.parent_edge;

    if (e < 0 ||e >= static_cast<int>(tree_flow.signed_edge_flow.size())) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Tree flow is missing a vertex-leaf edge.");
    }

    return tree_flow.signed_edge_flow[e];
}

Result<ElectricalEmbeddingResult> TreeFlowElectricalEmbedder::embed(const TreeFlowResult& tree_flow) const {
    if (tree_flow.signed_edge_flow.size() != tree_.edges.size()) {
        return makeErrorMessage(ErrorCode::InvalidDemand,"Tree flow size does not match the tree.");
    }

    auto profiles_result = buildBoundaryProfiles();

    if (!profiles_result) {
        return getError(profiles_result);
    }

    const auto& profiles = *profiles_result;

    ElectricalEmbeddingResult result;

    result.signed_edge_flow.assign(graph_.getNumDirectedEdges(),0.0);
    result.edge_usage.assign(graph_.getNumDirectedEdges(),0.0);

    /*
     * Hierarchy levels are ordered finest -> root, so this is the
     * bottom-up order from the proof.
     */
    for (const auto& level : hierarchy_.levels) {
        for (const auto& cluster : level.clusters) {
            std::unordered_map<int, int> global_to_local;

            global_to_local.reserve(cluster.original_vertices.size());

            for (int local = 0;local <static_cast<int>(cluster.original_vertices.size());++local) {
                global_to_local.emplace(cluster.original_vertices[local],local);
            }

            std::vector<double> local_demand(cluster.original_vertices.size(),0.0);

            /*
             * Input contribution from children.
             *
             * Finest hierarchy clusters have original-vertex tree
             * leaves as children.
             */
            if (cluster.level == 0) {
                for (const int vertex : cluster.original_vertices) {
                    auto flow_result = vertexTreeFlow(vertex,tree_flow);

                    if (!flow_result) {
                        return getError(flow_result);
                    }

                    local_demand[global_to_local.at(vertex)] += *flow_result;
                }
            } else {
                for (const int child_id :cluster.children) {
                    const auto child_it = profiles.find(child_id);

                    if (child_it == profiles.end()) {
                        return makeErrorMessage(ErrorCode::InvalidGraph,"Missing child boundary profile.");
                    }

                    auto child_flow_result = clusterTreeFlow(child_id,tree_flow);

                    if (!child_flow_result) {
                        return getError(child_flow_result);
                    }

                    const double child_flow = *child_flow_result;

                    if (std::abs(child_flow) <=kEmbeddingEpsilon) {
                        continue;
                    }

                    const auto& profile = child_it->second;

                    if (profile.total_capacity <= 0.0) {
                        return makeErrorMessage(ErrorCode::InvalidGraph,"Nonzero child flow has zero boundary.");
                    }

                    for (const auto& [vertex, capacity] : profile.vertex_capacity) {
                        const auto local_it = global_to_local.find(vertex);

                        if (local_it == global_to_local.end()) {
                            return makeErrorMessage(ErrorCode::InvalidGraph,"Child boundary vertex is not contained in parent cluster.");
                        }

                        local_demand[local_it->second] +=child_flow * capacity /profile.total_capacity;
                    }
                }
            }

            /*
             * Output contribution to the parent boundary:
             *
             *     -x_H mu_H
             */
            auto outgoing_result = clusterTreeFlow(cluster.id,tree_flow);

            if (!outgoing_result) {
                return getError(outgoing_result);
            }

            const double outgoing_flow =*outgoing_result;

            const auto& own_profile = profiles.at(cluster.id);

            if (std::abs(outgoing_flow) >
                kEmbeddingEpsilon) {
                if (own_profile.total_capacity <= 0.0) {
                    return makeErrorMessage(ErrorCode::InvalidGraph,"Non-root cluster has nonzero outgoing flow but zero boundary capacity."
                    );
                }

                for (const auto& [vertex, capacity] : own_profile.vertex_capacity) {
                    const auto local_it = global_to_local.find(vertex);

                    if (local_it == global_to_local.end()) {
                        return makeErrorMessage(ErrorCode::InvalidGraph,"Cluster boundary vertex is not in cluster.");
                    }

                    local_demand[local_it->second] -=outgoing_flow * capacity / own_profile.total_capacity;
                }
            }

            const double local_sum = std::accumulate(local_demand.begin(),local_demand.end(),0.0);

            if (std::abs(local_sum) > kEmbeddingEpsilon) {
                return makeErrorMessage(ErrorCode::InvalidDemand,"Induced demand is not balanced in cluster " +std::to_string(cluster.id) +". Sum=" +std::to_string(local_sum) );
            }

            ClusterEmbeddingInfo info;
            info.cluster_id = cluster.id;
            info.level = cluster.level;
            info.outgoing_tree_flow = outgoing_flow;
            info.local_demand_l1 = l1Norm(local_demand);

            bool trivial = true;

            for (const double value : local_demand) {
                if (std::abs(value) > kEmbeddingEpsilon) {
                    trivial = false;
                    break;
                }
            }

            if (!trivial) {
                auto local_flow_result = electrical_solver_.routeDemand(graph_,cluster,local_demand);

                if (!local_flow_result) {
                    return getError(local_flow_result);
                }

                const auto& local_flow = *local_flow_result;

                if (local_flow.edge_flow.size() !=cluster.induced_edges.size()) {
                    return makeErrorMessage(ErrorCode::SolverFailed,"Local electrical flow has invalid size.");
                }

                for (int e = 0; e <static_cast<int>(cluster.induced_edges.size());++e) {
                    const int global_edge = cluster.induced_edges[e];

                    if (global_edge < 0 || global_edge >= graph_.getNumDirectedEdges()) {
                        return makeErrorMessage(ErrorCode::InvalidGraph, "Cluster contains invalid edge ID.");
                    }

                    const double flow = local_flow.edge_flow[e];

                    if (!std::isfinite(flow)) {
                        return makeErrorMessage(ErrorCode::SolverFailed,"Electrical embedding produced non-finite flow.");
                    }

                    result.signed_edge_flow[global_edge] += flow;

                    result.edge_usage[global_edge] += std::abs(flow);
                }

                info.electrical_solve_performed = true;
                info.local_max_congestion = local_flow.max_congestion;
                info.local_conservation_error = local_flow.max_conservation_error;

                ++result.electrical_solves;
            }

            result.cluster_info.push_back(info);
        }
    }

    /*
     * Validate the final original-graph flow.
     *
     * The expected divergence at original vertex v is exactly the
     * signed tree flow leaving the tree leaf representing v.
     */
    std::vector<double> actual_divergence(graph_.getNumNodes(),0.0);

    for (int e = 0;e < graph_.getNumDirectedEdges();++e) {
        const double signed_flow = result.signed_edge_flow[e];
        const double usage = result.edge_usage[e];

        if (!std::isfinite(signed_flow) || !std::isfinite(usage)) {
            return makeErrorMessage( ErrorCode::SolverFailed,"Final embedded flow is non-finite.");
        }

        if (std::abs(signed_flow) <=kEmbeddingEpsilon &&usage <= kEmbeddingEpsilon) {
            continue;
        }

        const auto [u, v] = graph_.getEdgeEndpoints(e);

        const double capacity = graph_.getEdgeCapacity(e);

        if (!std::isfinite(capacity) || capacity <= 0.0) {
            return makeErrorMessage( ErrorCode::InvalidGraph,"Used original edge has invalid capacity.");
        }

        actual_divergence[u] += signed_flow;
        actual_divergence[v] -= signed_flow;

        const double usage_congestion = usage / capacity;
        const double net_congestion = std::abs(signed_flow) / capacity;

        if (usage_congestion > result.max_congestion) {
            result.max_congestion = usage_congestion;
            result.bottleneck_edge = e;
        }

        result.max_net_congestion = std::max(result.max_net_congestion,net_congestion);
    }

    for (int vertex = 0;vertex < graph_.getNumNodes();++vertex) {
        auto expected_result =vertexTreeFlow(vertex,tree_flow);

        if (!expected_result) {
            return getError(expected_result);
        }

        const double error = std::abs(actual_divergence[vertex] -*expected_result);

        result.max_conservation_error = std::max(result.max_conservation_error,error);
    }

    if (result.max_conservation_error >kEmbeddingEpsilon) {
        return makeErrorMessage(ErrorCode::SolverFailed,"Final graph embedding violates conservation. Error=" +std::to_string(result.max_conservation_error));
    }

    return result;
}