//
// Created by Mert Biyikli on 13.07.26.
//

#include "algorithms/semi_oblivious/expander_hierarchy/tree_sparsifier.h"

#include <iostream>
#include <queue>

template<typename Callback>
Result<void>
forEachUndirectedEdge(
    const optimized::Graph<EdgeData>& graph,
    Callback&& callback
) {
    /*
     * optimized::Graph stores both orientations of each undirected edge.
     * Keep only the canonical u < v representative.
     */
    for (int edge_id = 0;
         edge_id < graph.getNumDirectedEdges();
         ++edge_id) {
        const auto [u, v] =
            graph.getEdgeEndpoints(edge_id);

        if (u == v) {
            continue;
        }

        if (u > v) {
            continue;
        }

        if (u < 0 || v < 0 ||
            u >= graph.getNumNodes() ||
            v >= graph.getNumNodes()) {
            return makeErrorMessage(
                ErrorCode::InvalidGraph,
                "Original graph contains an invalid edge endpoint."
            );
        }

        const double capacity =
            graph.edgeData(edge_id).capacity;

        if (!std::isfinite(capacity) || capacity < 0.0) {
            return makeErrorMessage(
                ErrorCode::InvalidGraph,
                "Original graph contains a non-finite or negative capacity."
            );
        }

        callback(edge_id, u, v, capacity);
    }

    return {};
}

Result<std::vector<double>>
computeWeightedDegrees(const optimized::Graph<EdgeData>& graph) {
    std::vector<double> weighted_degree(
        graph.getNumNodes(),
        0.0
    );

    auto result = forEachUndirectedEdge(
        graph,
        [&](const int,
            const int u,
            const int v,
            const double capacity) {
            weighted_degree[u] += capacity;
            weighted_degree[v] += capacity;
        }
    );

    if (!result) {
        return std::unexpected(result.error());
    }

    return weighted_degree;
}

Result<std::unordered_map<int, double>> computeClusterBoundaryCapacities(const optimized::Graph<EdgeData>& graph,const HierarchyResult& hierarchy) {
    std::unordered_map<int, double> boundary_capacity;

    /*
     * At every hierarchy level, clusters partition the original
     * vertices. For each original edge crossing two clusters, its
     * capacity contributes once to the boundary of each cluster.
     */
    for (const auto& level : hierarchy.levels) {
        std::vector<int> vertex_owner(graph.getNumNodes(),-1);

        for (const auto& cluster : level.clusters) {
            if (cluster.id < 0) {
                return makeErrorMessage(ErrorCode::InvalidGraph,"Hierarchy contains an invalid cluster ID.");
            }

            const auto [_, inserted] = boundary_capacity.try_emplace(cluster.id,0.0);

            if (!inserted) {
                return makeErrorMessage(ErrorCode::InvalidGraph,"Hierarchy contains duplicate cluster ID " +std::to_string(cluster.id) + ".");
            }

            for (const int v : cluster.original_vertices) {
                if (v < 0 ||  v >= graph.getNumNodes()) {
                    return makeErrorMessage(ErrorCode::InvalidGraph,"Hierarchy cluster contains an invalid vertex.");
                }

                if (vertex_owner[v] != -1) {
                    return makeErrorMessage(ErrorCode::InvalidGraph,"Hierarchy level is not a vertex partition.");
                }

                vertex_owner[v] = cluster.id;
            }
        }

        for (int v = 0; v < graph.getNumNodes(); ++v) {
            if (vertex_owner[v] < 0) {
                return makeErrorMessage(ErrorCode::InvalidGraph,"Vertex " + std::to_string(v) +" is missing from hierarchy level " +std::to_string(level.level) + ".");
            }
        }

        auto result = forEachUndirectedEdge(graph,[&](const int,const int u,const int v,const double capacity) {
                const int cluster_u = vertex_owner[u];
                const int cluster_v = vertex_owner[v];

                if (cluster_u == cluster_v) {
                    return;
                }

                boundary_capacity[cluster_u] += capacity;
                boundary_capacity[cluster_v] += capacity;
            }
        );

        if (!result) {
            return getError(result);
        }
    }

    return boundary_capacity;
}

int addClusterNode(TreeSparsifier& tree,const HierarchyCluster& cluster) {
    const int node_id = static_cast<int>(tree.nodes.size());

    TreeSparsifierNode node;
    node.id = node_id;
    node.kind = TreeNodeKind::Cluster;
    node.cluster_id = cluster.id;
    node.hierarchy_level = cluster.level;

    tree.nodes.push_back(std::move(node));
    tree.cluster_to_node.emplace(cluster.id, node_id);

    return node_id;
}

int addOriginalVertexNode(TreeSparsifier& tree,const int original_vertex) {
    const int node_id = static_cast<int>(tree.nodes.size());

    TreeSparsifierNode node;
    node.id = node_id;
    node.kind = TreeNodeKind::OriginalVertex;
    node.original_vertex = original_vertex;
    node.hierarchy_level = -1;

    tree.nodes.push_back(std::move(node));

    return node_id;
}

Result<void> addTreeEdge(TreeSparsifier& tree,const int child,const int parent,const double capacity) {
    if (child < 0 || parent < 0 || child >= static_cast<int>(tree.nodes.size()) || parent >= static_cast<int>(tree.nodes.size())) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Cannot add tree edge with invalid endpoint.");
    }

    if (child == parent) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Tree edge cannot be a self-loop.");
    }

    if (!std::isfinite(capacity) || capacity <= 0.0) {
        return makeErrorMessage( ErrorCode::InvalidGraph,"Tree edge must have finite positive capacity.");
    }

    if (tree.nodes[child].parent.has_value()) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Tree node already has a parent.");
    }

    const int edge_id = static_cast<int>(tree.edges.size());

    tree.edges.push_back(TreeSparsifierEdge{
            .id = edge_id,
            .child = child,
            .parent = parent,
            .capacity = capacity
        }
    );

    tree.nodes[child].parent = parent;
    tree.nodes[child].parent_edge = edge_id;
    tree.nodes[parent].children.push_back(child);

    return {};
}

Result<void> assignDepthsAndValidateTree(TreeSparsifier& tree) {
    if (tree.root < 0 ||
        tree.root >= static_cast<int>(tree.nodes.size())) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Tree sparsifier has no valid root.");
    }

    if (tree.nodes[tree.root].parent.has_value()) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Tree root cannot have a parent.");
    }

    if (tree.edges.size() + 1 != tree.nodes.size()) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Tree sparsifier does not contain |V|-1 edges.");
    }

    std::vector<char> visited(tree.nodes.size(), false);
    std::queue<int> queue;

    tree.nodes[tree.root].depth = 0;
    visited[tree.root] = true;
    queue.push(tree.root);

    std::size_t visited_count = 0;

    while (!queue.empty()) {
        const int current = queue.front();
        queue.pop();

        ++visited_count;

        for (const int child :
             tree.nodes[current].children) {
            if (child < 0 ||
                child >= static_cast<int>(tree.nodes.size())) {
                return makeErrorMessage(
                    ErrorCode::InvalidGraph,
                    "Tree contains an invalid child reference."
                );
            }

            if (visited[child]) {
                return makeErrorMessage(
                    ErrorCode::InvalidGraph,
                    "Cycle detected in tree sparsifier."
                );
            }

            visited[child] = true;
            tree.nodes[child].depth =
                tree.nodes[current].depth + 1;

            queue.push(child);
        }
    }

    if (visited_count != tree.nodes.size()) {
        return makeErrorMessage(
            ErrorCode::InvalidGraph,
            "Tree sparsifier is disconnected."
        );
    }

    return {};
}


Result<TreeSparsifier> TreeSparsifierBuilder::build(const optimized::Graph<EdgeData> &graph, const HierarchyResult &hierarchy) const {
    auto degrees = computeWeightedDegrees(graph);
    auto boundary = computeClusterBoundaryCapacities(graph, hierarchy);

    if (!degrees) {
        return getError(degrees);
    }
    if (!boundary) {
        return getError(boundary);
    }
    auto weighted_degree = degrees.value();
    auto boundary_capacity = boundary.value();

    TreeSparsifier tree;
    tree.vertex_to_leaf.assign(graph.getNumNodes(), -1);

    /*
     * First create one tree node for every hierarchy cluster.
     */
    for (const auto& level : hierarchy.levels) {
        for (const auto& cluster : level.clusters) {
            if (tree.cluster_to_node.contains(cluster.id)) {
                return makeErrorMessage(ErrorCode::InvalidGraph,"Duplicate cluster ID while building tree.");
            }
            addClusterNode(tree, cluster);
        }
    }

    const auto root_it = tree.cluster_to_node.find(*hierarchy.root);

    if (root_it == tree.cluster_to_node.end()) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Hierarchy root cluster is missing.");
    }

    tree.root = root_it->second;

    /*
     * Connect every non-root cluster to its parent.
     *
     * Capacity:
     *
     *     c_T(S, parent(S)) = c_G(delta(S))
     */
    for (const auto& level : hierarchy.levels) {
        for (const auto& cluster : level.clusters) {
            if (cluster.id == *hierarchy.root) {
                if (cluster.parent.has_value()) {
                    return makeErrorMessage(ErrorCode::InvalidGraph,"Root cluster unexpectedly has a parent.");
                }
                continue;
            }

            if (!cluster.parent.has_value()) {
                return makeErrorMessage(ErrorCode::InvalidGraph,"Non-root cluster " +std::to_string(cluster.id) +" has no parent.");
            }

            const auto child_it = tree.cluster_to_node.find(cluster.id);
            const auto parent_it = tree.cluster_to_node.find(*cluster.parent);

            if (child_it == tree.cluster_to_node.end() ||
                parent_it == tree.cluster_to_node.end()) {
                return makeErrorMessage(ErrorCode::InvalidGraph,"Cluster-to-tree-node mapping is incomplete.");
            }

            const auto capacity_it = boundary_capacity.find(cluster.id);

            if (capacity_it == boundary_capacity.end()) {
                return makeErrorMessage(ErrorCode::InvalidGraph,"Boundary capacity is missing for cluster " +std::to_string(cluster.id) + ".");
            }

            auto added = addTreeEdge(
                tree,
                child_it->second,
                parent_it->second,
                capacity_it->second
            );

            if (!added) {
                return getError(added);
            }
        }
    }

    /*
     * A tree flow sparsifier must have the original vertices as leaves.
     *
     * Attach every original vertex v below its finest cluster using:
     *
     *     c_T(v, leafCluster(v)) = c_G(delta({v}))
     *                              = weighted degree of v.
     */
    for (int v = 0;v < graph.getNumNodes();++v) {
        const int finest_cluster = hierarchy.vertex_to_leaf[v];
        const auto cluster_it = tree.cluster_to_node.find(finest_cluster);

        if (cluster_it == tree.cluster_to_node.end()) {
            return makeErrorMessage(ErrorCode::InvalidGraph,"Finest cluster for vertex " +std::to_string(v) +" does not exist.");
        }

        const double degree = weighted_degree[v];

        if (!std::isfinite(degree) || degree <= 0.0) {
            return makeErrorMessage(ErrorCode::InvalidGraph,"Original vertex " +std::to_string(v) +" has non-positive weighted degree.");
        }

        const int leaf = addOriginalVertexNode(tree, v);

        tree.vertex_to_leaf[v] = leaf;

        auto added = addTreeEdge(
            tree,
            leaf,
            cluster_it->second,
            degree
        );

        if (!added) {
            return getError(added);
        }
    }

    auto validation =
        assignDepthsAndValidateTree(tree);

    if (!validation) {
        return getError(validation);
    }

    return tree;
}

void TreeSparsifier::print() const {
    std::cout << "TreeSparsifier:" << std::endl;
    std::cout << "Root: " << root << std::endl;
    std::cout << "Nodes: " << nodes.size() << std::endl;
    std::cout << "Edges: " << edges.size() << std::endl;

    for (const auto& node : nodes) {
        std::cout << "Node ID: " << node.id
                  << ", Kind: " << static_cast<int>(node.kind)
                  << ", Depth: " << node.depth
                  << ", Parent: "
                  << (node.parent.has_value() ? std::to_string(*node.parent) : "None")
                  << ", Children: [";
        for (const auto& child : node.children) {
            std::cout << child << " ";
        }
        std::cout << "]" << std::endl;
    }

    for (const auto& edge : edges) {
        std::cout << "Edge ID: " << edge.id
                  << ", Child: " << edge.child
                  << ", Parent: " << edge.parent
                  << ", Capacity: " << edge.capacity
                  << std::endl;
    }
}
