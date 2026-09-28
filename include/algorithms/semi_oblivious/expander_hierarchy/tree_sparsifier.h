//
// Created by Mert Biyikli on 13.07.26.
//

#ifndef OBLIVIOUSROUTING_TREE_SPARSIFIER_H
#define OBLIVIOUSROUTING_TREE_SPARSIFIER_H


#include "core/errors.h"
#include "data_structures/graph/Igraph.h"
#include "data_structures/graph/graph.h"
#include <vector>

#include "hierarchy_results.h"

enum class TreeNodeKind {
    OriginalVertex,
    Cluster
};

struct TreeSparsifierNode {
    int id = -1;
    TreeNodeKind kind = TreeNodeKind::Cluster;

    // Set only for cluster nodes.
    std::optional<int> cluster_id;

    // Set only for original-vertex leaves.
    std::optional<int> original_vertex;

    // Hierarchy level for cluster nodes.
    // Original vertices use -1.
    int hierarchy_level = -1;

    std::optional<int> parent;
    std::optional<int> parent_edge;

    std::vector<int> children;

    // Root has depth 0.
    int depth = -1;

    [[nodiscard]]
    bool isLeaf() const noexcept {
        return kind == TreeNodeKind::OriginalVertex;
    }
};

struct TreeSparsifierEdge {
    int id = -1;

    // Store the orientation child -> parent for convenience.
    int child = -1;
    int parent = -1;

    double capacity = 0.0;
};


/**
 * The Tree Sparsifier is build by an expander hierarchy. For each non-leaf cluster, we add
 * a parent node that is connected to each node in the cluster. The edge capacity
 * is cap(child, parent) = \sum_{u \in N(v)} c(v, u) where v is the original vertex and N(v) is the set of neighbors of v in the original graph.
 */
struct TreeSparsifier {
    int root = -1;

    std::vector<TreeSparsifierNode> nodes;
    std::vector<TreeSparsifierEdge> edges;

    // Original graph vertex -> corresponding tree leaf.
    std::vector<int> vertex_to_leaf;

    // Hierarchy cluster -> corresponding tree node.
    std::unordered_map<int, int> cluster_to_node;


    bool empty() const noexcept {
        return nodes.empty();
    }

    const TreeSparsifierNode*
    findNode(const int id) const noexcept {
        if (id < 0 || id >= static_cast<int>(nodes.size())) {
            return nullptr;
        }

        return &nodes[id];
    }

    const TreeSparsifierEdge* findEdge(const int id) const noexcept {
        if (id < 0 || id >= static_cast<int>(edges.size())) {
            return nullptr;
        }

        return &edges[id];
    }

    void print() const;
};

class TreeSparsifierBuilder {
public:
    Result<TreeSparsifier> build(const optimized::Graph<EdgeData>& graph,const HierarchyResult& hierarchy) const;
};
#endif //OBLIVIOUSROUTING_TREE_SPARSIFIER_H