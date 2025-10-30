//
// Created by Mert Biyikli on 25.10.25.
//

#ifndef OBLIVIOUSROUTING_CKR_TREE_DECOMPOSER_H
#define OBLIVIOUSROUTING_CKR_TREE_DECOMPOSER_H

#include "ckr_partition.h"
#include <vector>

struct TreeNode; // Forward declaration

void print_tree(TreeNode* node, int indent = 0);

struct TreeNode {
    int id;
    std::vector<int> members; // original graph nodes
    std::vector<TreeNode*> children;
    TreeNode* parent;
    double radius = 0.0;
};


class TreeDecomposer {
public:
    TreeNode* decompose(Graph& G, double delta, const std::vector<int>& global_node_ids, int level = 0) {
        const int n = G.getNumNodes();
        if (n == 1) {
            TreeNode* leaf = new TreeNode{level, {0}, {}};
            return leaf;
        }


        CKRPartion ckr;
        ckr.init(G);
        std::vector<int> X(n);
        std::iota(X.begin(), X.end(),0);
        auto clusters = ckr.computePartition(X, delta);


        std::unordered_map<int, std::vector<int>> cluster_to_nodes;
        for (int i = 0; i < clusters.size(); ++i) {
            cluster_to_nodes[clusters[i]].push_back(i);
        }


        TreeNode* root = new TreeNode{level};

        for (auto& [cid, nodes] : cluster_to_nodes) {
            std::unordered_map<int, int> index_map;
            std::vector<int> local_to_global;
            for (int i = 0; i < nodes.size(); ++i) {
                index_map[nodes[i]] = i;
                local_to_global.push_back(global_node_ids[nodes[i]]);
            }

            Graph subgraph(nodes.size());
            for (int u : nodes) {
                for (int v : G.neighbors(u)) {
                    if (index_map.count(v)) {
                        int i = index_map[u], j = index_map[v];
                        if (i > j) continue;
                        subgraph.addEdge(i, j, G.getEdgeCapacity(u, v));
                    }
                }
            }

            TreeNode* child = decompose(subgraph, delta * 0.5, local_to_global, level + 1);
            child->members = local_to_global;  // These are now true original graph node IDs
            root->children.push_back(child);
        }



        return root;
    }
};


#endif //OBLIVIOUSROUTING_CKR_TREE_DECOMPOSER_H