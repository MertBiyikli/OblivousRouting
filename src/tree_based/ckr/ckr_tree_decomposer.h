//
// Created by Mert Biyikli on 25.10.25.
//

#ifndef OBLIVIOUSROUTING_CKR_TREE_DECOMPOSER_H
#define OBLIVIOUSROUTING_CKR_TREE_DECOMPOSER_H

#include "ckr_partition.h"
#include <vector>
#include <unordered_map>
#include <memory>
#include "../raecke_oracle_iteration.h"
#include "../raecke_tree.h"

void print_tree(std::shared_ptr<ITreeNode> node, int indent = 0);

class TreeNode : public std::enable_shared_from_this<TreeNode>, public ITreeNode {
public:
    int id{};
    std::vector<int> members; // original graph nodes
    std::vector<std::shared_ptr<TreeNode>> children;
    std::weak_ptr<TreeNode> parent;
    double radius = 0.0;

    virtual std::vector<std::shared_ptr<ITreeNode>> getChildren() override{
        std::vector<std::shared_ptr<ITreeNode>> child_ptrs;
        for (const auto& child : children) {
            child_ptrs.push_back(child);
        }
        return child_ptrs;
    }
    virtual std::shared_ptr<ITreeNode> getParent() const override{
        return parent.lock();
    }
    virtual const std::vector<int>& getMembers() const override {
        return members;
    }
};


class TreeDecomposer {
public:
    std::shared_ptr<TreeNode> decompose(GraphADJ& G, double delta, const std::vector<int>& global_node_ids, int level = 0) {
        const int n = G.getNumNodes();
        if (n == 1) {
            auto leaf = std::make_shared<TreeNode>();
            leaf->id = level;
            leaf->members.push_back(global_node_ids[0]);
            leaf->radius = 0.0;
            return leaf;
        }


        CKRPartition ckr;
        ckr.init(G);
        std::vector<int> X(n);
        std::iota(X.begin(), X.end(),0);
        auto clusters = ckr.computePartition(X, delta);


        std::unordered_map<int, std::vector<int>> cluster_to_nodes;
        for (int i = 0; i < clusters.size(); ++i) {
            cluster_to_nodes[clusters[i]].push_back(i);
        }


        std::shared_ptr<TreeNode> root = std::make_shared<TreeNode>();
        root->radius = delta;

        for (auto& [cid, nodes] : cluster_to_nodes) {
            std::unordered_map<int, int> index_map;
            std::vector<int> local_to_global;
            for (int i = 0; i < nodes.size(); ++i) {
                index_map[nodes[i]] = i;
                local_to_global.push_back(global_node_ids[nodes[i]]);
            }

            GraphADJ subgraph(nodes.size());
            for (int u : nodes) {
                for (int v : G.neighbors(u)) {
                    if (index_map.count(v)) {
                        int i = index_map[u], j = index_map[v];
                        if (i > j) continue;

//#ifdef DEBUG
                        double cap = G.getEdgeCapacity(u, v);
                        double dist = G.getEdgeDistance(u, v);
                        std::cout << "Adding edge in subgraph: " << i << " - " << j << " | cap: " << cap << ", dist: " << dist << std::endl;
//#endif
                        subgraph.addEdge(i, j, G.getEdgeCapacity(u, v), G.getEdgeDistance(u, v));
                    }
                }
            }

            std::shared_ptr<TreeNode> child = decompose(subgraph, delta * 0.5, local_to_global, level + 1);
            child->members = local_to_global;  // These are now true original graph node IDs
            child->radius = delta * 0.5;
            root->children.push_back(std::move(child));
        }

        return root;
    }
};


#endif //OBLIVIOUSROUTING_CKR_TREE_DECOMPOSER_H