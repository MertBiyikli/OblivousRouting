//
// Created by Mert Biyikli on 25.10.25.
//

#include "ckr_tree_decomposer.h"



void print_tree(TreeNode* node, int indent) {
    for (int i = 0; i < indent; ++i) std::cout << " ";
    std::cout << "Level " << node->id << " | Nodes: ";
    for (int u : node->members) std::cout << u << " ";
    std::cout << std::endl;
    for (auto child : node->children) {
        print_tree(child, indent + 1);
    }
}