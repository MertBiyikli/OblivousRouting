//
// Created by Mert Biyikli on 25.10.25.
//

#include "ckr_tree_decomposer.h"



void print_tree(const std::shared_ptr<ITreeNode>& node, int indent) {
    for (int i = 0; i < indent; ++i) std::cout << " ";
    std::cout << "Center " << node->center << " | Nodes: ";
    for (int u : node->getMembers()) std::cout << u << " ";
    std::cout << std::endl;
    for (auto child : node->getChildren()) {
        print_tree(child, indent + 1);
    }
}