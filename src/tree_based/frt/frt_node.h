//
// Created by Mert Biyikli on 15.12.25.
//

#ifndef OBLIVIOUSROUTING_FRT_NODE_H
#define OBLIVIOUSROUTING_FRT_NODE_H

#include "../raecke_oracle_iteration.h"

class EfficientFRTTreeNode : public ITreeNode {
public:
    int id = -1;
    std::vector<int> members; // original graph nodes
    std::weak_ptr<EfficientFRTTreeNode> parent;
    std::vector<std::shared_ptr<EfficientFRTTreeNode>> children;

    virtual std::shared_ptr<ITreeNode> getParent() const override {
        return parent.lock();
    }

    std::vector<std::shared_ptr<ITreeNode>> getChildren() override {
        std::vector<std::shared_ptr<ITreeNode>> child_ptrs;
        for (std::shared_ptr<EfficientFRTTreeNode>& child : children) {
            child_ptrs.push_back(child);
        }
        return child_ptrs;
    }

    const std::vector<int> &getMembers() const override {
        return members;
    }

};

#endif //OBLIVIOUSROUTING_FRT_NODE_H