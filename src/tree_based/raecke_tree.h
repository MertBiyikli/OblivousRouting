//
// Created by Mert Biyikli on 09.12.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_TREE_H
#define OBLIVIOUSROUTING_RAECKE_TREE_H
#include <memory>
// Tree interface, abstract base class
class ITreeNode {
public:
    int center = -1;
    virtual ~ITreeNode() = default;
    virtual std::vector<std::shared_ptr<ITreeNode>> getChildren() = 0;
    virtual std::shared_ptr<ITreeNode> getParent() const = 0;
    virtual const std::vector<int>& getMembers() const = 0;
};


class ITree {
public:
    virtual ~ITree() = default;

    // Root node of the decomposition
    virtual ITreeNode* getRoot() const = 0;
};

#endif //OBLIVIOUSROUTING_RAECKE_TREE_H